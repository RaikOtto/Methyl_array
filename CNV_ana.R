# prep 

#annot_match = match( rownames(samples),  rownames(IlluminaHumanMethylationEPICanno.ilm10b2.hg19::Other ) , nomatch = 0)
#genes = IlluminaHumanMethylationEPICanno.ilm10b2.hg19::Other$UCSC_RefGene_Name[annot_match]
#genes = sapply( genes, FUN = function(vec){return(
#  paste( as.character( unique(as.character(unlist(str_split(vec,pattern = ";"))))), collapse = ";", sep ="" )
#)} )
#saveRDS(genes,"~/Koop_Klinghammer/Misc/HGNC_EPIC_name_mapping.RData")

#entrez_gene_set = AnnotationDbi::select(
#  org.Hs.eg.db::org.Hs.eg.db,
#  hgnc_gene_set,
#  c("ENTREZID"),
#  "SYMBOL"
#)
#saveRDS(entrez_gene_set,"~/Koop_Klinghammer/Misc/ENTREZ_EPIC_name_mapping.RData")

#tx = AnnotationDbi::select(
#    TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene,
#    entrez_gene_set$ENTREZID,
#    columns = "TXNAME",
#    keytype = "GENEID"
#)$TXNAME

#tx = GenomicFeatures::transcripts(
#    TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
#)$tx_name
#saveRDS(tx,"~/Koop_Klinghammer/Misc/Tx.RData")

## use z-transformation

workflow <- "C" #B, A, C

## normalize samples and controls sdasd
## together

rgSet = rgSet[,1:2]

normData = minfi::getCN(minfi::preprocessIllumina(rgSet))
colnames(normData)
dim(normData)

switch(
    workflow,
    A = {
        ctrlAll <- normData[,]
        ctrl <- apply(ctrlAll, MARGIN = 1, FUN= median)
    },
    B = {
        ctrlAll <- normData[,]
        ctrlAll[is.infinite(ctrlAll)] <- NA
        ctrlAll <- scale(ctrlAll)
        ctrl <- apply(ctrlAll, 1, "median")
    },
    C = {ctrl = normData[,]}
)

## samples
switch(workflow,
       A = {
         #without z-transformation
         samples <- normData[,]
       },
       B = {
         #with z-transformation
         samples <- normData[,]
         samples[is.infinite(samples)] <- NA
         samples <- scale(samples)
       },
       C = {
         #conumee path, Illumina
         samples <- normData[,]
       },
       {
         stop("Invalid Workflow! Please set workflow to either 'A', 'B' or 'C'.")}
)

candidates <- "genes" #segments, bins, transcripts, genes

candidatesDATA <- NULL
if (candidates == "segments") {
  
    if (workflow == "C") {
      #calculate segments with conumee
      candidatesDATA <- cnAnalysis450k::runConumee(samples, ctrl)
    } else {
      candidatesDATA <-
        cnAnalysis450k::findSegments(samples[, , drop = FALSE], ctrl, ctrlAll)
    }
} else if (candidates == "bins") {
    if (workflow == "C") {
      
      candidatesDATA = cnAnalysis450k::runConumee(samples, ctrl, what = "bins")
      
    } else {
      
      candidatesDATA = cnAnalysis450k::createBinsFast(samples[, 1:3, drop = FALSE], ctrl, ctrlAll, binsize = 500000)
    }
} else if (candidates == "transcripts" || candidates == "genes") {
  
    hgnc_gene_set = readRDS("~/Koop_Klinghammer/Misc/HGNC_EPIC_name_mapping.RData")
    entrez_gene_set = readRDS("~/Koop_Klinghammer/Misc/ENTREZ_EPIC_name_mapping.RData")
    tx = readRDS( "~/Koop_Klinghammer/Misc/Tx.RData" )

    if (workflow == "C") {
      
        candidatesDATA = cnAnalysis450k::runConumee(samples, ctrl, what = "transcripts", tx)
        
    } else {
      
        #candidatesDATA = cnAnalysis450k::getTxValues(samples, ctrl, ctrlAll, tx, output = "diff")
        candidatesDATA <- cnAnalysis450k::getTxValuesFast(samples, ctrl, ctrlAll, tx)
    }
} else {
  stop("Invalid candidate selection! Set to 'segments', 'genes', 
                     'transcripts', 'bins'.")
}

candidatesDATA

candidatesMATRIX <- NULL
if (workflow == "C") {
    
  candidatesMATRIX = cnAnalysis450k::createConumeeMatrix(candidatesDATA)
  
} else {
    if (candidates == "segments") {
    
        candidatesMATRIX = cnAnalysis450k::createSegmentMatrix(candidatesDATA, p.select = 0.05)
    
    } else if (candidates == "bins") {
        
        candidatesMATRIX = cnAnalysis450k::createBinMatrix(candidatesDATA, pval=1)
    
    } else  {
    
        candidatesMATRIX = candidatesDATA$data[which(candidatesDATA$p.val <= 0.05), ]
  }
}


saveRDS(candidatesMATRIX,"~/Koop_Klinghammer/Data/candidatesMatrix.RData")
candidatesMATRIX[1:5,1:5]
dim(candidatesMATRIX)

pheatmap::pheatmap(candidatesMATRIX)

defineCutoffs <- "auto" # man, auto, skip
LOSSEFFECT <- 0.15
GAINEFFECT <- 0.1
PROXIMITY <- c(2, 1)

candidatesCUT <- NULL
candidatesFINAL <- NULL
switch(defineCutoffs,
       man = {
         ##manual thresholds
         #set the thresholds as visually appropriate
         ### CHANGE HERE
         candidatesFINAL <-
           cnAnalysis450k::segmentDataAbs(
             candidatesMATRIX,
             upper = 0.78,
             lower = 1.21,
             ylim = c(0, 2)
           )
       },
       auto = {
         ##auto thresholds
         if (workflow == "C") {
           #conumee
           if (candidates == "segments") {
             candidatesCUT <-
               cnAnalysis450k::findCutoffs(
                 candidatesMATRIX[complete.cases(candidatesMATRIX * 0), , 
                                  drop =FALSE], proximity = PROXIMITY, ignoreNAs = T )
           } else {
             candidatesCUT <-
               cnAnalysis450k::findCutoffs(
                 candidatesMATRIX[complete.cases(candidatesMATRIX * 0), ,
                                  drop = FALSE], ignoreNAs = T)
           }
           candidatesFINAL <-
             cnAnalysis450k::segmentData(
               candidatesMATRIX[complete.cases(candidatesMATRIX * 0), ,
                                drop = FALSE],
               candidatesCUT,
               effectsize = c(LOSSEFFECT, GAINEFFECT), ignoreNAs = T)
         } else {
           if (candidates == "segments") {
             candidatesCUT <-
               cnAnalysis450k::findCutoffs(candidatesMATRIX, proximity = PROXIMITY, ignoreNAs = T)
           } else {
             candidatesCUT <- cnAnalysis450k::findCutoffs(candidatesMATRIX, ignoreNAs = T)
           }
           candidatesFINAL <-
             cnAnalysis450k::segmentData(candidatesMATRIX,
                                         candidatesCUT,
                                         effectsize = c(LOSSEFFECT, GAINEFFECT))
         }
       },
       skip = {
         ##no thresholds
         candidatesFINAL <- candidatesMATRIX
       },
       {
         stop("Invalid cutoff strategy! Please set defineCutoffs to
                'man', 'auto' or 'skip'.")
       })

## should Fishers' p-values be calculated?
calc <- TRUE #TRUE, FALSE
##if so, should they be aggregated?
binF <- FALSE #TRUE, FALSE


###################
fisherVal <- NULL
group <- NULL #change below
if (calc) {
  if (defineCutoffs == "skip") {
    calc <- FALSE
    stop("Please define & apply cutoffs to data!")
  }
  
  ##give group assignment: CHANGE HERE!
  group <- c(rep(1, length(candidatesFINAL[1, ]) - 2), rep(2, 2))
  ##CHECK ORDER !!!
  
  ## calculate fisher values
  fisherVal <-
    cnAnalysis450k::calcFisher(candidatesFINAL, group, bin = binF)
}

## Colors of legends
anno_color <- NULL
if (binF && calc) {
  # set colors for fisher p-val categories
  anno_color <-
    list(
      Fisher = c(
        "<0.001" = "black",
        "<0.01" = "green",
        "<0.05" = "yellow",
        "<0.1" = "lightgray",
        ">0.1" = "white"
      )
    )
}

## Row annotation
anno_row <- NULL
if (candidates == "bin" || candidates == "segments") {
  anno_row <-
    data.frame(id = rownames(candidatesFINAL),
               Chromosome = do.call(rbind, strsplit(rownames(
                 candidatesFINAL
               ), ":"))[, 1])
  anno_row$Chromosome <-
    factor(anno_row$Chromosome, levels = paste("chr", 1:22, sep = ""))
  rownames(anno_row) <- anno_row[, 1]
} else if (candidates == "transcripts" || candidates == "genes") {
  require(Homo.sapiens)
  gen_dat <-
    AnnotationDbi::select(Homo.sapiens,
                          rownames(candidatesFINAL),
                          c("SYMBOL", "TXCHROM"),
                          "TXNAME")
  anno_row <-
    data.frame(
      id = rownames(candidatesFINAL),
      Chromosome = gen_dat$TXCHROM,
      Symbol = gen_dat$SYMBOL
    )
  rownames(anno_row) <- anno_row[, 1]
}
if (calc && !is.null(anno_row)) {
  anno_row$Fisher <- factor(fisherVal)
  anno_row$Fisher[is.na(anno_row$Fisher)]= 1
  anno_row$Fisher = as.double(as.character(anno_row$Fisher))
}

## Col annotation
anno_col <- NULL
if (!is.null(group)) {
  anno_col <- data.frame(group)
  rownames(anno_col) <- colnames(candidatesFINAL)
}


##output
pheatmap::pheatmap(
  candidatesFINAL,
  annotation_row = anno_row[, -1, drop = FALSE],
  annotation_col = anno_col[],
  annotation_colors = anno_color,
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  show_rownames = FALSE
)

## get significant canididates list
if (calc) {
  anno_row[which(fisherVal < 0.1), ]
}
