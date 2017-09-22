MSet <- preprocessRaw(rgSet) 
RSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
beta <- getBeta(RSet)

GRset <- mapToGenome(RSet)
beta <- getBeta(GRset)
M <- getM(GRset)
CN <- getCN(GRset)

plotQC(getQC(MSet))

GRset.funnorm = readRDS("~/Koop_Klinghammer/Data/GRset.funnorm.RData")
#GRset.funnorm <- preprocessFunnorm(rgSet)
#saveRDS(GRset.funnorm, "~/Koop_Klinghammer/Data/GRset.funnorm.RData")

# QC

qc <- getQC(mSetSw)
plotQC(qc, id=colnames(mSetSw))
densityPlot(mSetSw, sampGroups = subtype)
densityBeanPlot(mSet, sampGroups = subtype)
controlStripPlot(rgSet, controls="BISULFITE CONVERSION II")
qcReport(rgSet, pdf= "~/Koop_Klinghammer/Results/qcReport.pdf")

GRset <- mapToGenome(mSetSw)
predictedSex <- getSex(GRset, cutoff = -2)$predictedSex
head(predictedSex)
plotSex(getSex(GRset, cutoff = -2), id = s_names)

# Dif exp

pheno = subtype
designMatrix <- model.matrix(~ pheno)
dmrs <- bumphunter(GRset.funnorm, design = designMatrix,cutoff = 0.2, B=0, type="Beta")
dmrs$coef

ab <- compartments(grset.quantile, chr="chr14", resolution=100-1000)

#
par(mfrow=c(1,2), cex=1.25)
densityByProbeType(mSet[,1], main = "Raw")
densityByProbeType(mSetSw[,1], main = "SWAN")

# filtering outlier

detP <- detectionP(rgSet)
keep <- rowSums(detP < 0.01) == ncol(rgSet)
mSetSw <- mSetSw[keep,]

# extract values

meth <- getMeth(mSetSw)
unmeth <- getUnmeth(mSetSw)
Mval <- log2((meth + 100)/(unmeth + 100))

beta <- getBeta(mSetSw)
dim(Mval)

#

par(mfrow=c(1,1))

vis_match = match(colnames(Mval), paste(targets$Slide, targets$Array, sep = "_"), nomatch = 0)
vis_names = targets$Sample_Name[vis_match]
cohort_info

pdf("~/Koop_Klinghammer/Results/Methyl/Unsupervised_similarity.pdf",onefile = F, paper="a4r",width = 28, height = 18)
    plotMDS(Mval, labels = vis_names, col = as.integer(factor(subtype)))
    legend("topleft",legend=c("BA","CL","MS","Not_classified"),pch=16,cex=1.2,col=1:4)
dev.off()

#

group  = factor( subtype,levels=c("BA","CL","MS","Not_classified"))
id     = factor(targets$Sample_Name)
design = model.matrix(~ group)
design

#

fit = lmFit(Mval,design)
fit = eBayes(fit)
summary(decideTests(fit))

top = topTable(fit,coef = 3)
top

cpgs <- rownames(top)
beta <- getBeta(mSetSw)

par(mfrow=c(2,2))

for(i in 1:3){
    stripchart(
        beta[ rownames(beta) == cpgs[i], ] ~ design[,1],
        method = "jitter",
        group.names = c("BA","CL","MS","Not_classified")[i],
        pch = 16,
        cex = 1.5,
        col = c(4,2),
        ylab = "Beta values",
        vertical = TRUE,
        cex.axis = 1.5,
        cex.lab = 1.5
    )
    title(cpgs[i],cex.main=1.5)
}

# get M-values for ALL probes

CASE = "CL"
CTRL = "BA"

cohort_match = match(colnames(GRset.funnorm), names(subtype[subtype %in% c(CASE,CTRL)]), nomatch = 0)
subset_data = mSetSw[,cohort_match != 0]
subset_data_rg = rgSet[,cohort_match != 0]

grp = subtype[ cohort_match != 0]
des <- model.matrix(~ grp)

meth <- getMeth(subset_data)
unmeth <- getUnmeth(subset_data)
Mval <- log2((meth + 100)/(unmeth + 100))
beta <- getBeta(subset_data)
dim(Mval)

fit <- lmFit(Mval,des)
fit <- eBayes(fit)
summary(decideTests(fit))

top = topTable(fit,coef=2, sort.by = "logFC", p.value = .05, number = 300)
sigCpGs = rownames(top)

library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19 )

annot_match = match( sigCpGs,  rownames(IlluminaHumanMethylationEPICanno.ilm10b2.hg19::Other ) , nomatch = 0)
hgnc_names = IlluminaHumanMethylationEPICanno.ilm10b2.hg19::Other$UCSC_RefGene_Name[annot_match]
hgnc_names = sapply( hgnc_names, FUN = function(vec){return(
    paste( as.character( unique(as.character(unlist(str_split(vec,pattern = ";"))))), collapse = ";", sep ="" )
)} )

res_tab = as.data.frame(cbind(hgnc_names,sigCpGs,top[c("logFC","adj.P.Val")]))
res_tab = res_tab[order(res_tab$logFC, decreasing = T),]

m_data = getM( mSetSw )

boxplot(
    m_data[ 
        rownames(m_data) == res_tab$sigCpGs[1],
        colnames(m_data) %in% names(subtype)[subtype %in% CASE]
    ],
    m_data[ 
      rownames(m_data) == res_tab$sigCpGs[1],
      colnames(m_data) %in% names(subtype)[subtype %in% CTRL]
    ],
    main = res_tab$sigCpGs[1],
    names= c(CTRL,CASE)   
)
res_tab$logFC = round(res_tab$logFC,1)

write.table(
    res_tab,
    paste( c( "~/Koop_Klinghammer/Results/Dif_Meth/Dif_methy_CASE_",CASE,"_CTRL_",CTRL,".tsv"), collapse = "", sep=""),
    quote = F,
    sep = "\t",
    row.names = F
)
