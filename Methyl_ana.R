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


# bumphunter

# Dif exp

pheno = subtype
designMatrix <- model.matrix(~ pheno)
dmrs <- bumphunter(GRset.funnorm, design = designMatrix,cutoff = 0.2, B=0, type="Beta")
dmrs$coef

ab <- compartments(grset.quantile, chr="chr14", resolution=100-1000)
