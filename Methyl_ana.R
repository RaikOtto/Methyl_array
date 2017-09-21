rgSet <- read.metharray.exp(targets = targets)
mSet <- preprocessRaw(rgSet)
mSetSw <- SWAN(mSet,verbose=TRUE)

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
beta <- getBeta(mset_reduced)
dim(Mval)

#

s_names = targets$Sample_Name[ match( colnames(Mval), paste(targets$Slide, targets$Array, sep ="_") , nomatch = 0) ]
subtype_match = match(s_names, cohort_info$ID, nomatch = 0)
subtype = cohort_info$Group[subtype_match]
par(mfrow=c(1,1))
pdf("~/Koop_Klinghammer/Results/Methyl/Unsupervised_similarity.pdf",onefile = F, paper="a4r",width = 28, height = 18)
    plotMDS(Mval, labels=targets$Sample_Name, col = as.integer(factor(subtype)))
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
par(mfrow=c(2,2))
for(i in 1:3){
    stripchart(
        beta[ rownames(beta) == cpgs[i], ] ~ design[,1],
        method = "jitter",
        group.names=c("BA","CL","MS","Not_classified"),pch=16,cex=1.5,col=c(4,2),ylab="Beta values",
        vertical=TRUE,cex.axis=1.5,cex.lab=1.5
    )
    title(cpgs[i],cex.main=1.5)
}

# get M-values for ALL probes

meth <- getMeth(mSet)
unmeth <- getUnmeth(mSet)
M <- log2((meth + 100)/(unmeth + 100))

grp <- factor( subtypes,levels=c("MS","BA","CL"))
des <- model.matrix(~grp)
des

INCs <- getINCs(rgSet)
head(INCs)

Mc <- rbind(M,INCs)
ctl <- rownames(Mc) %in% rownames(INCs)
table(ctl)

rfit1 <- RUVfit(data=Mc, design=des, coef=2, ctl=ctl) # Stage 1 analysis
rfit2 <- RUVadj(rfit1)

top1 <- topRUV(rfit2, num=Inf)
head(top1)

ctl <- rownames(M) %in% rownames(top1[top1$p.ebayes.BH > 0.5,])
table(ctl)

rfit1 <- RUVfit(data=M, design=des, coef=2, ctl=ctl) # Stage 2 analysis
rfit2 <- RUVadj(rfit1)

# Look at table of top results
topRUV(rfit2)
table(rfit2$p.ebayes.BH < 0.01)

beta <- getBeta(mSet)
beta_norm <- rowMeans(beta[,des[,2]==0])
beta_can <- rowMeans(beta[,des[,2]==1])
Delta_beta <- beta_can - beta_norm
sigDM <- rfit2$p.ebayes.BH < 0.01 & abs(Delta_beta) > 0.25
table(sigDM)

topCpGs<-topRUV(rfit2,number=10000)
sigCpGs <- rownames(topCpGs)
sigCpGs[1:10]

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
gst <- gometh(sig.cpg=sigCpGs, all.cpg=rownames(rfit2), collection="GO")

topCpGs<-topRUV(rfit2,number=10000)
sigCpGs <- rownames(topCpGs)

## annotation

genes <- toTable(org.Hs.egSYMBOL2EG)
set1 <- sample(genes$gene_id,size=80)
set2 <- sample(genes$gene_id,size=100)
set3 <- sample(genes$gene_id,size=30)
genesets <- list(set1,set2,set3)
gsa <- gsameth(sig.cpg=sigCpGs, all.cpg=rownames(rfit2), collection=genesets)

topGSA(gsa)

#

source("https://bioconductor.org/biocLite.R")
biocLite("IlluminaHumanMethylation450k.db")

library(IlluminaHumanMethylation450k.db)
CpG_annotation <- as.list(IlluminaHumanMethylation450kSYMBOL[mappedkeys(IlluminaHumanMethylation450kSYMBOL)])

mapToGenome("cg07155336")
require(minfiData)
GMsetEx.sub <- mapToGenome(MsetEx.sub)

GRanges(GMsetEx.sub@assays)

biocLite('FDb.InfiniumMethylation.hg19')
library('FDb.InfiniumMethylation.hg19')

hm450 <- get450k()
probenames <- c("cg16392865", "cg00395291", "cg09310185", "cg21749424")
probes <- hm450[sigCpGs]
getNearestTSS(probes)
