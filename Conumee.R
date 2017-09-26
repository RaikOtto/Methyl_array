library("minfi")
library("conumee")
#RGsetTCGA <- read.450k.exp(base = "~/conumee_analysis/jhu-usc.edu_BRCA.HumanMethylation450.Level_1.8.8.0")  # adjust path
RGsetTCGA <- read.450k.url()  # use default parameters for vignette examples
MsetTCGA <- preprocessIllumina(RGsetTCGA)
MsetTCGA
## class: MethylSet 
## dim: 485512 2 
## metadata(0):
## assays(2): Meth Unmeth
## rownames(485512): cg00050873 cg00212031 ... ch.22.47579720R
##   ch.22.48274842R
## rowData names(0):
## colnames(2): 6042324037_R05C02 6042324037_R06C01
## colData names(0):
## Annotation
##   array: IlluminaHumanMethylation450k
##   annotation: ilmn12.hg19
## Preprocessing
##   Method: Illumina, bg.correct = TRUE, normalize = controls, reference = 1
##   minfi version: 1.22.0
##   Manifest version: 0.4.0