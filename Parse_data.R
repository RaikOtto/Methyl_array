library("cnAnalysis450k")
library("org.Hs.eg.db")
library("stringr")
library("missMethyl")
library("limma")
library("minfi")
library("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
library("grid")     ## Need to attach (and not just load) grid package
draw_colnames_45 <- function (coln, gaps, ...) {coord = pheatmap:::find_coordinates(length(coln), gaps);x = coord$coord - 0.5 * coord$size;  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 45, gp = gpar(...));  return(res)}
assignInNamespace(x="draw_colnames", value="draw_colnames_45",ns=asNamespace("pheatmap"))

targets = read.metharray.sheet( base = file.path("~/Koop_Klinghammer/Misc/"), pattern = "PDXsamplesheet130516.csv")
targets[,1:9]
raw_files = list.files("~/Koop_Klinghammer/Data/",pattern =".idat",recursive = T,full.names = T)
raw_file_names = list.files("~/Koop_Klinghammer/Data/",pattern =".idat",recursive = T,full.names = F)
raw_file_names = str_extract_all(raw_file_names, regex("/.*_.*_"))
raw_file_names = str_replace_all(raw_file_names,regex("(_$)|(^/)"), "")

raw_match = match(
  paste( targets$Slide,targets$Array, sep ="_" ),
  raw_file_names,
  nomatch = 0
)
targets$Basename = raw_files[raw_match]

cohort_info = read.table("~/Koop_Klinghammer/Misc/samples_methylierung.tsv",sep ="\t", header = T, stringsAsFactors = F)
meta_match= match( targets$Sample_Name, cohort_info$ID , nomatch = 0)
subtype = cohort_info$Group[meta_match  ]
subtype[is.na(subtype) ] = "Not_classified"

targets_all = targets
targets = targets[meta_match != 0,]

s_names =  cohort_info$ID[meta_match ]
names(subtype) = s_names

dim(targets)
#rgSet <- read.metharray.exp(targets = targets)
pData(rgSet)$subtype = subtype

mSetSw = readRDS("~/Koop_Klinghammer/Data/rg_Swan_normalized.RData")
mSet = readRDS("~/Koop_Klinghammer/Data/mSet.RData")
#saveRDS(mSet, "~/Koop_Klinghammer/Data/mSet.RData")
#saveRDS(mSetSw, "~/Koop_Klinghammer/Data/rg_Swan_normalized.RData")

#mSet <- preprocessRaw(rgSet)
#mSetSw <- SWAN(mSet,verbose=TRUE)
colnames(rgSet) = s_names
colnames(mSetSw) = s_names
colnames(mSet) = s_names

###

RSet <- ratioConvert(mSet, what = "both", keepCN = TRUE)
beta <- getBeta(RSet)

GRset <- mapToGenome(RSet)
beta <- getBeta(GRset)
M <- getM(GRset)
CN <- getCN(GRset)

#plotQC(getQC(MSet))

GRset.funnorm = readRDS("~/Koop_Klinghammer/Data/GRset.funnorm.RData")
#GRset.funnorm <- preprocessFunnorm(rgSet)
#saveRDS(GRset.funnorm, "~/Koop_Klinghammer/Data/GRset.funnorm.RData")

# filtering outlier

detP <- detectionP(rgSet)
keep <- rowSums(detP < 0.01) == ncol(rgSet)
mSetSw_raw = mSetSw
mSetSw <- mSetSw_raw[keep,]

# extract values

meth <- getMeth(mSetSw)
unmeth <- getUnmeth(mSetSw)
Mval <- log2((meth + 100)/(unmeth + 100))

beta <- getBeta(mSetSw)
dim(Mval)
