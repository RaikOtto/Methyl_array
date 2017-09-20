#Function "centroid2expr": apply the predictor centroid to validation expression data matrix
library("ggbiplot")
library("rgl")
library("pheatmap")
#library("PoiClaClu")
library("stringr")
library("biomaRt")
library("GEOquery")

keck_mode = TRUE

balanced.centroid = read.table( "~/Koop_Klinghammer/Misc//balanced.centroid.txt", header=TRUE, row.names=1, sep="\t",stringsAsFactors = F)
balanced.centroid_importance = sort(rowSums(abs(balanced.centroid)), decreasing = T)
balanced.centroid = balanced.centroid[ match(names(balanced.centroid_importance),rownames(balanced.centroid)),]

# GSE40774

gse = getGEO("GSE23036",GSEMatrix = T)#GSE23036 # GSE3292 <- both hgus # GSE73339 # GSE40774
#pure_data = exprs(gse$)
annot_table = read.table("~/Koop_Klinghammer/Misc/GPL13497-9755.txt",sep ="\t", header = T, stringsAsFactors = F, comment.char = "#", fill = T)
hgnc_genes = annot_table$GENE_SYMBOL

colnames(pure_data) = str_replace(colnames(pure_data), pattern  ="^X", "");colnames(pure_data) = str_replace(colnames(pure_data), pattern  =".CEL", "")
col_var = apply(as.matrix(pure_data),FUN = function(vec){return (var(as.double(vec) ))},MARGIN = 2)
row_var = apply(as.matrix(pure_data),FUN = function(vec){return (var(as.double(vec) ))},MARGIN = 1)
pure_data = pure_data[row_var > 0, col_var > 0]

### GEO Klinghammer

#pure_data = read.table("~/Koop_Klinghammer/upload/Normalized_expresssion_15_06_2017.tsv",sep="\t", header = T, stringsAsFactors = F, comment.char = "!", row.names = 1)
pure_data = read.table("~/Koop_Klinghammer/Misc/GSE3292_series_matrix.txt",sep="\t", header = T, stringsAsFactors = F, comment.char = "!", row.names = 1)
colnames(pure_data) = str_replace(colnames(pure_data), pattern  ="^X", "");colnames(pure_data) = str_replace(colnames(pure_data), pattern  =".CEL", "")
col_var = apply(as.matrix(pure_data),FUN = function(vec){return (var(as.double(vec) ))},MARGIN = 2)
row_var = apply(as.matrix(pure_data),FUN = function(vec){return (var(as.double(vec) ))},MARGIN = 1)
pure_data = pure_data[row_var > 0, col_var > 0]
pure_data = pure_data[, which(!(colnames(pure_data) %in% c("10922","10980A")))]

library(hgu133plus2.db)

rm = rowMeans(pure_data)
summary(rm)
pure_data = pure_data[rm >= summary(rm)[2],]
col_var = apply(as.matrix(pure_data),FUN = function(vec){return (var(as.double(vec) ))},MARGIN = 2)
row_var = apply(as.matrix(pure_data),FUN = function(vec){return (var(as.double(vec) ))},MARGIN = 1)
summary(row_var)
pure_data = pure_data[row_var > summary(row_var)[2], col_var > 0]
pure_data[1:5,1:5]

# time for variance filtering

hgnc_genes = as.character(unlist(mget( rownames(pure_data), hgu133plus2SYMBOL )))
pure_data = pure_data[ !is.na(hgnc_genes),];
hgnc_genes = hgnc_genes[ !is.na(hgnc_genes)  ]
hgnc_list = hgnc_genes
hgnc_list_uni = unique(hgnc_list[!is.na(hgnc_list)])

max_list <<- c()
for (gene in hgnc_list_uni){
  var_match = which( hgnc_list == gene )
  row_var_max = which.max( apply(as.matrix(pure_data[ var_match,]),FUN = function(vec){return (var(as.double(vec) ))},MARGIN = 1) )
  max_list = c( max_list, var_match[row_var_max])
}
expr = pure_data[max_list,]
rownames(expr) = hgnc_list[max_list]
dim(expr)
expr[1:5,1:5]

#ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")
#agilent.df = getBM(attributes = c("hgnc_symbol","agilent_wholegenome_4x44k_v2"), filters=c("agilent_wholegenome_4x44k_v2"),values=rownames(pure_data), mart=ensembl)

#affy.df = getBM(attributes = c("hgnc_symbol","affy_hg_u133_plus_2"), filters=c("affy_hg_u133_plus_2"),values=rownames(pure_data), mart=ensembl)
#new_match = match( affy.df$affy_hg_u133_plus_2, rownames(pure_data), nomatch = 0)
#new_mat = pure_data[new_match,]; rownames(new_mat) = 1:nrow(new_mat)
#hgnc_genes = affy.df$hgnc_symbol

# filter for unwanted samples

###

source("~/Koop_Klinghammer/Scripts/Classification_scripts.R")
expr = expr[ which( str_to_upper( rownames(expr) ) %in% str_to_upper( rownames(balanced.centroid) ) ) , ]

pub_cor <<- matrix( as.double(), ncol = length( colnames( balanced.centroid )  ) )
expr2bc = centroid2expr( balanced.centroid[,], expr )
groups = as.character( expr2bc$correlation[,2] )
source("~/Koop_Klinghammer/Scripts/setup_meta_data.R")
colnames(expr2bc$correlation) = c("ID","Group","Cor","PValue")
q_value = p.adjust(as.double(as.character(as.data.frame(expr2bc$correlation)[,4])),"BH")
expr2bc$correlation = cbind( expr2bc$correlation , q_value)
expr2bc$correlation[order(expr2bc$correlation[,2]),]

write.table(expr2bc$correlation,"~/koop_klinghammer/Misc/GSE3292_all_genes.tsv",sep ="\t",quote =F,row.names =F)
### pca

source("~/Koop_Klinghammer/Scripts/setup_meta_data.R")
meta_data$Cetuximab = as.double(meta_data$Cetuximab )
meta_data$Cetuximab[ is.na(meta_data$Cetuximab)] = -1.0
colfunc<-colorRampPalette(c("black","green","white","red"))
cetuximab_colors = list( Cetuximab = colfunc(nrow(meta_data)) )
cetuximab_colors = as.character(unlist(cetuximab_colors))
cetuximab_colors[ is.na(meta_data$Cetuximab)] = "#000000"
cetuximab_colors = list(Cetuximab = cetuximab_colors)

meta_data$Subtype = as.character(expr2bc$correlation[match(rownames(meta_data),expr2bc$correlation[,1]),2])

#pdf("~/Dropbox/PhD/koop_klinghammer/gene expression paper/Grafiken/Vorschlaege/pHeatmap_300_sig_genes_vorschlag_raik.pdf", onefile = F)
df = as.data.frame(groups); rownames(df) = colnames(expr)
pheatmap::pheatmap(
  cor(expr[ rownames(expr) %in% rownames(balanced.centroid)[],]),
  annotation_col = df["groups"],
  method = "mean",
  clustering_method = "single",
  #annotation_colors = cetuximab_colors,
  cluster_cols = T
)
#dev.off()
#
eset_pca = prcomp(
  t(expr[ rownames(expr) %in% rownames(balanced.centroid),]),
  center = TRUE,
  scale = TRUE
)

g = ggbiplot(
  eset_pca,
  obs.scale = 1, 
  var.scale = 1, 
  labels.size = 4,
  alpha = 1,
  groups = groups,
  ellipse = TRUE, 
  circle = TRUE,
  var.axes = F
)
g = g + scale_color_discrete(name = '')
g = g + theme(legend.direction = 'horizontal', legend.position = 'top')
g
#pdf("~/Dropbox/PhD/Koop_Klinghammer/gene expression paper/Grafiken/final/Classification_without_10980A_and_11303_based_on_821_sig_genes.pdf")
  print(g)
#dev.off()

# heatmap 3
library("devtools")
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
require(RColorBrewer)

cols = colorRampPalette(brewer.pal(10, "RdBu"))(256)

distCor <- function(x) as.dist(1-cor(t(x)))
hclustAvg <- function(x) hclust(x, method="average")

if (keck_mode)
  pdf("~/Dropbox/PhD/Koop_Klinghammer/gene expression paper/Grafiken/final/Classification.pdf")

heatmap.3( 
  pub_cor[ order(pub_cor[,1], decreasing = T),],
  Rowv = TRUE,
  Colv = TRUE,
  #dendrogram = "column",
  trace = "none",
  #scale = "row",
  col = rev(cols),
  hclustfun = hclustAvg
)
dev.off()

### exclude 66/73 genes

mice_t = read.table("~/Dropbox/PhD/Koop_Klinghammer/Lit/mice_stroma_genes.tab",sep = "\t", header = T, stringsAsFactors = F)
intersect_genes = as.character( mice_t$überlappende.Gene )
intersect_genes = intersect_genes[ intersect_genes != ""]

names_input = str_replace( colnames( pure_data )[-c(1,2)], pattern = "(^X)", "" )
match_names = match( names_input, str_replace_all(cohort_t$ID,pattern = ".CEL",""))
groups = cohort_t$Subtype[match_names]

match =match(intersect_genes,pure_data$HGNC)
sub_t = pure_data_mat[match,]
colnames(sub_t) = colnames( pure_data )[ -c(1,2)]
colnames(sub_t) = str_replace( colnames(sub_t), pattern = ".CEL", ""  )
rownames( sub_t ) = pure_data$HGNC[ match ]
avg_t = sub_t - rowMeans(sub_t)

# top bar start

subtype_top_bar = as.matrix(
  c(  
    colorRampPalette(
      colors = c("blue")
    )( length( groups ) )
  ) 
)

subtype_top_bar[ groups == "BA",  ] = "yellow"
subtype_top_bar[ groups == "CL",  ] = "green"
subtype_top_bar[ !( groups %in% c("CL","BA","MS") ),  ] = "black"
#subtype_top_bar = subtype_top_bar[ groups %in% c("","CL","NA"),  ] = "yellow"

colnames(subtype_top_bar) = c("Subtype")

# top bar end
pdf("~/Dropbox/PhD/Koop_Klinghammer/gene expression paper/Grafiken/Absolute_exp_74_probes_mice_stroma.pdf");
heatmap.3( 
  sub_t,
  Rowv = TRUE,
  Colv = TRUE,
  #dendrogram = "both",
  trace = "none",
  #scale = "row",
  col = rev(cols),
  hclustfun = hclustAvg,
  main = "Absolute Expression 74 probes",
  ColSideColors = subtype_top_bar
)
dev.off()

pdf("~/Dropbox/PhD/Koop_Klinghammer/gene expression paper/Grafiken/Average_exp_74_probes_mice_stroma.pdf");
heatmap.3( 
  avg_t,
  Rowv = TRUE,
  Colv = TRUE,
  #dendrogram = "both",
  trace = "none",
  #scale = "row",
  col = rev(cols),
  hclustfun = hclustAvg,
  main = "Absolute Expression - average expression 74 probes",
  ColSideColors = subtype_top_bar
)
dev.off()

avg_t_plot = avg_t[, !is.na(groups) & (groups != "")]
groups_pca = groups[ !is.na(groups) & (groups != "") ]

rownames(avg_t_plot) = 1:dim(avg_t_plot)[1]

avg_t_pca = prcomp(
  t(avg_t_plot),
  center = TRUE,
  scale. = TRUE
)


g = ggbiplot( 
  avg_t_pca, 
  obs.scale = 1, 
  var.scale = 1, 
  groups = as.character( groups_pca ),
  ellipse = TRUE,
  circle = TRUE,
  var.axes = F,
  labels = colnames(avg_t_plot)
)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
pdf("~/Dropbox/PhD/Koop_Klinghammer/gene expression paper/Grafiken/pca_74_probes.pdf")
print(g)
dev.off()

# distance analyses
library("PoiClaClu")
library("pheatmap")
library("RColorBrewer")

cols = colorRampPalette(brewer.pal(10, "RdBu"))(256)

dist_n = dist(t(result))
dist_mat = as.matrix(dist_n)

df = data.frame(
    "Subtype" = groups
)

annotation = data.frame( "Sub_type" = factor( groups ) )
rownames(annotation) = colnames( result )
Var1        <- rep( "darkgreen", length(groups))
names(Var1) <- rep( "Sub_type", length(groups) )
anno_colors <- list(Var1 = Var1)

poisd <- PoissonDistance(t(result))
pois_m = as.matrix(poisd$dd)
colnames( pois_m ) = rownames(t(result))
rownames( pois_m ) = colnames( result  )

pheatmap(
    pois_m,
    #clustering_distance_rows=dist_n,
    #clustering_distance_cols=dist_n,
    #color = (cols),
    #annotation_col=anno_colors,
    annotation = annotation
)

