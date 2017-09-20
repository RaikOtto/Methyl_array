#Function "centroid2expr": apply the predictor centroid to validation expression data matrix
library("stringr")

# pure data pre-processing

pure_data = read.table(
  "~/Downloads/Expression_HNSC.txt",
  header = T,
  sep = "\t",
  stringsAsFactors = F
)

library("hgu133plus2.db")
hgnc_genes = as.character( unlist(mget( pure_data$NAME, hgu133plus2SYMBOL ) ) )
rownames( pure_data ) = pure_data$NAME

pure_data = pure_data[ ,c(-1,-2) ]

#colnames( pure_data ) = str_replace( colnames( pure_data ), pattern = "(^X)", ""  )
#colnames( pure_data ) = str_replace_all( colnames( pure_data ), pattern = ".CEL", "" )

colnames( pure_data ) = c(
  "10110",
  "10114",
  "10159",
  "10309",
  "10379",
  "10511",
  "10621",
  "10632",
  "10913",
  "10922",
  "10924",
  "10927",
  "10960",
  "10980B",
  "11097",
  "11142",
  "11143",
  "11218",
  "11269A",
  "11269B",
  "11303",
  "11437A",
  "11437B",
  "11452",
  "11482",
  "9619",
  "9876",
  "9897"
  )

# load centroid data

balanced.centroid = read.table(
  "~/Koop_Klinghammer//balanced.centroid.txt",
  header=TRUE,
  row.names = 1,
  sep = "\t",
  stringsAsFactors = F
)

pure_data = pure_data[ str_to_upper( hgnc_genes ) %in% str_to_upper( rownames( balanced.centroid ) ), ]
hgnc_genes= hgnc_genes[ str_to_upper( hgnc_genes ) %in% str_to_upper( rownames( balanced.centroid ) ) ]

# qa check

str_to_upper( rownames( balanced.centroid ) )[ which(  !( str_to_upper( rownames( balanced.centroid ) ) %in% str_to_upper( hgnc_genes ) ) ) ]

### zero check

var_row_equal = function( i_vec )( return( var( as.double( i_vec ) ) < 0.1  )  )
zero_var = apply(
  pure_data,
  FUN = var_row_equal,
  MARGIN = 1
)
pure_data = pure_data[ ! as.logical( zero_var ),]
hgnc_genes = hgnc_genes[ ! as.logical( zero_var ) ]

### variance selection of most diverging probe for same genes

variance_selection = function(  elem, hgnc_genes ){
  
  mapping = hgnc_genes == as.character( elem )
  exprs_vals = matrix(
    as.double(
        unlist( 
            pure_data[ mapping,  ]
        )
    ),
    nrow = sum(mapping)
  )
  exprs_var = apply(
    exprs_vals,
    MARGIN = 1,
    FUN = var 
  )
  selector  = which( exprs_var == max( exprs_var )  )[1]
  
  max_var_gene = matrix(
    c(
      elem,
      as.character( exprs_vals[ selector, ] ) 
    ), 
    ncol = length( colnames( pure_data ) ) +1
  )
  
  return( max_var_gene )
}
# YME1L1

uni_genes = unique( hgnc_genes )
uni_genes = uni_genes[ uni_genes != ""  ]

result_mat = sapply(
  uni_genes,
  FUN = variance_selection,
  hgnc_genes
)
result_mat = t( result_mat )
result_mat = result_mat[, c(-1)]
result_mat = t( apply( result_mat, FUN = as.double, MARGIN=1 ) )
colnames(result_mat) = colnames(pure_data)

### cohorts

cohort_t = read.table(
  "~//Koop_Klinghammer//upload/Cel_files/cohorts.tab",
  header = T,
  sep = "\t",
  stringsAsFactors = F
)

match_names = match( 
  colnames(result_mat),
  str_replace_all( cohort_t$ID, pattern = ".CEL","")
)

classifier_groups = cohort_t$Subtype[ match_names ]

pure_data = pure_data[ , order( classifier_groups) ]

classifier_groups = classifier_groups[ order( classifier_groups ) ]

### mean_mat

mean_mat <<- matrix( as.double(), ncol = length(colnames(result_mat))  )

create_dif_mat = function( vector ){
  
  mean_mat <<- rbind(
    mean_mat,
    vector - mean(vector)
  )
}
apply( result_mat, FUN = create_dif_mat, MARGIN = 1  )
rownames(mean_mat) = rownames(result_mat  )
