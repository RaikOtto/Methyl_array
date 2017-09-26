# RUV

dif_var_ana = function(CASE,CTRL){

    print(c(CASE,CTRL))
    cohort_match = match(colnames(mSetSw), names(subtype[subtype %in% c(CASE,CTRL)]), nomatch = 0)
    cohort_match_case = colnames(mSetSw)[ match(colnames(mSetSw), names(subtype[subtype %in% c(CASE)]), nomatch = 0) ]
    cohort_match_ctrl = colnames(mSetSw)[ match(colnames(mSetSw), names(subtype[subtype %in% c(CTRL)]), nomatch = 0) ]
    subset_data = mSetSw[,cohort_match != 0]
    subset_data_rg = rgSet[,cohort_match != 0]
    #subset_data_rg = rgSet[,cohort_match != 0]
    
    grp = subtype[ cohort_match != 0]
    des <- model.matrix(~ 0 + grp)
    
    # get M-values for ALL probes
    
    meth <- getMeth(subset_data)
    unmeth <- getUnmeth(subset_data)
    Mval = log2((meth + 100)/(unmeth + 100))
    beta <- getBeta(subset_data)
    
    INCs <- getINCs(subset_data_rg)
    Mc <- rbind(Mval,INCs)
    ctl <- rownames(Mc) %in% rownames(INCs)
    
    rfit1 = RUVfit(data=Mc, design=des, coef=2, ctl=ctl)
    rfit2 = RUVadj(rfit1)
    
    top1 <- topRUV(rfit2, num=Inf)
    ctl <- rownames(Mval) %in% rownames(top1[top1$p.ebayes.BH > 0.5,])
    rfit1 <- RUVfit(data=Mval, design=des, coef=2, ctl=ctl) # Stage 2 analysis
    rfit2 <- RUVadj(rfit1)
    
    # variance ana
    
    fitvar = varFit(Mval, design = des, coef = 1:2)
    summary(decideTests(fitvar))
    topDV = topVar(fitvar, coef = 2, number = 200)
    summary( topDV )
    
    cpgsDV = rownames(topDV)
    annot_match = match( cpgsDV,  rownames(IlluminaHumanMethylationEPICanno.ilm10b2.hg19::Other ) , nomatch = 0)
    hgnc_names = IlluminaHumanMethylationEPICanno.ilm10b2.hg19::Other$UCSC_RefGene_Name[annot_match]
    hgnc_names = sapply( hgnc_names, FUN = function(vec){return(
      paste( as.character( unique(as.character(unlist(str_split(vec,pattern = ";"))))), collapse = ";", sep ="" )
    )} )
    
    var_case = round( apply( Mval[ match(res_tab$cpgsDV, rownames(Mval),nomatch = 0), which( colnames(Mval) %in% cohort_match_case)] , FUN =var, MARGIN = 1), 1)
    var_ctrl = round( apply( Mval[ match(res_tab$cpgsDV, rownames(Mval),nomatch = 0), which( colnames(Mval) %in% cohort_match_ctrl)] , FUN =var, MARGIN = 1), 1)
    dif_var = var_ctrl - var_case
    
    res_tab = as.data.frame(cbind(hgnc_names,cpgsDV,dif_var,var_case,var_ctrl,topDV[c("SampleVar","LogVarRatio","DiffLevene","Adj.P.Value")]))
    res_tab = res_tab[order(res_tab$dif_var, decreasing = T),]
    res_tab$SampleVar = round(res_tab$SampleVar,1);res_tab$LogVarRatio = round(res_tab$LogVarRatio,1);res_tab$DiffLevene = round(res_tab$DiffLevene,1)
    
    pdf(
      paste( c("~/Koop_Klinghammer/Results/Dif_Meth/Differential_variance/CASE_",CASE,"_CTRL_",CTRL,".pdf"), sep ="",collapse = ""),
      onefile = F, paper="a4r",width = 28, height = 18
    ) 
    par(mfrow=c(2,2))
    for( i in 1:4){
        
        
        boxplot(
            Mval[
                rownames(Mval) == cpgsDV[i],
                colnames(Mval) %in% names(subtype)[subtype %in% CASE]
            ],
            Mval[
                rownames(Mval) == cpgsDV[i],
                colnames(Mval) %in% names(subtype)[subtype %in% CTRL]
            ],
            main = paste( cpgsDV[i] ,hgnc_names[i], sep ="_"),
            names= c(CTRL,CASE)   
        )
        
    }
    dev.off()
    par(mfrow=c(1,1))
    
    write.table(
      res_tab,
      paste( c( "~/Koop_Klinghammer/Results/Dif_Meth/Differential_variance/Dif_var_methy_CASE_",CASE,"_CTRL_",CTRL,".tsv"), collapse = "", sep=""),
      quote = F,
      sep = "\t",
      row.names = F
    )

}

entry_list = list( data.frame( CASE = "MS", CTRL = "BA" ),data.frame( "CASE" = "MS", "CTRL" = "CL" ), data.frame("CASE" = "CL", "CTRL"="BA") )

for( i in 1:length(entry_list)){
  
    l = entry_list[[i]]
    print(l)
    dif_var_ana( CASE = as.character(l$CASE), CTRL = as.character(l$CTRL))
}
