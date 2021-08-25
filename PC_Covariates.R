#From expression matrix (exp) and meta data (meta), compute covariate vs PC significance, create results. 
#The significance is computed using correlation method for numeric covariates, and ANOVA for categorical covariates.
#Output files (Excel file and PDF plots) will be saved in output directory with prefix as specified by the user (out_prefix). If out_prefix=NULL, will only return results as R list object. 
#exp is expression matrix (logTPM, logCPM, etc), meta is meta data. Column names (samples) of exp must match row names of meta. 
#PC_cutoff: select which principle components to be used in analysis. Default 5 will select components that explain more than 5% of variance in the data.
#FDR_cutoff will choose the covariate-PC pairs that pass this cutoff.
#N_col: number of columns for facet plots
#Besides the output files, return R list object including the following:
#data.all (PC scores and meta data combined), selVar_All (significant covariate-PC pairs)
#sel_dataN (Numeric covariate results), sel_dataC (Categorical covariate results)
Covariate_PC_Analysis<-function(exp, meta, out_prefix, PC_cutoff=5, FDR_cutoff=0.1, N_col=3, PCA_plots="ranked") {
  require(tidyverse);  require(cowplot); require(openxlsx)
  if (!is.null(out_prefix) ) {
    out_dir=dirname(out_prefix)
    if (!dir.exists(out_dir)) {dir.create(out_dir)}
  }
  res1=Compute_Corr_Anova(exp, meta, PC_cutoff=PC_cutoff)
  sel_dataN<-get_PC_meta_plot(res1, 'numeric', FDR_cutoff, N_col=N_col)
  selVarN=sel_dataN$selVar
  if (!is.null(selVarN) && !is.null(out_prefix)) { #output correlation results
    graphH=ceiling(nrow(selVarN)/3)*4
    graphW=min(nrow(selVarN)*4+1, 12)
    ggsave(str_c(out_prefix, "_Significant_Numeric_Covariates.pdf"), sel_dataN$plot, width=graphW, height=graphH)
  }  
  sel_dataC<-get_PC_meta_plot(res1, 'categorical', FDR_cutoff, N_col=N_col)
  selVarC=sel_dataC$selVar
  if (!is.null(selVarC) && !is.null(out_prefix)) { #output anova results
    graphH=ceiling(nrow(selVarC)/3)*4
    graphW=min(nrow(selVarC)*4+1, 12)
    ggsave(str_c(out_prefix, "_Significant_Categorical_Covariates.pdf"), sel_dataC$plot, width=graphW, height=graphH)
  }  
  
  if (!is.null(selVarC)) {selVarC<-selVarC%>%arrange(fdr)}
  if (!is.null(selVarN)) {selVarN<-selVarN%>%arrange(fdr)}
  
  PC_info<-data.frame(PC=colnames(res1$PC_scores), Per=res1$percentVar[1:ncol(res1$PC_scores)])%>%mutate(PC_new=str_c(PC, " (", Per, "%)"))
  PC_Scores=res1$PC_scores; PC_Scores=data.frame(Sample=rownames(PC_Scores), PC_Scores)
  MetaData=res1$meta; MetaData=data.frame(Sample=rownames(MetaData), MetaData)
  if (!is.null(out_prefix) ) {
    dat=list(Sig_Cat_Anova=selVarC, Sig_Num_corr=selVarN, All_Cat_Anova=res1$anova, All_Num_Corr=res1$corr, PC_Scores=PC_Scores, MetaData=MetaData, 
             PC_Perentage=PC_info[, 1:2])
    write.xlsx(dat, str_c(out_prefix, file='_Covariate_PC_Results.xlsx'))
  }
  
  #now create PCA plots
  data.all<-cbind(res1$PC_scores, res1$meta)
  if (is.null(selVarN)) {selVarN1=NULL
  } else { selVarN1<-selVarN%>%mutate(NewText=str_c(covar, " vs. ", PC, ": ", text ))%>%
    mutate(Type="Numeric")%>%dplyr::select(PC, covar, Type, NewText, pvalue, fdr) }
  if (is.null(selVarC)) {selVarC1=NULL
  } else { selVarC1<-selVarC%>%mutate(NewText=str_c(covar, " vs. ", PC, " ", text ))%>%
    mutate(Type="Categorical")%>%dplyr::select(PC,  covar, Type,NewText, pvalue, fdr) }
  selVar_All<-rbind(selVarC1, selVarN1)
  
  if (!is.null(selVar_All)) {
    selVar_All<-selVar_All%>%arrange(fdr) #sort by FDR
    if (nrow(selVar_All)>0 && !is.null(out_prefix) ) {
      pdf(str_c(out_prefix, "_PCA_Plots.pdf"), width=8, height=9)
      
      if (PCA_plots=="ranked") { #plot PCA plots in order of FDR
        for (i in 1:nrow(selVar_All)) {
          x0=selVar_All$PC[i]
          if (x0=="PC1") {y0="PC2"} else {y0=x0; x0="PC1"}
          x=sym(x0); y=sym(y0); color_by=sym(selVar_All$covar[i])
          p<-ggplot(data.all, aes(x=!!x, y=!!y, col=!!color_by))+geom_point()+
            labs(x=PC_info$PC_new[PC_info$PC==x0], y=PC_info$PC_new[PC_info$PC==y0])+ theme_half_open()
          #additional text to add  
          p<-add_sub(p, selVar_All$NewText[i], x=0.2, hjust=0)
          print(ggdraw(p))
        }
      }  else { #plot PCA plots for each variate, try to minimize plots
        var_list=as.character(sort(unique(selVar_All$covar)))
        for (i in 1:length(var_list) ) {
          var=as.character(var_list[i])
          PCs<-selVar_All%>%filter(covar==var)%>%dplyr::select(PC)%>%unlist%>%as.character()
          if (length(PCs)==1) {
            PCs=sort(c(PCs, ifelse(PCs=="PC1", "PC2", "PC1")))
          }
          for  (j in 2:length(PCs)) {
            x=sym(PCs[1]); y=sym(PCs[j]); color_by=sym(var)
            p<-ggplot(data.all, aes(x=!!x, y=!!y, col=!!color_by))+geom_point()+
              labs(x=PC_info$PC_new[PC_info$PC==PCs[1]], y=PC_info$PC_new[PC_info$PC==PCs[j]])+ theme_half_open()
            #additional text to add  
            more_text<-selVar_All%>%filter(covar==var)%>%dplyr::select(NewText)%>%unlist()%>%paste(collapse="\n")
            p<-add_sub(p, str_c(more_text, "\n(This plot shows ", PCs[1], " in X and ", PCs[j], " in Y)"), x=0.2, hjust=0)
            print(ggdraw(p))
          }
        }
      }
      dev.off()
    }
  }

  if (!is.null(out_prefix) ) {cat("Please check output files at:", out_dir, "\n")}
  if (!is.null(selVar_All)) {
    names(selVar_All)[c(2, 4, 5, 6)]=c("Covariate", "Significance", "P-value", "FDR")
  }

  return(list(data.all=data.all, selVar_All=selVar_All, sel_dataN=sel_dataN, sel_dataC=sel_dataC, ncol=N_col))
}

#Plot one covariate vs. one principle component. Can be useful when there are too many categories for the default faceted plot.
#Input: res is from function Covariate_PC_Analysis. add_text will add a line below the plot to show if the covariate vs PC is significant. 
plot_covariate_PC<-function(res, pc, var, out_file, width=10, height=8, add_text=TRUE) {
  data.all=res$data.all
  selVar=res$selVar_All
  if (!(var %in% names(data.all))) {cat("Covariate ", var, " not in MetaData. Please check the spelling of covariate.\n", sep=""); return(NULL)}
  if (!(pc %in% names(data.all))) {cat(pc, " not in principle component scores. Please check the spellin.\n", sep=""); return(NULL)}
  Num_names=names(select_if(data.all, is.numeric))
  if (var %in% Num_names) { #numeric covariate
    p<-ggplot(data.all, aes(x=!!sym(var), y=!!sym(pc)) )+geom_point()+
      stat_summary(fun.data= mean_cl_normal) + geom_smooth(method='lm')+theme_half_open()
  } else {
    p<-ggplot(data.all,  aes(x=!!sym(var), y=!!sym(pc)) )+geom_boxplot()+geom_jitter(alpha=0.7, width=0.1)+theme_half_open() +
      theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5))
  }
  if (add_text) {
    info<-selVar%>%filter(PC==pc, covar==var)
    if ( nrow(info)>0) {text_info=info$NewText[1]} else {text_info=str_c(var, " vs. ", pc, " not significant.")}
    p<-add_sub(p, text_info, x=0.2, hjust=0)
  }
  pdf(out_file, width=width, height=height)
  print(ggdraw(p))
  dev.off()
  cat("Please check output at: ", out_file)
}


#compute coraviate vs PC significance. For numerica, use correlation, for categorical, use ANOVA.
#exp is expression matrix (logTPM, logCPM, etc), meta is meta data. Column names (samples) of exp must match row names of meta. 
#PC_cutoff: select which principle components to be used in analysis. Default 5 will select components that explain more than 5% of variance in the data.
Compute_Corr_Anova<-function(exp, meta, PC_cutoff=5) {
  require(tidyverse); require(psych); require(broom)
  meta=data.frame(meta)
  sel=match(colnames(exp), rownames(meta))
  meta=meta[sel,,drop=F ]
  stopifnot(identical(colnames(exp), rownames(meta)))
  exp[is.na(exp)]=0
  pca <- 	prcomp(t(exp),rank. = 10, scale = FALSE)
  percentVar <- 	round((pca$sdev)^2/sum(pca$sdev^2), 3) * 100
  scores <- as.data.frame(pca$x)
  Npc=max(2,min(sum(percentVar>PC_cutoff), 10))
  
  cor_mat=NULL
  anova_results=NULL
  
  meta_num=select_if(meta, is.numeric)
  if (ncol(meta_num)>0) {
    cov_cor <- corr.test(scores[, 1:Npc, drop=F],
                         data.matrix(meta_num),
                         use = 'pairwise.complete.obs',
                         method = "kendall",
                         adjust = "none")
    all_cor_vals <- cov_cor[["r"]]
    all_cor_p <- cov_cor[["p"]]
    cor_mat <- reshape2::melt(all_cor_p, varnames = c("PC", "covar"))
    colnames(cor_mat)[colnames(cor_mat) == "value"] <- "pvalue"
    cor_mat$r <- reshape2::melt(all_cor_vals)[["value"]]
    cor_mat$fdr <- p.adjust(cor_mat$pvalue, method = "fdr") 
  }
  
  sel_col=which(!(colnames(meta) %in% colnames(meta_num)))
  meta_cat=meta[, sel_col, drop=F]
  if (ncol(meta_cat)>0) {
    for (i in 1:ncol(meta_cat)) {
      if (!is.factor(meta_cat[i, 1])) {meta_cat[, i]=as.factor(meta_cat[, i])} #turn all into factors, even logical
    }  
    nelem<-function(x) nlevels(factor(x))
    ( nL<-apply(meta_cat, 2, nelem) )
    meta_cat=meta_cat[, nL>1 & nL<nrow(meta_cat), drop=F]
    if (ncol(meta_cat)>0) {
      data.df=cbind(scores[, 1:Npc, drop=F], meta_cat)
      anova_results=NULL
      for (i in 1:ncol(meta_cat)) {
        for (j in 1:Npc) {
          test_formula=str_c(colnames(scores)[j], "~", colnames(meta_cat)[i])
          res.aov <- aov(as.formula(test_formula), data = data.df)
          if ("p.value" %in% names(tidy(res.aov))) {
            pvalue=tidy(res.aov)$p.value[1]
          } else {pvalue=1}
          PWinfo=NA
          if (pvalue<0.05) {  #run pairwise
            PW<-tidy(TukeyHSD(res.aov))
            PW.sig<-PW%>%filter(adj.p.value<0.05)
            if (nrow(PW.sig)>0) {
              names(PW.sig)[names(PW.sig)=="contrast"]="comparison"  #older version comparison, new version contrast
              PWinfo<-paste(paste(PW.sig$comparison, format.pval(PW.sig$adj.p.value, digits=2)), collapse="; ")
            }
          }
          result1=data.frame(PC=colnames(scores)[j], covar=colnames(meta_cat)[i], pvalue, PairWise=PWinfo)
          anova_results<-rbind(anova_results, result1)
        }
      }
      anova_results$fdr=p.adjust(anova_results$pvalue, method = "fdr") 
    }
  }
return(list(PC_scores=scores, meta=meta, percentVar=percentVar, corr=cor_mat, anova=anova_results))
}

#select significant covariate vs PC pairs, and create plot. 
#var_type can be "numeric" or "categorical"
#N_col: number of columns for facet plots
get_PC_meta_plot<-function(res, var_type, FDR_cutoff=0.1, N_col=3) {
  require(tidyverse);require(cowplot)
  if (var_type=="numeric") {
    if (is.null(res$corr)) {return(NULL)}
    selVar=res$corr%>%filter(fdr<FDR_cutoff)
  } else {
    if (is.null(res$anova)) {return(NULL)}
    selVar<-res$anova%>%filter(fdr<FDR_cutoff)
  }
  if (nrow(selVar)==0) {return(NULL)}
  meta=res$meta
  scores=res$PC_scores
  data.df=NULL
  for (i in 1:nrow(selVar)) {
    N1=match(selVar$PC[i], colnames(scores) )
    N2=match(selVar$covar[i], colnames(meta))
    df1=data.frame(PC=selVar$PC[i], covar=selVar$covar[i], Value=meta[, N2], Score=scores[, N1])
    data.df=rbind(data.df, df1)
  }
  
  ##make plot
  if (var_type=="numeric") {
    p<-ggplot(data.df, aes(x=Value, y=Score, color=covar) )+geom_point()+
      stat_summary(fun.data= mean_cl_normal) + geom_smooth(method='lm')+
      facet_wrap(c("covar", "PC"), ncol=N_col,scales = "free")+theme_half_open() +background_grid()+panel_border() + 
      labs(x="Covariate Value", y="PC Scores", color="numeric\ncovariates")  
    
    selVar$text=str_c("r=", round(selVar$r*1000)/1000,  "; fdr=",format.pval(selVar$fdr, digits=2))
    p<-p+ geom_text( data = selVar, color="black", 
                  mapping = aes(x = -Inf, y = Inf, label = text),
                  hjust   = -0.1, vjust   = 1.5)
    
  } else {
    p<-ggplot(data.df, aes(x=Value, y=Score, color=covar) )+geom_boxplot()+geom_jitter(alpha=0.7, width=0.1)+
    facet_wrap(c("covar", "PC"), ncol=N_col,scales = "free")+theme_half_open() +background_grid()+panel_border() + 
    theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5)) + labs(x="Covariate Category", y="PC Scores", color="categorical\ncovariates")  
    selVar$text=str_c("ANOVA fdr: ",format.pval(selVar$fdr, digits=2))
    p<-p+ geom_text( data = selVar, color="black", 
                mapping = aes(x = -Inf, y = Inf, label = text),
                hjust   = -0.1, vjust   = 1.5)
  }
  return(list(plot=p, data.df=data.df, selVar=selVar))
}
    
  


