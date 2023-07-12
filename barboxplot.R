###########################################################################################################
## Proteomics Visualization R Shiny App
##
##This software belongs to Biogen Inc. All right reserved.
##
##@file: barplot.R
##@Developer : Benbo Gao (benbo.gao@Biogen.com)
##@Date : 5/16/2018
##@version 1.0
###########################################################################################################


observe({
  MetaData=all_metadata()
	DataIn = DataReactive()
	ProteinGeneName = DataIn$ProteinGeneName
	#ProteinGeneName = DataIn$data_results
	#DataIngenes <- ProteinGeneName %>% dplyr::select(UniqueID) %>% collect %>% .[["UniqueID"]] %>%	as.character()
	if (input$exp_label=="UniqueID") {
	  DataIngenes <- ProteinGeneName %>% dplyr::select(UniqueID) %>% collect %>% .[["UniqueID"]] %>%	as.character()
	} else
	{DataIngenes <- ProteinGeneName %>% dplyr::select(Gene.Name) %>% collect %>% .[["Gene.Name"]] %>%	as.character()}
	updateSelectizeInput(session,'sel_gene', choices= DataIngenes, server=TRUE)
	attributes=sort(setdiff(colnames(MetaData), c("sampleid", "Order", "ComparePairs") ))
	updateSelectInput(session, "colorby", choices=c("None", attributes), selected="group")
	updateSelectInput(session, "plotx", choices=attributes, selected="group")

})

observe({
	#DataIn = DataReactive()
	#groups = group_order()
	tests = all_tests()
	allgroups = all_groups()
	ProteinGeneName_Header = ProteinGeneNameHeader()
	#updateSelectizeInput(session,'sel_group', choices=allgroups, selected=groups)
	updateRadioButtons(session,'sel_geneid', inline = TRUE, choices=c(ProteinGeneName_Header[-1], "Gene.Name_UniqueID"), selected="Gene.Name")
	updateSelectizeInput(session,'expression_test',choices=tests, selected=tests[1])
})

#no longer needed, assign at heatmap.R
#output$selectGroupSampleExpression <- renderText({ paste("Selected ",length(group_order()), " out of ", length(all_groups()), " Groups, ",
#     length(sample_order()), " out of ", length(all_samples()), " Samples.", " (Update Selection at: QC Plots->Groups and Samples.)", sep="")})

#linear value parameters
observe( {
  if (input$exp_plot_Y_scale=='Linear') {
    expU = exp_unit()
    N_log = as.numeric(str_replace(str_extract(expU,"log\\d+"),"log",""))
    small_value=as.numeric(str_replace(str_split_fixed(expU, "\\+", 2)[2], "\\)", ""))
    unit=str_replace_all(str_extract(expU, "\\(.+\\+"), "(\\(|\\+)", "")
    if (!is.na(N_log) & !is.na(small_value)) {
      updateTextInput(session, "linear_base", value=N_log)
      updateTextInput(session, "linear_small_value", value=small_value)
      updateTextInput(session, "Ylab", value=unit)
    }
  } else if (input$exp_plot_Y_scale=='Log') {updateTextInput(session, "Ylab", value=exp_unit())}
})

observe({
	DataIn = DataReactive()
	results_long = DataIn$results_long
	expression_test =input$expression_test
	expression_fccut = log2(as.numeric(input$expression_fccut))
	expression_pvalcut =  as.numeric(input$expression_pvalcut)
	numperpage = as.numeric(input$numperpage)

	if (input$expression_psel == "Padj") {
		filteredgene = results_long %>%
		dplyr::filter(abs(logFC) > expression_fccut & Adj.P.Value < expression_pvalcut) %>%
		dplyr::filter(test == expression_test)
	} else {
		filteredgene = results_long %>%
		dplyr::filter(abs(logFC) > expression_fccut & P.Value < expression_pvalcut) %>%
		dplyr::filter(test == expression_test)
	}

	output$expfilteredgene <- renderText({paste("Selected Genes:",nrow(filteredgene),sep="")})
	updateSelectInput(session,'sel_page', choices= seq_len(ceiling(nrow(filteredgene)/numperpage)))
})


DataExpReactive <- reactive({

	validate(need(length(group_order())>0,"Please select group(s)."))

	DataIn = DataReactive()
	data_long = DataIn$data_long
	results_long = DataIn$results_long
	ProteinGeneName = DataIn$ProteinGeneName
	sel_group=group_order() #input$sel_group
	sel_gene=input$sel_gene
	genelabel=input$sel_geneid
	#group_order(input$sel_group)
	sel_samples=sample_order()

	if (input$exp_subset == "Select") {
	  validate(need(length(input$sel_gene)>0,"Please select a gene."))
	  if (input$exp_label=="UniqueID") {
	    tmpids = ProteinGeneName[unique(na.omit(c(apply(ProteinGeneName,2,function(k) match(sel_gene,k))))),]
	  } else { #Gene.Name can be duplicate
	    tmpids <- ProteinGeneName %>% dplyr::filter (Gene.Name %in% sel_gene)
	  }
	  tmpids=tmpids$UniqueID
	}
	if (input$exp_subset == "Upload Genes") {
	    exp_list <- input$exp_list
	    if(grepl("\n",exp_list)) {
	      exp_list <-  stringr::str_split(exp_list, "\n")[[1]]
	    } else if(grepl(",",exp_list)) {
	      exp_list <-  stringr::str_split(exp_list, ",")[[1]]
	    }
	    exp_list <- gsub(" ", "", exp_list, fixed = TRUE)
	    exp_list <- unique(exp_list[exp_list != ""])
	    validate(need(length(exp_list)>0, message = "Please input at least 1 valid genes."))
	    tmpids <- dplyr::filter(ProteinGeneName, (UniqueID %in% exp_list) | (Protein.ID %in% exp_list) | (toupper(Gene.Name) %in% toupper(exp_list)))  %>%
	      dplyr::select(UniqueID) %>% 	collect %>%	.[["UniqueID"]] %>%	as.character()
	    validate(need(length(tmpids)>0, message = "Please input at least 1 valid genes."))
	}
	if (input$exp_subset == "Geneset") {
	  req(input$geneset_list_exp)
	  exp_list <- input$geneset_list_exp
	  if(grepl("\n",exp_list)) {
	    exp_list <-  stringr::str_split(exp_list, "\n")[[1]]
	  } else if(grepl(",",exp_list)) {
	    exp_list <-  stringr::str_split(exp_list, ",")[[1]]
	  }
	  exp_list <- gsub(" ", "", exp_list, fixed = TRUE)
	  exp_list <- unique(exp_list[exp_list != ""])
	  tmpids <- dplyr::filter(ProteinGeneName, (UniqueID %in% exp_list) | (Protein.ID %in% exp_list) | (toupper(Gene.Name) %in% toupper(exp_list)))  %>%
	    dplyr::select(UniqueID) %>% 	collect %>%	.[["UniqueID"]] %>%	as.character()
	  validate(need(length(tmpids)>0, message = "Please input at least 1 valid genes."))
	}
	if (length(tmpids)>100) {cat("show only first 100 genes in exprssion plot.\n"); tmpids=tmpids[1:100]}

	data_long_tmp = filter(data_long, UniqueID %in% tmpids, group %in% sel_group, sampleid %in% sel_samples) %>%
	filter(!is.na(expr)) %>% as.data.frame()
	data_long_tmp<-data_long_tmp%>%mutate(Gene.Name_UniqueID=str_c(Gene.Name, "_", UniqueID))
	data_long_tmp$labelgeneid = data_long_tmp[,match(genelabel,colnames(data_long_tmp))]
	data_long_tmp$group = factor(data_long_tmp$group,levels = sel_group)
	if (input$exp_plot_Y_scale=='Linear') {
	  data_long_tmp<-data_long_tmp%>%mutate(expr=input$linear_base^(expr-input$linear_small_value))
	}
	result_long_tmp = filter(results_long, UniqueID %in% tmpids) %>%  as.data.frame()
	Ng=length(unique(data_long_tmp$Gene.Name)); Nuid=length(unique(data_long_tmp$UniqueID))
	search_gene_info<-str_c("Displaying ", Ng, " Gene.Names from ", Nuid, " UniqueIDs.")
	gene_multi_uid<-result_long_tmp%>%distinct(Gene.Name, UniqueID)%>%group_by(Gene.Name)%>%dplyr::count()%>%dplyr::filter(n>1)
	if (nrow(gene_multi_uid)>0) {search_gene_info<-str_c(search_gene_info, 
	       "\nPlease note some gene names map to mulitiple UniqueIDs, we recommend using Gene.Name_UniqueID as Gene Label to separate the UniqueIDs in the plot.")}
	output$geneSearchInfo<-renderText({search_gene_info})
	#browser() #debug

	return(list("data_long_tmp"=data_long_tmp,"result_long_tmp"= result_long_tmp, "tmpids"=tmpids))

})

output$dat_dotplot <- DT::renderDT(server=FALSE, {
	data_long_tmp <- DataExpReactive()$data_long_tmp
	data_long_tmp <- data_long_tmp %>%dplyr::select(-labelgeneid, -Gene.Name_UniqueID)
	data_long_tmp[,sapply(data_long_tmp,is.numeric)] <- signif(data_long_tmp[,sapply(data_long_tmp,is.numeric)],3)
	#data_long_tmp <- data_long_tmp[,-7]
	DT::datatable(data_long_tmp,  extensions = 'Buttons',  options = list(
	  dom = 'lBfrtip', pageLength = 15,
	  buttons = list(
	    list(extend = "csv", text = "Download Page", filename = "Page_results",
	         exportOptions = list(modifier = list(page = "current"))),
	    list(extend = "csv", text = "Download All", filename = "All_Results",
	         exportOptions = list(modifier = list(page = "all")))
	  )
	  ))
})

output$res_dotplot <- DT::renderDT(server=FALSE,{
  result_long_tmp <- DataExpReactive()$result_long_tmp
  result_long_tmp[,sapply(result_long_tmp,is.numeric)] <- signif(result_long_tmp[,sapply(result_long_tmp,is.numeric)],3)
  DT::datatable(result_long_tmp,  extensions = 'Buttons',  options = list(
    dom = 'lBfrtip', pageLength = 15,
    buttons = list(
      list(extend = "csv", text = "Download Page", filename = "Page_results",
           exportOptions = list(modifier = list(page = "current"))),
      list(extend = "csv", text = "Download All", filename = "All_Results",
           exportOptions = list(modifier = list(page = "all")))
    )
  ))
})

boxplot_out <- eventReactive(input$plot_exp,  {
  barcol = input$barcol
  sel_group=group_order() #input$sel_group
  #group_order(sel_group)
  DataIn = DataReactive()
  colorby=sym(input$colorby)
  Val_colorby=input$colorby
  MetaData=DataIn$MetaData
  plotx=sym(input$plotx)
  ncol=input$exp_plot_ncol
  data_long_tmp <- DataExpReactive()$data_long_tmp
  if (Val_colorby!="None" & Val_colorby!="group" ) { #add coloyby column
    data_long_tmp<-data_long_tmp%>%left_join(MetaData%>%dplyr::select(sampleid, !!colorby))
  } else {
    data_long_tmp$None="None"
  }
  if (input$plotx!="group" ) { #add plotx column
    data_long_tmp<-data_long_tmp%>%left_join(MetaData%>%dplyr::select(sampleid, !!plotx))
  }

  if (input$SeparateOnePlot == "Separate") {

    p <- ggplot(data_long_tmp,aes(x=!!plotx,y=expr,fill=!!colorby)) +
      facet_wrap(~ labelgeneid, scales = "free", ncol = ncol)
    if (input$plotformat == "boxplot") {
      p <- p + geom_boxplot() +
        stat_summary(aes(group=!!colorby), fun=mean, geom="point", shape=18,size=3, color = "red", position = position_dodge(width=0.8))
    }
    if (input$plotformat == "violin") {
      p <- p + geom_violin(trim = FALSE) +
        stat_summary(fun=mean, geom="point",shape=18,size=3,color = "red",position = position_dodge(width=0.8))
    }
    if (input$plotformat == "barplot") {
      p <- p + stat_summary(fun.data=mean_se, position=position_dodge(0.8), geom="errorbar",aes(width=0.5)) +
        stat_summary(fun=mean, position=position_dodge(0.8), geom="bar")
    }
    if (input$plotformat == "line") {
      p <- p + stat_summary(aes(color=!!colorby), fun=mean, geom="point",shape=18, size=3) +
        stat_summary(aes(y = expr, group=!!colorby, color=!!colorby), fun=mean, geom="line")+
        stat_summary(fun.data=mean_se, geom="errorbar",aes(width=0.3, color=!!colorby))
    }

    if (input$IndividualPoint == "YES")
      #browser()
      p <- p +   geom_jitter(aes(fill=!!colorby), shape=21, size=2, color="black", position = position_jitterdodge(jitter.width=0.25))
        #geom_dotplot(binaxis='y', stackdir='center', dotsize = 0.5,  position = position_dodge(width=0.8))
    if (Val_colorby!="None" ) {
      #browser()
      N_color<-data_long_tmp%>%dplyr::select(!!colorby)%>%unlist%>%unname%>%as.character%>%unique%>%length
      use_color=colorRampPalette(brewer.pal(8, input$colpalette))(N_color)
      if (input$plotformat == "line") {
        p <- p +scale_color_manual(values=use_color)+ scale_fill_manual(values =use_color)
      } else {p <- p + scale_fill_manual(values =use_color)}

    } else {
      p <- p + scale_fill_manual(values=barcol) #+scale_color_manual(values=rep(barcol,length(sel_group)))
    }

    p <- p + theme_bw(base_size = 14) + ylab(input$Ylab) + xlab(input$Xlab) +guides(fill = guide_legend(override.aes = list(shape = NA) ) )+
      theme (plot.margin = unit(c(1,1,1,1), "cm"),
             text = element_text(size=input$expression_axisfontsize),
             axis.text.x = element_text(angle = input$Xangle, hjust=0.5, vjust=0.5),
             strip.text.x = element_text(size=input$expression_titlefontsize))
    if (Val_colorby=="None" ) {p <- p +	 theme (legend.position="none") }
  }
 # browser()
  if (input$SeparateOnePlot == "OnePlot") {
    data_long_tmp1 <- ddply(data_long_tmp, c("UniqueID", input$plotx), summarise,
                           N    = sum(!is.na(expr)),
                           mean = mean(expr, na.rm=TRUE),
                           sd   = sd(expr, na.rm=TRUE),
                           se   = sd / sqrt(N)
    )

    data_long_tmp1 <- data_long_tmp1 %>%left_join(data_long_tmp%>%filter(!duplicated(UniqueID))%>%transmute(UniqueID, Gene.Name=labelgeneid)   )


    pd <- position_dodge(0.1) # move them .05 to the left and right
    p <-	ggplot(data_long_tmp1, aes(x=!!plotx, y=mean, group=Gene.Name))

    if (input$plotformat == "line") {
      p <- p + geom_errorbar(aes(ymin=mean-se, ymax=mean+se, color = Gene.Name),size=1, width=.2, position=pd) +
        geom_line(position=pd, size = 1, aes(color = Gene.Name)) +
        geom_point(position=pd, size=3, shape=21, fill="white")
    } else {
      p <- p + geom_bar(aes(fill= Gene.Name), position=position_dodge(), stat="identity", colour="black", size=.3) +
        geom_errorbar(aes(ymin=mean-se, ymax=mean+se), size=.3, width=.2, position=position_dodge(.9))
    }

    p <- p + theme_bw(base_size = 14) + ylab(input$Ylab) + xlab(input$Xlab) +scale_fill_discrete(name=input$sel_geneid)+
      theme (plot.margin = unit(c(1,1,1,1), "cm"),
             text = element_text(size=input$expression_axisfontsize),
             axis.text.x = element_text(angle = input$Xangle, hjust=0.5, vjust=0.5),
             strip.text.x = element_text(size=input$expression_titlefontsize))
  }
  if (input$exp_plot_Y_range=="Manual") {
    p <- p + ylim(input$exp_plot_Ymin, input$exp_plot_Ymax)
  }
  p
})

graph_height_boxplot=eventReactive(input$plot_exp,  {
  D_exp<-DataExpReactive()
  graph_height=800
  if (input$SeparateOnePlot=="Separate") {
    graph_height=max(800, ceiling(length(D_exp$tmpids)/3)*300 )
  }
  return(graph_height)
})


output$plot.exp=renderUI({
  graph_height=graph_height_boxplot()
  if (is.null(graph_height)){graph_height=800}
  plotOutput("boxplot", height =graph_height)
})

output$boxplot <- renderPlot({
  withProgress(message = 'Making Expression Plot of selected genes...', value = 0, {
    p_boxplot=boxplot_out()
    print(p_boxplot)
  })
})




observeEvent(input$boxplot, {
	saved.num <- length(saved_plots$boxplot) + 1
	saved_plots$boxplot[[saved.num]] <- boxplot_out()
})

observeEvent(input$plot_browsing, {
  plot_exp_control(plot_exp_control()+1)
})

browsing_out <- eventReactive(plot_exp_control(),{
	validate(need(length(group_order())>0,"Please select group(s)."))
	barcol = input$barcol
	DataIn = DataReactive()
	data_long = DataIn$data_long
	results_long = DataIn$results_long
	ProteinGeneName = DataIn$ProteinGeneName
	colorby=sym(input$colorby)
	Val_colorby=input$colorby
	MetaData=DataIn$MetaData
	plotx=sym(input$plotx)
	genelabel=input$sel_geneid
	sel_group=group_order() #input$sel_group
	sel_samples=sample_order()
	#group_order(sel_group)
	expression_test = input$expression_test
	expression_fccut =log2(as.numeric(input$expression_fccut))
	expression_pvalcut = as.numeric(input$expression_pvalcut)
	numperpage = as.numeric(input$numperpage)
	ncol=input$exp_plot_ncol
	nrow=ceiling(numperpage/ncol)

	sel_page = as.numeric(input$sel_page)-1
	startslice = sel_page * 6 + 1
	endslice = startslice + numperpage -1
	if (input$browsing_gene_order=="P value") {results_long<-results_long%>%arrange(P.Value)
	} else {results_long<-results_long%>%arrange(desc(abs(logFC)))}

	if (input$expression_psel == "Padj") {
		sel_gene = results_long %>% filter(test %in% expression_test & abs(logFC) > expression_fccut & Adj.P.Value < expression_pvalcut) %>%
		dplyr::slice(startslice:endslice) %>%
		dplyr::select(UniqueID) %>%
		collect %>% .[["UniqueID"]] %>% as.character()
	} else {
		sel_gene = results_long %>% filter(test %in% expression_test & abs(logFC) > expression_fccut & P.Value < expression_pvalcut) %>%
		dplyr::slice(startslice:endslice) %>%
		dplyr::select(UniqueID) %>%
		collect %>% .[["UniqueID"]] %>% as.character()
	}

	tmpids = ProteinGeneName[unique(na.omit(c(apply(ProteinGeneName,2,function(k) match(sel_gene,k))))),]

	data_long_tmp = filter(data_long, UniqueID %in% tmpids$UniqueID, group %in% sel_group, sampleid %in% sel_samples) %>%
	filter(!is.na(expr)) %>% as.data.frame()
	if (Val_colorby!="None" & Val_colorby!="group" ) { #add coloyby column
	  data_long_tmp<-data_long_tmp%>%left_join(MetaData%>%dplyr::select(sampleid, !!colorby))
	} else {
	  data_long_tmp$None="None"
	}
	if (input$plotx!="group" ) { #add plotx column
	  data_long_tmp<-data_long_tmp%>%left_join(MetaData%>%dplyr::select(sampleid, !!plotx))
	}
#	browser() #debug
	data_long_tmp<-data_long_tmp%>%mutate(Gene.Name_UniqueID=str_c(Gene.Name, "_", UniqueID))
	data_long_tmp$labelgeneid = data_long_tmp[,match(genelabel,colnames(data_long_tmp))]
	data_long_tmp$group = factor(data_long_tmp$group,levels = sel_group)
	validate(need(nrow(data_long_tmp)>0, message = "Please select at least one valid gene to plot."))
	#browser() #debug
	data_long_tmp$labelgeneid=factor(data_long_tmp$labelgeneid, levels=unique(data_long_tmp$labelgeneid))
	if (input$exp_plot_Y_scale=='Linear') {
	  data_long_tmp<-data_long_tmp%>%mutate(expr=input$linear_base^(expr-input$linear_small_value))
	}

	p <- ggplot(data_long_tmp,aes(x=!!plotx,y=expr,fill=!!colorby)) +
	  facet_wrap(~ labelgeneid, scales = "free",nrow = nrow, ncol = ncol)

	if (input$plotformat == "boxplot") {
		p <- p + geom_boxplot() +
		stat_summary(aes(group=!!colorby), fun=mean, geom="point", shape=18,size=3, color = "red", position = position_dodge(width=0.8))
	}
	if (input$plotformat == "violin") {
		p <- p + geom_violin(trim = FALSE) +
		stat_summary(fun=mean, geom="point",shape=18,size=3,color = "red",position = position_dodge(width=0.8))
	}
	if (input$plotformat == "barplot") {
		p <- p + stat_summary(fun.data=mean_se, position=position_dodge(0.8), geom="errorbar",aes(width=0.5)) +
		stat_summary(fun=mean, position=position_dodge(0.8), geom="bar")
	}
	if (input$plotformat == "line") {
		p <- p + stat_summary(aes(color=!!colorby), fun=mean, geom="point",shape=18, size=3) +
		stat_summary(aes(y = expr, group=!!colorby, color=!!colorby), fun=mean, geom="line")+
		stat_summary(fun.data=mean_se, geom="errorbar",aes(width=0.3, color=!!colorby))
	}

	if (input$IndividualPoint == "YES")
	  p <- p +   geom_jitter(aes(fill=!!colorby), shape=21, size=2, color="black", position = position_jitterdodge(jitter.width=0.25))
	#geom_dotplot(binaxis='y', stackdir='center', dotsize = 0.5,  position = position_dodge(width=0.8))
	if (Val_colorby!="None" ) {
	    #browser()
	    N_color<-data_long_tmp%>%dplyr::select(!!colorby)%>%unlist%>%unname%>%as.character%>%unique%>%length
	    use_color=colorRampPalette(brewer.pal(8, input$colpalette))(N_color)
			if (input$plotformat == "line") {
			  p <- p +scale_color_manual(values=use_color)+ scale_fill_manual(values =use_color)
			} else {p <- p + scale_fill_manual(values =use_color)}

	} else {
		p <- p + scale_fill_manual(values=barcol) #+scale_color_manual(values=rep(barcol,length(sel_group)))
	}
	p <- p +	theme_bw(base_size = 14) + ylab(input$Ylab) + xlab(input$Xlab) +guides(fill = guide_legend(override.aes = list(shape = NA) ) )+
	  theme (plot.margin = unit(c(1,1,1,1), "cm"),
	         text = element_text(size=input$expression_axisfontsize),
	         axis.text.x = element_text(angle = input$Xangle, hjust=0.5, vjust=0.5),
	         strip.text.x = element_text(size=input$expression_titlefontsize))
	if (Val_colorby=="None" ) {
	p <- p +	 theme (legend.position="none")
	}

	if (input$exp_plot_Y_range=="Manual") {
	  p <- p + ylim(input$exp_plot_Ymin, input$exp_plot_Ymax)
	}
	p

})

output$browsing <- renderPlot({
  ptm <- proc.time()
  withProgress(message = 'Drawing Expression Plot...\nIt may take a while', value = 0, {
  print(browsing_out()) })
  cat("plotted expression plot",(proc.time() - ptm)[["elapsed"]], "\n")
})

observeEvent(input$browsing, {
	saved.num <- length(saved_plots$browsing) +1
	saved_plots$browsing[[saved.num]] <- browsing_out()
})
