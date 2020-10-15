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
	DataIn = DataReactive()
	#ProteinGeneName = DataIn$ProteinGeneName
	ProteinGeneName = DataIn$data_results
	DataIngenes <- ProteinGeneName %>% dplyr::select(UniqueID) %>% collect %>% .[["UniqueID"]] %>%	as.character()
	updateSelectizeInput(session,'sel_gene', choices= DataIngenes, server=TRUE)
})

observe({
	DataIn = DataReactive()
	groups = group_order()
	tests = DataIn$tests
	allgroups = DataIn$groups
	ProteinGeneName = DataIn$ProteinGeneName
	updateSelectizeInput(session,'sel_group', choices=allgroups, selected=groups)
	updateRadioButtons(session,'sel_geneid', inline = TRUE, choices=colnames(ProteinGeneName)[-1])
	updateSelectizeInput(session,'expression_test',choices=tests, selected=tests[1])
})


observe({
	DataIn = DataReactive()
	results_long = DataIn$results_long
	expression_test =input$expression_test
	expression_fccut = as.numeric(input$expression_fccut)
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

	validate(need(length(input$sel_group)>0,"Please select group(s)."))

	DataIn = DataReactive()
	data_long = DataIn$data_long
	results_long = DataIn$results_long
	ProteinGeneName = DataIn$ProteinGeneName
	sel_group=input$sel_group
	sel_gene=input$sel_gene
	genelabel=input$sel_geneid
	group_order(input$sel_group)

	if (input$exp_subset == "Select") {
	  validate(need(length(input$sel_gene)>0,"Please select a gene."))	
	  tmpids = ProteinGeneName[unique(na.omit(c(apply(ProteinGeneName,2,function(k) match(sel_gene,k))))),]
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
	
	data_long_tmp = filter(data_long, UniqueID %in% tmpids, group %in% sel_group) %>%
	filter(!is.na(expr)) %>% as.data.frame()
	data_long_tmp$labelgeneid = data_long_tmp[,match(genelabel,colnames(data_long_tmp))]
	data_long_tmp$group = factor(data_long_tmp$group,levels = sel_group)
	
	result_long_tmp = filter(results_long, UniqueID %in% tmpids) %>%  as.data.frame()

	return(list("data_long_tmp"=data_long_tmp,"result_long_tmp"= result_long_tmp, "tmpids"=tmpids))
	     
})

output$plot.exp=renderUI({
  req(DataExpReactive())
  graph_height=800
  if (input$SeparateOnePlot=="Separate") {
    graph_height=max(800, ceiling(length(DataExpReactive()$tmpids)/3)*300 )
    #cat("nrow: ", length(DataExpReactive()$tmpids), "\nheight: ", graph_height, "\n")
  }
  plotOutput("boxplot", height =graph_height)
})


output$dat_dotplot <- DT::renderDataTable({
	data_long_tmp <- DataExpReactive()$data_long_tmp
	data_long_tmp[,sapply(data_long_tmp,is.numeric)] <- signif(data_long_tmp[,sapply(data_long_tmp,is.numeric)],3)
	data_long_tmp <- data_long_tmp[,-7]
	DT::datatable(data_long_tmp,  extensions = 'Buttons',  options = list(
	  dom = 'lBfrtip', buttons = c('csv', 'excel', 'print'), pageLength = 15))
})

output$res_dotplot <- DT::renderDataTable({
  result_long_tmp <- DataExpReactive()$result_long_tmp
  result_long_tmp[,sapply(result_long_tmp,is.numeric)] <- signif(result_long_tmp[,sapply(result_long_tmp,is.numeric)],3)

  DT::datatable(result_long_tmp,  extensions = 'Buttons',  options = list(
    dom = 'lBfrtip', buttons = c('csv', 'excel', 'print'), pageLength = 15))

})

boxplot_out <- reactive({
  barcol = input$barcol
  sel_group=input$sel_group
  group_order(sel_group)
  
  data_long_tmp <- DataExpReactive()$data_long_tmp

  if (input$SeparateOnePlot == "Separate") {

    p <- ggplot(data_long_tmp,aes(x=group, y=expr, fill=factor(group))) +
      facet_wrap(~ labelgeneid, scales = "free", ncol = 3)
    if (input$plotformat == "boxplot") {
      p <- p + geom_boxplot() +
        stat_summary(fun=mean, geom="point", shape=18,size=4,color = "red")
    } else 	if (input$plotformat == "violin") {
      p <- p + geom_violin(trim = FALSE) +
        stat_summary(fun=mean, geom="point",shape=18,size=4,color = "red")
    } else if (input$plotformat == "barplot") {
      p <- p + stat_summary(fun.data=mean_se, position=position_dodge(0.8), geom="errorbar",aes(width=0.5)) +
        stat_summary(fun=mean, position=position_dodge(0.8), geom="bar")
    } else if (input$plotformat == "line") {
      p <- p + stat_summary(fun=mean, geom="point",shape=18, size=4, color = "red") +
        stat_summary(aes(y = expr, group=1), fun=mean, colour="red", geom="line", group=1)+
        stat_summary(fun.data=mean_se, position=position_dodge(0.8), geom="errorbar",aes(width=0.5))
    }
    if (input$IndividualPoint == "YES")
      p <- p +  geom_dotplot(binaxis='y', stackdir='center', position = "dodge", dotsize = 0.5)
    
    if (input$ColPattern == "Palette") {
      p <- p + scale_fill_manual(values = colorRampPalette(brewer.pal(8, input$colpalette))(length(sel_group)))
    } else {
      p <- p + scale_fill_manual(values=rep(barcol, length(sel_group)))
    } 

    p <- p + theme_bw(base_size = 14) + ylab(input$Ylab) + xlab(input$Xlab) +
      theme (plot.margin = unit(c(1,1,1,1), "cm"),
             text = element_text(size=input$expression_axisfontsize),
             axis.text.x = element_text(angle = input$Xangle, hjust=0.5, vjust=0.5),
             legend.position="none", 
             strip.text.x = element_text(size=input$expression_titlefontsize))
    
  }
 # browser()
  if (input$SeparateOnePlot == "OnePlot") {
    
    library(plyr)
    
    data_long_tmp1 <- ddply(data_long_tmp, c("UniqueID", "group"), summarise,
                           N    = sum(!is.na(expr)),
                           mean = mean(expr, na.rm=TRUE),
                           sd   = sd(expr, na.rm=TRUE),
                           se   = sd / sqrt(N)
    )
    
    data_long_tmp1 <- data_long_tmp1 %>%left_join(data_long_tmp%>%filter(!duplicated(UniqueID))%>%transmute(UniqueID, Gene.Name=labelgeneid)   )


    pd <- position_dodge(0.1) # move them .05 to the left and right
    p <-	ggplot(data_long_tmp1, aes(x=group, y=mean, group=Gene.Name))  
    
    if (input$plotformat == "line") {
      p <- p + geom_errorbar(aes(ymin=mean-se, ymax=mean+se, color = Gene.Name),size=1, width=.2, position=pd) +
        geom_line(position=pd, size = 1, aes(color = Gene.Name)) +
        geom_point(position=pd, size=3, shape=21, fill="white")
    } else {
      p <- p + geom_bar(aes(fill= Gene.Name), position=position_dodge(), stat="identity", colour="black", size=.3) + 
        geom_errorbar(aes(ymin=mean-se, ymax=mean+se), size=.3, width=.2, position=position_dodge(.9))
    }
    
    p <- p + theme_bw(base_size = 14) + ylab(input$Ylab) + xlab(input$Xlab) +
      theme (plot.margin = unit(c(1,1,1,1), "cm"),
             text = element_text(size=input$expression_axisfontsize),
             axis.text.x = element_text(angle = input$Xangle, hjust=0.5, vjust=0.5),
             strip.text.x = element_text(size=input$expression_titlefontsize))
  }

  p
  
  
})

output$boxplot <- renderPlot({
boxplot_out()
})

observeEvent(input$boxplot, {
	saved.num <- length(saved_plots$boxplot) + 1
	saved_plots$boxplot[[saved.num]] <- boxplot_out()
})

browsing_out <- reactive({
	validate(need(length(input$sel_group)>0,"Please select group(s)."))
	barcol = input$barcol
	DataIn = DataReactive()
	data_long = DataIn$data_long
	results_long = DataIn$results_long
	ProteinGeneName = DataIn$ProteinGeneName

	genelabel=input$sel_geneid
	sel_group=input$sel_group
	group_order(sel_group)
	expression_test = input$expression_test
	expression_fccut =as.numeric(input$expression_fccut)
	expression_pvalcut = as.numeric(input$expression_pvalcut)
	numperpage = as.numeric(input$numperpage)

	sel_page = as.numeric(input$sel_page)-1
	startslice = sel_page * 6 + 1
	endslice = startslice + numperpage -1
	if (input$expression_psel == "Padj") {
		sel_gene = results_long %>% filter(test %in% expression_test & abs(logFC) > expression_fccut & Adj.P.Value < expression_pvalcut) %>%
		dplyr::arrange(P.Value) %>%
		dplyr::slice(startslice:endslice) %>%
		dplyr::select(UniqueID) %>%
		collect %>% .[["UniqueID"]] %>% as.character()
	} else {
		sel_gene = results_long %>% filter(test %in% expression_test & abs(logFC) > expression_fccut & P.Value < expression_pvalcut) %>%
		dplyr::arrange(P.Value) %>%
		dplyr::slice(startslice:endslice) %>%
		dplyr::select(UniqueID) %>%
		collect %>% .[["UniqueID"]] %>% as.character()
	}

	tmpids = ProteinGeneName[unique(na.omit(c(apply(ProteinGeneName,2,function(k) match(sel_gene,k))))),]

	data_long_tmp = filter(data_long, UniqueID %in% tmpids$UniqueID, group %in% sel_group) %>%
	filter(!is.na(expr)) %>% as.data.frame()

	data_long_tmp$labelgeneid = data_long_tmp[,match(genelabel,colnames(data_long_tmp))]
	data_long_tmp$group = factor(data_long_tmp$group,levels = sel_group)
	validate(need(nrow(data_long_tmp)>0, message = "Please select at least one valid gene to plot."))
  #browser() #debug
	if(numperpage==4) {
		p <- ggplot(data_long_tmp,aes(x=group,y=expr,fill=group)) +
		facet_wrap(~ labelgeneid, scales = "free",nrow = 2, ncol = 2)
	}
	if(numperpage==6) {
		p <- ggplot(data_long_tmp,aes(x=group,y=expr,fill=group)) +
		facet_wrap(~ labelgeneid, scales = "free", nrow = 2, ncol = 3)
	}
	if(numperpage==9) {
		p <- ggplot(data_long_tmp,aes(x=group,y=expr,fill=group)) +
		facet_wrap(~ labelgeneid, scales = "free", nrow = 3, ncol = 3)
	}

	if (input$plotformat == "boxplot") {
		p <- p + geom_boxplot() +
		stat_summary(fun=mean, geom="point", shape=18,size=4,color = "red")
	}
	if (input$plotformat == "violin") {
		p <- p + geom_violin(trim = FALSE) +
		stat_summary(fun=mean, geom="point",shape=18,size=4,color = "red")
	}
	if (input$plotformat == "barplot") {
		p <- p + stat_summary(fun.data=mean_se, position=position_dodge(0.8), geom="errorbar",aes(width=0.5)) +
		stat_summary(fun=mean, position=position_dodge(0.8), geom="bar")
	}
	if (input$plotformat == "line") {
		p <- p + stat_summary(fun=mean, geom="point",shape=18, size=4, color = "red") +
		stat_summary(aes(y = expr, group=1), fun=mean, colour="red", geom="line", group=1)+
		stat_summary(fun.data=mean_se, position=position_dodge(0.8), geom="errorbar",aes(width=0.5))
	}

	if (input$IndividualPoint == "YES")
	p <- p +  geom_dotplot(binaxis='y', stackdir='center', dotsize = 0.5)

	if (input$ColPattern == "Palette") {
			#p <- p + scale_fill_brewer(palette=input$colpalette)
			p <- p + scale_fill_manual(values = colorRampPalette(brewer.pal(8, input$colpalette))(length(sel_group)))
	} else {
		p <- p + scale_fill_manual(values=rep(barcol,length(sel_group)))
	}

	p <- p +	theme_bw(base_size = 14) + ylab(input$Ylab) + xlab(input$Xlab) +
	  theme (plot.margin = unit(c(1,1,1,1), "cm"),
	         text = element_text(size=input$expression_axisfontsize),
	         axis.text.x = element_text(angle = input$Xangle, hjust=0.5, vjust=0.5),
	         legend.position="none", 
	         strip.text.x = element_text(size=input$expression_titlefontsize))
	p

})

output$browsing <- renderPlot({
	browsing_out()
})

observeEvent(input$browsing, {
	saved.num <- length(saved_plots$browsing) +1
	saved_plots$browsing[[saved.num]] <- browsing_out()
})


