###########################################################################################################
## Proteomics Visualization R Shiny App
##
##This software belongs to Biogen Inc. All right reserved.
##
##@file: venn.R
##@Developer : Benbo Gao (benbo.gao@Biogen.com)
##@Date : 5/16/2018
##@version 1.0
###########################################################################################################


observe({
	DataIn = DataReactive()
	tmptests = DataIn$tests
	ntest <- length(tmptests)
	if (ntest >= 5) 	{
		ntest = 5
		tmptests = c(tmptests, "Empty List")
	}
	if (ntest < 5) 	{
		emptylist = 5- ntest
		tmptests = c(tmptests, rep("Empty List", emptylist))
	}
	for (i in 1:5){
		venn_test <- paste("venn_test",i,sep="")
		updateSelectizeInput(session, venn_test, choices=tmptests, selected=tmptests[i])
	}
})

DataVennReactive <- reactive({
	DataIn = DataReactive()
	results_long = DataIn$results_long
	venn_fccut = log2(as.numeric(input$venn_fccut))
	venn_pvalcut = as.numeric(input$venn_pvalcut)
	
	if (input$venn_psel == "Padj") {
		results_long <-  results_long %>% dplyr::filter(abs(logFC) > venn_fccut & Adj.P.Value < venn_pvalcut)
	} else{
		results_long <-  results_long %>% dplyr::filter(abs(logFC) > venn_fccut & P.Value < venn_pvalcut)
	}
	
	if (input$venn_updown == "Up") {
	  results_long <-  results_long %>% dplyr::filter(logFC > 0)
	} 
	
	if (input$venn_updown == "Down") {
	  results_long <-  results_long %>% dplyr::filter(logFC < 0)
	} 
	
	vennlist <- list()
	fill <- list() 
		
	if (input$venn_test1 != "Empty List") {
		list1 = results_long %>% 
		filter(test == input$venn_test1) %>%
		dplyr::select(UniqueID) %>%	collect %>%	.[["UniqueID"]] %>% as.character()
		fill[[input$venn_test1]] <- input$col1
		vennlist[[input$venn_test1]]  <- list1
	}

	if (input$venn_test2 != "Empty List") {
		list2 = results_long %>% 
		filter(test == input$venn_test2) %>%
		dplyr::select(UniqueID) %>%	collect %>%	.[["UniqueID"]] %>% as.character()
		fill[[input$venn_test2]] <-input$col2
		vennlist[[input$venn_test2]]  <- list2
	}

	if (input$venn_test3 != "Empty List") {
		list3 = results_long %>%
		filter(test == input$venn_test3) %>%
		dplyr::select(UniqueID) %>%	collect %>%	.[["UniqueID"]] %>% as.character()
		fill[[input$venn_test3]] <-input$col3
		vennlist[[input$venn_test3]]  <- list3
	}

	if (input$venn_test4 != "Empty List") {
		list4 = results_long %>% 
		filter(test == input$venn_test4) %>%
		dplyr::select(UniqueID) %>%	collect %>%	.[["UniqueID"]] %>% as.character()
		fill[[input$venn_test4]] <-input$col4
		vennlist[[input$venn_test4]]  <- list4
	}

	if (input$venn_test5 != "Empty List") {
		list5 = results_long %>% 
		filter(test == input$venn_test5) %>%
		dplyr::select(UniqueID) %>%	collect %>%	.[["UniqueID"]] %>% as.character()
		fill[[input$venn_test5]] <-input$col5
		vennlist[[input$venn_test5]]  <- list5
	}

	return(venndata = list("vennlist"=vennlist, "fillcor"=fill) )
})

vennDiagram_out <- reactive({
	print("drawing Venn diagram")
	venndata <- DataVennReactive()
	vennlist <- venndata$vennlist
	fillcor <- unlist(venndata$fillcor)
	SetNum = length(vennlist)
	futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
	venn.plot <- venn.diagram(x = vennlist,
		fill=fillcor, margin=input$margin,
		lty=input$lty, lwd=input$lwd, alpha=input$alpha,
		cex=input$cex, cat.cex=input$catcex,
		fontface = input$fontface, cat.fontface=input$catfontface,
		main = input$title, main.cex = input$maincex, main.pos = c(0.5, 1.1), main.fontface = "bold",
	filename = NULL)
	
	return(venn.plot)
})

output$vennDiagram <- renderPlot({
	grid.draw(vennDiagram_out())
})

#show all DEGs from selected comparisons
output$venn_DEG_Data <- DT::renderDataTable({
  venndata <- DataVennReactive()
  vennlist <- venndata$vennlist
  allIDs=unique(unlist(vennlist))
  dataIn=DataReactive()
  data_results=dataIn$data_results
  all_names=names(data_results)
  tests=names(vennlist)
  selCol=NULL
  for (i in 1:length(tests)) {
    sel_i=which(str_detect(all_names, regex(str_c("^", tests[i]), ignore_case=T)))
    if (length(sel_i)>0) {selCol=c(selCol, sel_i)}
  }
  name_col=which(all_names %in% c("UniqueID", "Gene.Name") )
  sel_row=which(data_results$UniqueID %in% allIDs)
  #browser()#debug
  DEG_outdata=data_results[sel_row, c(name_col, selCol)]
  DEG_outdata[,sapply(DEG_outdata,is.numeric)] <- signif(DEG_outdata[,sapply(DEG_outdata,is.numeric)],3)
  DT::datatable(DEG_outdata,extensions = 'Buttons',  options = list(
    dom = 'lBfrtip', buttons = c('csv', 'excel', 'print'), pageLength = 20), rownames= FALSE)
})

observeEvent(input$vennDiagram, {
	saved.num <- length(saved_plots$vennDiagram) + 1
	saved_plots$vennDiagram[[saved.num]] <- vennDiagram_out()
})

observeEvent(input$venn_DEG_data, {
  venndata <- DataVennReactive()
  vennlist <- venndata$vennlist
  allIDs=unique(unlist(vennlist))
  dataIn=DataReactive()
  data_results=dataIn$data_results
  all_names=names(data_results)
  tests=names(vennlist)
  selCol=NULL
  for (i in 1:length(tests)) {
    sel_i=which(str_detect(all_names, regex(str_c("^", tests[i]), ignore_case=T)))
    if (length(sel_i)>0) {selCol=c(selCol, sel_i)}
  }
  name_col=which(all_names %in% c("UniqueID", "Gene.Name") )
  sel_row=which(data_results$UniqueID %in% allIDs)
  #browser()#debug
  DEG_outdata=data_results[sel_row, c(name_col, selCol)]
  saved_table$DEG_outdata_Venn <- DEG_outdata
})


output$SvennDiagram <- renderPlot({
	print("drawing Venn diagram 2")
	venndata <- DataVennReactive()
	vennlist <- venndata$vennlist
	venn(vennlist, show.plot = TRUE, intersections = FALSE)
})

output$vennHTML <- renderText({

	venndata <- DataVennReactive()
	vennlist <- venndata$vennlist

	v.table <- venn(vennlist, show.plot = FALSE, intersections = TRUE)
	intersect <- attr(v.table,"intersections")
	
	venndata <- DataVennReactive()
	vennlist <- venndata$vennlist
	v.table <- venn(vennlist,show.plot = FALSE, intersections = TRUE)
	intersect <- attr(v.table,"intersections")
	htmlstr <- "  <br>"
	for (i in 1:length(intersect)) {
		if(input$vennlistname == "Gene"){
			intersectlist <- toString(sapply(strsplit(intersect[[i]],split= "\\_"),'[',1))
		} else if (input$vennlistname == "AC") {
			intersectlist <- toString(sapply(strsplit(intersect[[i]],split= "\\_"),'[',2))
		} else if (input$vennlistname == "UniqueID") {
			intersectlist <- toString(intersect[[i]])
		}
		htmlstr <-  paste(htmlstr,"<p><b><font color=red>", names(intersect[i]),"</font></b>:",intersectlist, sep="")
	}
	htmlstr
})



