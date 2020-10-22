###########################################################################################################
## Proteomics Visualization R Shiny App
##
##This software belongs to Biogen Inc. All right reserved.
##
##@file: vennprojects.R
##@Developer : Benbo Gao (benbo.gao@Biogen.com)
##@Date : 5/16/2018
##@version 1.0
###########################################################################################################


observe({
	#data_sets <- list.files(path = "./data", pattern = "\\.RData$", full.names = FALSE) %>% gsub("\\.RData$","",.)
	data_sets  <- c("empty",projects)
	for (i in 1:length(data_sets)){
		dataset <- paste("dataset",i,sep="")
		updateSelectizeInput(session, dataset, choices=data_sets, selected="empty")
	}
})

observe({
	if(input$dataset1 != "empty" & input$dataset1 != "") {
		RDataFile <- paste("data/",input$dataset1,".RData", sep="")
		load(RDataFile)
		tests  <- as.character(MetaData$ComparePairs[MetaData$ComparePairs!=""])
		tests <-  gsub("-", "vs", tests)
		if(length(tests)==0) {
			tests = unique(as.character(results_long$test))
		}

		updateSelectizeInput(session, "vennP_test1", choices=tests, selected=tests[1])
	}
})

observe({
	if(input$dataset2 != "empty" & input$dataset2 != "") {
		RDataFile <- paste("data/",input$dataset2,".RData", sep="")
		load(RDataFile)
		tests  <- as.character(MetaData$ComparePairs[MetaData$ComparePairs!=""])
		tests <-  gsub("-", "vs", tests)
		if(length(tests)==0) {
			tests = unique(as.character(results_long$test))
		}

		updateSelectizeInput(session, "vennP_test2", choices=tests, selected=tests[1])
	}
})

observe({
	if(input$dataset3 != "empty" & input$dataset3 != "") {
		RDataFile <- paste("data/",input$dataset3,".RData", sep="")
		load(RDataFile)
		tests  <- as.character(MetaData$ComparePairs[MetaData$ComparePairs!=""])
		tests <-  gsub("-", "vs", tests)
		if(length(tests)==0) {
			tests = unique(as.character(results_long$test))
		}

		updateSelectizeInput(session, "vennP_test3", choices=tests, selected=tests[1])
	}
})

observe({
	if(input$dataset4 != "empty" & input$dataset4 != "") {
		RDataFile <- paste("data/",input$dataset4,".RData", sep="")
		load(RDataFile)
		tests  <- as.character(MetaData$ComparePairs[MetaData$ComparePairs!=""])
		tests <-  gsub("-", "vs", tests)
		if(length(tests)==0) {
			tests = unique(as.character(results_long$test))
		}

		updateSelectizeInput(session, "vennP_test4", choices=tests, selected=tests[1])
	}
})

observe({
	if(input$dataset5 != "empty" & input$dataset5 != "") {
		RDataFile <- paste("data/",input$dataset5,".RData", sep="")
		load(RDataFile)
		tests  <- as.character(MetaData$ComparePairs[MetaData$ComparePairs!=""])
		tests <-  gsub("-", "vs", tests)
		if(length(tests)==0) {
			tests = unique(as.character(results_long$test))
		}

		updateSelectizeInput(session, "vennP_test5", choices=tests, selected=tests[1])
	}

})

DataVennPReactive <- reactive({
	vennP_fccut =log2(input$vennP_fccut)
	vennP_pvalcut = input$vennP_pvalcut
	vennlist <- list()
	fill <- list()

	if (input$vennP_test1 != "Empty List" & input$vennP_test1 != "") {
		RDataFile <- paste("data/",input$dataset1,".RData", sep="")
		load(RDataFile)
		if (input$vennP_psel == "Padj") {
			results_long <-  results_long %>% dplyr::filter(abs(logFC) > vennP_fccut & Adj.P.Value < vennP_pvalcut)
		} else{
			results_long <-  results_long %>% dplyr::filter(abs(logFC) > vennP_fccut & P.Value < vennP_pvalcut)
		}



		list1 = results_long %>%
		filter(test == input$vennP_test1) %>%
		dplyr::left_join(.,ProteinGeneName,by="UniqueID") %>%
		dplyr::select(Gene.Name) %>%	collect %>%	.[["Gene.Name"]] %>% as.character()%>% unique()

		listname1 <- paste(input$dataset1,input$vennP_test1,sep="\n" )
		fill[[listname1]] <- input$col1
		vennlist[[listname1]]  <- list1
	}

	if (input$dataset2 != "empty" & input$dataset2 != ""& input$vennP_test2 != "Empty List" & input$vennP_test2 != "") {
		RDataFile <- paste("data/",input$dataset2,".RData", sep="")
		load(RDataFile)

		if (input$vennP_psel == "Padj") {
			results_long <-  results_long %>% dplyr::filter(abs(logFC) > vennP_fccut & Adj.P.Value < vennP_pvalcut)
		} else{
			results_long <-  results_long %>% dplyr::filter(abs(logFC) > vennP_fccut & P.Value < vennP_pvalcut)
		}

		list2 = results_long %>%
		filter(test == input$vennP_test2) %>%
		dplyr::left_join(.,ProteinGeneName,by="UniqueID") %>%
		dplyr::select(Gene.Name) %>%	collect %>%	.[["Gene.Name"]] %>% as.character()%>% unique()
		listname2 <- paste(input$dataset2,input$vennP_test2,sep="\n" )
		fill[[listname2]] <- input$col2
		vennlist[[listname2]]  <- list2
	}

	if (input$dataset3 != "empty" & input$dataset3 != ""& input$vennP_test3 != "Empty List" & input$vennP_test3 != "") {
		RDataFile <- paste("data/",input$dataset3,".RData", sep="")
		load(RDataFile)

		if (input$vennP_psel == "Padj") {
			results_long <-  results_long %>% dplyr::filter(abs(logFC) > vennP_fccut & Adj.P.Value < vennP_pvalcut)
		} else{
			results_long <-  results_long %>% dplyr::filter(abs(logFC) > vennP_fccut & P.Value < vennP_pvalcut)
		}

		list3 = results_long %>%
		filter(test == input$vennP_test3) %>%
		dplyr::left_join(.,ProteinGeneName,by="UniqueID") %>%
		dplyr::select(Gene.Name) %>%	collect %>%	.[["Gene.Name"]] %>% as.character()%>% unique()
		listname3 <- paste(input$dataset3,input$vennP_test3,sep="\n" )
		fill[[listname3]] <- input$col3
		vennlist[[listname3]]  <- list3
	}

	if (input$dataset4 != "empty" & input$dataset4 != ""& input$vennP_test4 != "Empty List" & input$vennP_test4 != "") {
		RDataFile <- paste("data/",input$dataset4,".RData", sep="")
		load(RDataFile)

		if (input$vennP_psel == "Padj") {
			results_long <-  results_long %>% dplyr::filter(abs(logFC) > vennP_fccut & Adj.P.Value < vennP_pvalcut)
		} else{
			results_long <-  results_long %>% dplyr::filter(abs(logFC) > vennP_fccut & P.Value < vennP_pvalcut)
		}

		list4 = results_long %>%
		filter(test == input$vennP_test4) %>%
		dplyr::left_join(.,ProteinGeneName,by="UniqueID") %>%
		dplyr::select(Gene.Name) %>%	collect %>%	.[["Gene.Name"]] %>% as.character()%>% unique()
		listname4 <- paste(input$dataset4,input$vennP_test4,sep="\n" )
		fill[[listname4]] <- input$col4
		vennlist[[listname4]]  <- list4
	}

	if (input$dataset5 != "empty" & input$dataset5 != ""& input$vennP_test5 != "Empty List" & input$vennP_test5 != "") {
		RDataFile <- paste("data/",input$dataset5,".RData", sep="")
		load(RDataFile)

		if (input$vennP_psel == "Padj") {
			results_long <-  results_long %>% dplyr::filter(abs(logFC) > vennP_fccut & Adj.P.Value < vennP_pvalcut)
		} else{
			results_long <-  results_long %>% dplyr::filter(abs(logFC) > vennP_fccut & P.Value < vennP_pvalcut)
		}

		list5 = results_long %>%
		filter(test == input$vennP_test5) %>%
		dplyr::left_join(.,ProteinGeneName,by="UniqueID") %>%
		dplyr::select(Gene.Name) %>%	collect %>%	.[["Gene.Name"]] %>% as.character()%>% unique()
		listname5 <- paste(input$dataset5,input$vennP_test5,sep="\n" )
		fill[[listname5]] <- input$col5
		vennlist[[listname5]]  <- list5
	}

	return(venndata = list("vennlist"=vennlist, "fillcor"=fill) )
})

output$vennPDiagram <- renderPlot({
	print("drawing Venn diagram")
	venndata <- DataVennPReactive()
	vennlist <- venndata$vennlist
	validate(need(length(vennlist)>=1, message = "Select projects."))

	fillcor <- unlist(venndata$fillcor)
	SetNum = length(vennlist)
	futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
	venn.plot <- venn.diagram(x = vennlist,
		fill=fillcor,
		lty=input$vennPlty, lwd=input$vennPlwd, alpha=input$vennPalpha,
		cex=input$vennPcex, cat.cex=input$vennPcatcex, margin=input$vennPmargin,
		fontface = input$vennPfontface, cat.fontface=input$vennPcatfontface,
		main = input$vennPtitle, main.cex = input$vennPmaincex, main.pos = c(0.5, 1.1), main.fontface = "bold",
	filename = NULL);
	grid.draw(venn.plot);
})

output$SvennPDiagram <- renderPlot({
	print("drawing Venn diagram 2")
	venndata <- DataVennPReactive()
	vennlist <- venndata$vennlist
	validate(need(length(vennlist)>=1, message = "Select projects."))
	venn(vennlist, show.plot = TRUE, intersections = FALSE)
})

output$vennPHTML <- renderText({
	venndata <- DataVennPReactive()
	vennlist <- venndata$vennlist
	validate(need(length(vennlist)>=1, message = "Select projects."))
	v.table <- venn(vennlist,show.plot = FALSE, intersections = TRUE)
	intersect <- attr(v.table,"intersections")
	htmlstr <- "  <br>"
	for (i in 1:length(intersect)) {
		htmlstr <-  paste(htmlstr,"<p><b><font color=red>", names(intersect[i]),"</font></b>:",toString(intersect[[i]]), sep="")
	}
	htmlstr
})



