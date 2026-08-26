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
GetRDataFile <- function(Pname) {
  if (Pname %in% projects) {
    RDataFile <- paste("data/",Pname,".RData", sep="")
  } else {
    RDataFile <- paste("unlisted/",Pname,".RData", sep="")
  }
}

observe({
	#data_sets <- list.files(path = "./data", pattern = "\\.RData$", full.names = FALSE) %>% gsub("\\.RData$","",.)
	data_sets  <- c("empty",projects)
	if (!is.null(pub_projects)) {
	  data_sets<-c(data_sets, pub_projects)
	}
	for (i in 1:5){
		dataset <- paste("dataset",i,sep="")
		updateSelectizeInput(session, dataset, choices=data_sets, selected="empty")
	}
})

observe({
	if(input$dataset1 != "empty" & input$dataset1 != "") {
		RDataFile <- GetRDataFile(input$dataset1) #paste("data/",input$dataset1,".RData", sep="")
		load(RDataFile)
		tests  <- as.character(MetaData$ComparePairs[MetaData$ComparePairs!=""])
		comp_tests=as.character(unique(results_long$test));  if (!all(tests %in% comp_tests) ) { tests <-  gsub("-", "vs", tests) } 
		if(length(tests)==0) {
			tests = unique(as.character(results_long$test))
		}

		updateSelectizeInput(session, "vennP_test1", choices=tests, selected=tests[1])
	}
})

observe({
	if(input$dataset2 != "empty" & input$dataset2 != "") {
		RDataFile <- GetRDataFile(input$dataset2) #paste("data/",input$dataset2,".RData", sep="")
		load(RDataFile)
		tests  <- as.character(MetaData$ComparePairs[MetaData$ComparePairs!=""])
		comp_tests=as.character(unique(results_long$test)); if (!all(tests %in% comp_tests) ) { tests <-  gsub("-", "vs", tests) } 
		if(length(tests)==0) {
			tests = unique(as.character(results_long$test))
		}

		updateSelectizeInput(session, "vennP_test2", choices=tests, selected=tests[1])
	}
})

observe({
	if(input$dataset3 != "empty" & input$dataset3 != "") {
		RDataFile <- GetRDataFile(input$dataset3) #paste("data/",input$dataset3,".RData", sep="")
		load(RDataFile)
		tests  <- as.character(MetaData$ComparePairs[MetaData$ComparePairs!=""])
		comp_tests=as.character(unique(results_long$test)); if (!all(tests %in% comp_tests) ) { tests <-  gsub("-", "vs", tests) } 
		if(length(tests)==0) {
			tests = unique(as.character(results_long$test))
		}

		updateSelectizeInput(session, "vennP_test3", choices=tests, selected=tests[1])
	}
})

observe({
	if(input$dataset4 != "empty" & input$dataset4 != "") {
		RDataFile <- GetRDataFile(input$dataset4) #paste("data/",input$dataset4,".RData", sep="")
		load(RDataFile)
		tests  <- as.character(MetaData$ComparePairs[MetaData$ComparePairs!=""])
		comp_tests=as.character(unique(results_long$test)); if (!all(tests %in% comp_tests) ) { tests <-  gsub("-", "vs", tests) } 
		if(length(tests)==0) {
			tests = unique(as.character(results_long$test))
		}

		updateSelectizeInput(session, "vennP_test4", choices=tests, selected=tests[1])
	}
})

observe({
	if(input$dataset5 != "empty" & input$dataset5 != "") {
		RDataFile <- GetRDataFile(input$dataset5) #paste("data/",input$dataset5,".RData", sep="")
		load(RDataFile)
		tests  <- as.character(MetaData$ComparePairs[MetaData$ComparePairs!=""])
		comp_tests=as.character(unique(results_long$test)); if (!all(tests %in% comp_tests) ) { tests <-  gsub("-", "vs", tests) } 
		if(length(tests)==0) {
			tests = unique(as.character(results_long$test))
		}

		updateSelectizeInput(session, "vennP_test5", choices=tests, selected=tests[1])
	}

})

DataVennPReactive <- reactive({
	vennP_fccut =log2(input$vennP_fccut)
	vennP_pvalcut = input$vennP_pvalcut
	data_sets  <- c("empty",projects)
	if (!is.null(pub_projects)) {
	  data_sets<-c(data_sets, pub_projects)
	}
	
	vennlist <- list()
	fill <- list()

	if (input$vennP_test1 != "Empty List" & input$vennP_test1 != "") {
		RDataFile <-  GetRDataFile(input$dataset1) #paste("data/",input$dataset1,".RData", sep="")
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

		listname1 <- paste(names(data_sets)[data_sets==input$dataset1],input$vennP_test1,sep="\n" )
		fill[[listname1]] <- input$col1
		vennlist[[listname1]]  <- list1
	}

	if (input$dataset2 != "empty" & input$dataset2 != ""& input$vennP_test2 != "Empty List" & input$vennP_test2 != "") {
		RDataFile <-  GetRDataFile(input$dataset2) #paste("data/",input$dataset2,".RData", sep="")
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
		listname2 <- paste(names(data_sets)[data_sets==input$dataset2],input$vennP_test2,sep="\n" )
		fill[[listname2]] <- input$col2
		vennlist[[listname2]]  <- list2
	}

	if (input$dataset3 != "empty" & input$dataset3 != ""& input$vennP_test3 != "Empty List" & input$vennP_test3 != "") {
		RDataFile <- GetRDataFile(input$dataset3) # paste("data/",input$dataset3,".RData", sep="")
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
		listname3 <- paste(names(data_sets)[data_sets==input$dataset3],input$vennP_test3,sep="\n" )
		fill[[listname3]] <- input$col3
		vennlist[[listname3]]  <- list3
	}

	if (input$dataset4 != "empty" & input$dataset4 != ""& input$vennP_test4 != "Empty List" & input$vennP_test4 != "") {
		RDataFile <-  GetRDataFile(input$dataset4) #paste("data/",input$dataset4,".RData", sep="")
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
		listname4 <- paste(names(data_sets)[data_sets==input$dataset4],input$vennP_test4,sep="\n" )
		fill[[listname4]] <- input$col4
		vennlist[[listname4]]  <- list4
	}

	if (input$dataset5 != "empty" & input$dataset5 != ""& input$vennP_test5 != "Empty List" & input$vennP_test5 != "") {
		RDataFile <-  GetRDataFile(input$dataset5) #paste("data/",input$dataset5,".RData", sep="")
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
		listname5 <- paste(names(data_sets)[data_sets==input$dataset5],input$vennP_test5,sep="\n" )
		fill[[listname5]] <- input$col5
		vennlist[[listname5]]  <- list5
	}
 if (input$upperSymbols) {
   for (i in 1:length(vennlist)) {
     vennlist[[i]]=toupper(vennlist[[i]])
   }
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
	# Labels here are already "<dataset>\n<comparison>" (see DataVennPReactive
	# above) -- wrap each of those two lines separately so a long dataset name
	# or comparison name doesn't overlap adjacent labels, without merging the
	# two into one run-on line. Wrap width is user-adjustable (input$vennPcatwrap).
	wrap_label <- function(x, width) {
		parts <- strsplit(x, "\n", fixed = TRUE)[[1]]
		wrap_hard <- function(p) {
			chunks <- regmatches(p, gregexpr(paste0(".{1,", width, "}"), p, perl = TRUE))[[1]]
			paste(chunks, collapse = "\n")
		}
		paste(vapply(parts, wrap_hard, character(1)), collapse = "\n")
	}
	cat_names <- vapply(names(vennlist), wrap_label, character(1), width = input$vennPcatwrap)
	venn.plot <- venn.diagram(x = vennlist,
		category.names = cat_names,
		fill=fillcor,
		lty=input$vennPlty, lwd=input$vennPlwd, alpha=input$vennPalpha,
		cex=input$vennPcex, cat.cex=input$vennPcatcex, margin=input$vennPmargin,
		fontface = input$vennPfontface, cat.fontface=input$vennPcatfontface,
		main = input$vennPtitle, main.cex = input$vennPmaincex, main.pos = c(0.5, 1.1), main.fontface = "bold",
	filename = NULL);
	grid.newpage()
	# Reserve ~12% of height at the bottom so wrapped category labels
	# that extend below the diagram aren't cropped by the image boundary.
	pushViewport(viewport(x = 0.5, y = 0.56, width = 1, height = 0.88, clip = "off"))
	grid.draw(venn.plot);
	popViewport()
})

output$SvennPDiagram <- renderPlot({
	print("drawing Venn diagram 2")
	venndata <- DataVennPReactive()
	vennlist <- venndata$vennlist
	validate(need(length(vennlist)>=1, message = "Select projects."))
	venn(vennlist, show.plot = TRUE, intersections = FALSE)
})

#' Intersection Output as a data.table -- see VennIntersectReactive() in
#' venn.R for the current-project equivalent. vennlist's values here are
#' already Gene.Name strings (see DataVennPReactive above), so there's no
#' Gene.Name/AC/UniqueID choice to make like the single-project tab has.
VennPIntersectReactive <- reactive({
	venndata <- DataVennPReactive()
	vennlist <- venndata$vennlist
	validate(need(length(vennlist)>=2, message = "Select at least 2 projects/comparisons to see intersections."))
	v.table <- venn(vennlist,show.plot = FALSE, intersections = TRUE)
	intersect <- attr(v.table,"intersections")
	data.frame(
		Intersection = names(intersect),
		Count = lengths(intersect),
		Genes = vapply(intersect, toString, character(1)),
		stringsAsFactors = FALSE, check.names = FALSE
	)
})

output$vennP_intersect_table <- DT::renderDataTable({
	DT::datatable(VennPIntersectReactive(), rownames = FALSE, selection = "multiple",
		options = list(pageLength = 20, dom = "lfrtip"))
})

output$vennP_copy_btn <- renderUI({
	df <- VennPIntersectReactive()
	sel <- input$vennP_intersect_table_rows_selected
	req(length(sel) > 0)
	genes <- unique(unlist(strsplit(df$Genes[sel], ", ", fixed = TRUE)))
	rclipboard::rclipButton(
		"vennP_copy_genes_btn",
		label = paste0("Copy Selected Gene List (", length(genes), " genes)"),
		clipText = paste(genes, collapse = ","),
		icon = icon("copy", lib = "glyphicon"),
		class = "btn-primary btn-sm"
	)
})



