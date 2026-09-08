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
	req(tmptests)
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
	# Wrap long comparison names onto multiple lines so Venn diagram category
	# labels don't overlap -- wrap width is user-adjustable (input$catwrap).
	# Only affects the displayed labels; vennlist's own names (used for fill
	# color matching and downstream table/intersection lookups) are untouched.
	wrap_hard <- function(x, width) {
		chunks <- regmatches(x, gregexpr(paste0(".{1,", width, "}"), x, perl = TRUE))[[1]]
		paste(chunks, collapse = "\n")
	}
	cat_names <- vapply(names(vennlist), wrap_hard, character(1), width = input$catwrap)
	venn.plot <- venn.diagram(x = vennlist,
		category.names = cat_names,
		fill=fillcor, margin=input$margin,
		lty=input$lty, lwd=input$lwd, alpha=input$alpha,
		cex=input$cex, cat.cex=input$catcex,
		fontface = input$fontface, cat.fontface=input$catfontface,
		main = input$title, main.cex = input$maincex, main.pos = c(0.5, 1.1), main.fontface = "bold",
	filename = NULL)

	return(venn.plot)
})

output$vennDiagram <- renderPlot({
	grid.newpage()
	# Reserve ~12% of height at the bottom so wrapped category labels
	# that extend below the diagram aren't cropped by the image boundary.
	pushViewport(viewport(x = 0.5, y = 0.56, width = 1, height = 0.88, clip = "off"))
	grid.draw(vennDiagram_out())
	popViewport()
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

#' Intersection Output as a data.table: one row per non-empty Venn subset,
#' with a comma-joined gene list (respecting the Gene.Name/AC/UniqueID
#' choice, same as the old vennHTML text did) plus its size, so the DT UI
#' below can offer row selection and a "copy gene list" button.
VennIntersectReactive <- reactive({
	DataIn = DataReactive()
	ProteinGeneName = DataIn$ProteinGeneName

	venndata <- DataVennReactive()
	vennlist <- venndata$vennlist
	validate(need(length(vennlist) >= 2, "Select at least 2 comparisons to see intersections."))

	v.table <- venn(vennlist, show.plot = FALSE, intersections = TRUE)
	intersect <- attr(v.table,"intersections")

	gene_lists <- lapply(intersect, function(ids) {
		if (input$vennlistname == "Gene.Name") {
			ProteinGeneName %>% dplyr::filter(UniqueID %in% ids) %>% dplyr::pull(Gene.Name)
		} else if (input$vennlistname == "AC") {
			sapply(strsplit(ids, split = "\\_"), '[', 2)
		} else {
			ids
		}
	})

	data.frame(
		Intersection = names(intersect),
		Count = lengths(gene_lists),
		Genes = vapply(gene_lists, toString, character(1)),
		stringsAsFactors = FALSE, check.names = FALSE
	)
})

output$venn_intersect_table <- DT::renderDataTable({
	df <- VennIntersectReactive()
	# Size the Intersection column to the longest individual comparison name
	# actually selected (venn_test1..5), so short names get a compact column
	# and longer ones aren't hard-wrapped mid-word unnecessarily.
	col_width_px <- intersection_column_width_px(names(DataVennReactive()$vennlist))
	# Break the Intersection label at each ":" so a deep overlap name (e.g.
	# "A:B:C:D") wraps onto its own lines in a narrow column instead of
	# forcing the whole table wide. escape=-1 below leaves these <br> tags
	# (and only these) un-escaped.
	df$Intersection <- gsub(":", "<br>", df$Intersection, fixed = TRUE)

	DT::datatable(df, rownames = FALSE, selection = "multiple", escape = -1,
		options = list(
			pageLength = 20, dom = "lfrtip",
			columnDefs = list(
				list(targets = 0, width = paste0(col_width_px, "px"), className = "wrap-cell"),
				# Genes column: show only the first 50 genes with a "Show all"
				# link appended, generated client-side at draw time -- the
				# underlying data (read via df$Genes[sel] in venn_copy_btn
				# below) always stays the full, untruncated list; only what's
				# drawn on screen for type=="display" is shortened.
				list(targets = 2, render = venn_genes_show_all_js)
			)
		))
})

# clipText is baked into the button's HTML at render time (clipboard.js has
# no server round-trip), so the button itself has to be regenerated via
# renderUI every time the row selection changes.
output$venn_copy_btn <- renderUI({
	df <- VennIntersectReactive()
	sel <- input$venn_intersect_table_rows_selected
	req(length(sel) > 0)
	genes <- unique(unlist(strsplit(df$Genes[sel], ", ", fixed = TRUE)))
	rclipboard::rclipButton(
		"venn_copy_genes_btn",
		# "unique" called out explicitly: this count can differ from the
		# Count column / "Show all (N genes)" link, which are both based on
		# UniqueID membership -- two different UniqueIDs (e.g. distinct
		# probes/isoforms) can share the same Gene.Name, so the deduplicated
		# copy list can be shorter than the row's raw intersection size.
		label = paste0("Copy Selected Gene List (", length(genes), " unique genes)"),
		clipText = paste(genes, collapse = ","),
		icon = icon("copy", lib = "glyphicon"),
		class = "btn-primary btn-sm"
	)
})



