###########################################################################################################
## Proteomics Visualization R Shiny App
##
##This software belongs to Biogen Inc. All right reserved.
##
##@file: heatmap.R
##@Developer : Benbo Gao (benbo.gao@Biogen.com)
##@Date : 5/16/2018
##@version 1.0
###########################################################################################################


observe({
	DataIn = DataReactive()
	groups = group_order()
	samples = DataIn$MetaData$sampleid
	tests = DataIn$tests
	allgroups = DataIn$groups
	updateSelectizeInput(session,'heatmap_groups', choices=allgroups, selected=groups)
	updateSelectizeInput(session,'heatmap_samples', choices=samples, selected=samples)
	updateSelectizeInput(session,'heatmap_test',choices=tests, selected=tests[1])
})

observe({
	DataIn = DataReactive()
	tmpgroups = input$heatmap_groups
	tmpdat = DataIn$MetaData %>% filter(group %in% tmpgroups)
	tmpsamples = as.character(tmpdat$sampleid)
	updateSelectizeInput(session,'heatmap_samples', choices=tmpsamples, selected=tmpsamples)
})


DataHeatMapReactive <- reactive({
	DataIn = DataReactive()
	results_long = DataIn$results_long
	ProteinGeneName = DataIn$ProteinGeneName
	MetaData = DataIn$MetaData

	tmpgroups = input$heatmap_groups
	group_order(input$heatmap_groups)
	tmpsamples = input$heatmap_samples
	tmpkeep = which((MetaData$group %in% tmpgroups)&(MetaData$sampleid %in% tmpsamples))

	tmp_group = MetaData$group[tmpkeep]
	tmp_sampleid = MetaData$sampleid[tmpkeep]
	annotation = data.frame("group" = tmp_group)
	rownames(annotation) <- tmp_sampleid

	if(length(tmpkeep)>0) {
		y <- input$heatmap_groups
		x= MetaData$group[tmpkeep]
		z = MetaData$sampleid[tmpkeep]
		new_order <- as.character(z[order(match(x, y))])
		tmpdat  <- DataIn$data_wide %>% dplyr::select(new_order)
		
		
		tmpdat[is.na(tmpdat)] <- 0
		rownames(tmpdat) <-  rownames(DataIn$data_wide)
	}

	if (input$heatmap_subset == "subset") {
		heatmap_test = input$heatmap_test
		heatmap_fccut = as.numeric(input$heatmap_fccut)
		heatmap_pvalcut = as.numeric(input$heatmap_pvalcut)

		if (input$heatmap_psel == "Padj") {
			filteredGene = results_long %>% filter(test %in% heatmap_test & abs(logFC) > heatmap_fccut & Adj.P.Value < heatmap_pvalcut) %>%
			dplyr::select(UniqueID) %>% 	collect %>%	.[["UniqueID"]] %>%	as.character()
		} else {
			filteredGene = results_long %>% filter(test %in% heatmap_test & abs(logFC) > heatmap_fccut & P.Value < heatmap_pvalcut) %>%
			dplyr::select(UniqueID) %>% 	collect %>%	.[["UniqueID"]] %>%	as.character()
		}

		output$heatmapfilteredgene <- renderText({ paste("Selected Genes:",length(filteredGene),sep="")})

		if(length(filteredGene)>0) {
			tmpdat  <-  tmpdat[filteredGene,]
		}
	}

	if (input$heatmap_subset == "all") {
		#tmpdat=tmpdat[sample(1:nrow(tmpdat), input$maxgenes),]
		tmpdat <- tmpdat %>% sample_n(input$maxgenes)
		
	}

	if (input$heatmap_subset == "upload genes") {
		heatmap_list <- input$heatmap_list
		if(grepl("\n",heatmap_list)) {
			heatmap_list <-  stringr::str_split(heatmap_list, "\n")[[1]]
		} else if(grepl(",",heatmap_list)) {
			heatmap_list <-  stringr::str_split(heatmap_list, ",")[[1]]
		}

		heatmap_list <- gsub(" ", "", heatmap_list, fixed = TRUE)
		heatmap_list <- unique(heatmap_list[heatmap_list != ""])

		validate(need(length(heatmap_list)>2, message = "input gene list"))

		uploadlist <- dplyr::filter(ProteinGeneName, (UniqueID %in% heatmap_list) | (Protein.ID %in% heatmap_list) | (Gene.Name %in% heatmap_list))  %>%
		dplyr::select(UniqueID) %>% 	collect %>%	.[["UniqueID"]] %>%	as.character()
		tmpdat  <-  tmpdat[uploadlist,]
		
		
		
	}


	df <- data.matrix(tmpdat)
	return(list("df"=df, "annotation"=annotation))
})

pheatmap2_out <- reactive({
	withProgress(message = 'Making static heatmap:', value = 0, {
		DataHeatMap <- DataHeatMapReactive()
		data.in <- DataHeatMap$df
		annotation <- DataHeatMap$annotation

		cluster_rows = 	cluster_cols = FALSE
		if (input$dendrogram == "both" | input$dendrogram == "row")
		cluster_rows = TRUE
		if (input$dendrogram == "both" | input$dendrogram == "column")
		cluster_cols = TRUE

		cexRow = as.numeric(as.character(input$hyfontsizep))
		cexCol = as.numeric(as.character(input$hxfontsizep))

		labCol = TRUE
		labRow = TRUE

		if (cexRow  == 0 | nrow(data.in) > 50) {
			labRow = FALSE
			cexRow = 5
		}
		if (cexCol == 0) {
			labCol = FALSE
			cexCol  = 5
		}

		cutree_rows = input$cutreerows
		cutree_cols = input$cutreecols

		if (cutree_rows == 0)
		cutree_rows = NA
		if (cutree_cols == 0)
		cutree_cols = NA

		p <- pheatmap(data.in,
			color = colorpanel (32, low = input$lowColor,mid = input$midColor, high = input$highColor),
			kmeans_k = NA, breaks = NA, border_color = "grey60",
			cellwidth = NA, cellheight = NA,
			scale = input$scale,
			cluster_rows = cluster_rows,
			cluster_cols = cluster_cols,
			clustering_distance_rows = input$distanceMethod,
			clustering_distance_cols = input$distanceMethod,
			clustering_method = input$agglomerationMethod,
			cutree_rows = cutree_rows,
			cutree_cols = cutree_cols,
			legend = TRUE, legend_breaks = NA,
			legend_labels = NA, annotation_row = NA, annotation_col = annotation,
			annotation_names_row = TRUE, annotation_names_col = TRUE,
			drop_levels = TRUE,
			show_rownames = labRow,
			show_colnames = labCol,
			main = NA,
			fontsize = 10,
			fontsize_row = cexRow,
			fontsize_col = cexCol,
			display_numbers = F, number_format = "%.2f", number_color = "grey30",
			labels_row = NULL, labels_col = NULL, filename = NA,
		silent = FALSE)
		return(p)
	})
})


output$pheatmap2 <- renderPlot({
	grid.draw(pheatmap2_out()$gtable)
})

observeEvent(input$pheatmap2, {
	saved_plots$pheatmap2 <- pheatmap2_out()$gtable
}
)

staticheatmap_out <- reactive({
	withProgress(message = 'Making static heatmap:', value = 0, {
		DataHeatMap <- DataHeatMapReactive()
		data.in <- DataHeatMap$df
		annotation <- DataHeatMap$annotation

		cutree_rows = input$cutreerows
		cutree_cols = input$cutreecols
		if (cutree_rows == 0)
		cutree_rows = NULL
		if (cutree_cols == 0)
		cutree_cols = NULL

		if (input$dendrogram == "both" | input$dendrogram == "row")
		dend_r <- data.in %>% dist(method = input$distanceMethod) %>% hclust(method = input$agglomerationMethod) %>% as.dendrogram 
		#%>% ladderize %>%  color_branches(k=cutree_rows)
		if (input$dendrogram == "both" | input$dendrogram == "column")
		dend_c <- t(data.in) %>% dist(method = input$distanceMethod) %>% hclust(method = input$agglomerationMethod) %>% as.dendrogram 
		#%>% ladderize %>% color_branches(k=cutree_cols)


		cexRow = as.numeric(as.character(input$hyfontsize))
		cexCol = as.numeric(as.character(input$hxfontsize))

		labCol = labRow = NULL

		if (cexRow  == 0 | nrow(data.in) > 50) {
			labRow = FALSE
			cexRow = 0.2
		}

		if (cexCol == 0) {
			labCol = FALSE
			cexCol  = 0.2
		}

		p<-	heatmap.2(
			data.in,
			trace = "none",
			scale = input$scale,
			dendrogram = input$dendrogram,
			key = input$key,
			labRow = labRow,
			labCol = labCol,
			cexRow = cexRow,
			cexCol = cexCol,
			Rowv = if (input$dendrogram == "both" | input$dendrogram == "row") dend_r else FALSE,
			Colv = if (input$dendrogram == "both" | input$dendrogram == "column") dend_c else FALSE,
			col = colorpanel (32, low = input$lowColor,mid = input$midColor, high = input$highColor),
			srtCol = as.numeric(as.character(input$srtCol)),
			margins = c(input$bottom,input$right)
		)
		obj = recordPlot()
		return(obj)
	})
})


output$staticheatmap <- renderPlot({
	replayPlot(staticheatmap_out())
})

observeEvent(input$staticheatmap, {
	saved_plots$staticheatmap <- staticheatmap_out()
}
)

interactiveHeatmap <- eventReactive(input$action_heatmaps, {
	DataHeatMap <- DataHeatMapReactive()
	data.in <- DataHeatMap$df
	annotation <- DataHeatMap$annotation
	cutree_rows = input$cutreerows
	cutree_cols = input$cutreecols
	if (cutree_rows == 0)
	cutree_rows = NULL
	if (cutree_cols == 0)
	cutree_cols = NULL


	if (input$dendrogram == "both" | input$dendrogram == "row")
	dend_r <- data.in %>% dist(method = input$distanceMethod) %>% hclust(method = input$agglomerationMethod) %>% as.dendrogram %>% ladderize %>%   color_branches(k=cutree_rows)
	if (input$dendrogram == "both" | input$dendrogram == "column")
	dend_c <- t(data.in) %>% dist(method = input$distanceMethod) %>% hclust(method = input$agglomerationMethod) %>% as.dendrogram %>% ladderize %>% color_branches(k=cutree_cols)

	cexRow = as.numeric(as.character(input$hyfontsizei))
	cexCol = as.numeric(as.character(input$hxfontsizei))

	labCol = colnames(data.in)
	labRow = rownames(data.in)


	if (cexRow  == 0 | nrow(data.in) > 50) {
		labRow = NA
		cexRow = 0.2
	}

	if (cexCol == 0) {
		labCol = NA
		cexCol  = 0.2
	}

	hide_colorbar=FALSE
	if (input$key == "FALSE")
	hide_colorbar=TRUE

	heatmaply(data.in,
		dendrogram = input$dendrogram,
		colors=colorpanel (32, low = input$lowColor,mid = input$midColor, high = input$highColor),
		Rowv = if (input$dendrogram == "both" | input$dendrogram == "row") dend_r else FALSE,
		Colv = if (input$dendrogram == "both" | input$dendrogram == "column") dend_c else FALSE,
		labRow = labRow,
		labCol = labCol,
		cexRow = cexRow,
		cexCol = cexCol,
		srtCol = as.numeric(as.character(input$srtCol)),
		hide_colorbar = hide_colorbar
	) %>% layout(margin = list(l = input$l, b = input$b))
	}
)

output$interactiveheatmap <- renderPlotly({
	withProgress(message = 'Making interactive heatmap:', value = 0, {
		interactiveHeatmap()
	})
})

output$text <- renderText({ "Click Generate Interactive Heatmap to view. (Disabled. This function is slow)"})