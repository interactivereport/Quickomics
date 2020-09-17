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
  #cat("show all samples", length(samples), length(groups), "\n") #debug
  updateSelectizeInput(session,'heatmap_samples', choices=samples, selected=samples)
  updateSelectizeInput(session,'heatmap_test',choices=tests, selected=tests[1])
  ProteinGeneName = DataIn$ProteinGeneName
  updateRadioButtons(session,'heatmap_label', inline = TRUE, choices=colnames(ProteinGeneName)[-1])
})

observe({
  DataIn = DataReactive()
  tmpgroups = input$heatmap_groups
  tmpdat = DataIn$MetaData %>% filter(group %in% tmpgroups)
  tmpsamples = as.character(tmpdat$sampleid)
  #cat("update samples", length(tmpsamples), length(tmpgroups), "\n") #debug
  updateSelectizeInput(session,'heatmap_samples', choices=tmpsamples, selected=tmpsamples)
})


filteredGene=reactive({
  heatmap_test = input$heatmap_test
  heatmap_fccut = as.numeric(input$heatmap_fccut)
  heatmap_pvalcut = as.numeric(input$heatmap_pvalcut)
  DataIn = DataReactive()
  results_long = DataIn$results_long
  
  if (input$heatmap_psel == "Padj") {
    filteredGene = results_long %>% filter(test %in% heatmap_test & abs(logFC) > heatmap_fccut & Adj.P.Value < heatmap_pvalcut) %>%
      dplyr::select(UniqueID) %>% 	collect %>%	.[["UniqueID"]] %>%	as.character()
  } else {
    filteredGene = results_long %>% filter(test %in% heatmap_test & abs(logFC) > heatmap_fccut & P.Value < heatmap_pvalcut) %>%
      dplyr::select(UniqueID) %>% 	collect %>%	.[["UniqueID"]] %>%	as.character()
  }
  #cat("Selected Genes:",length(filteredGene), "\n") #debug
  return(filteredGene)
})

output$heatmapfilteredgene <- renderText({ paste("Selected Genes:",length(filteredGene()),sep="")})


DataHeatMapReactive <- reactive({
  validate(need(input$heatmap_groups, FALSE))
  validate(need(input$heatmap_samples, FALSE))
  DataIn = DataReactive()
  results_long = DataIn$results_long
  ProteinGeneName = DataIn$ProteinGeneName
  MetaData = DataIn$MetaData
  #cat("work on Data for Heatmap", date(), "\n") #debug
  tmpgroups = input$heatmap_groups
  #group_order(input$heatmap_groups)
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
    if(length(filteredGene())>0) {
      tmpdat  <-  tmpdat[filteredGene(),]
    }
  }
  
  if (input$heatmap_subset == "all") {
    if (nrow(tmpdat)>input$maxgenes) {
      if (input$heatmap_submethod=="Random") {
        tmpdat=tmpdat[sample(1:nrow(tmpdat), input$maxgenes),] #this will keep rownames
        #tmpdat <- tmpdat %>% sample_n(input$maxgenes) #this will remove rownames
      } else {
        dataSD=apply(tmpdat, 1, function(x) sd(x,na.rm=T))
        dataM=rowMeans(tmpdat)
        diff=dataSD/(dataM+median(dataM)) #SD/mean, added median value to penalized lower expressed genes
        tmpdat=tmpdat[order(diff, decreasing=TRUE)[1:input$maxgenes], ]	    
      }
    }
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
  
  if (nrow(tmpdat)>5000 ) {tmpdat=tmpdat[sample(1:nrow(tmpdat), 5000),]; cat("Reduce data pionts to 5K\n")} #Use at most 5000 genes so the App won't crash
  
  df <- data.matrix(tmpdat)
  #use selected gene label
  sel=match(rownames(df), ProteinGeneName$UniqueID)
  selCol=match(input$heatmap_label, names(ProteinGeneName))
  
  if (sum(is.na(sel))==0 & sum(is.na(selCol)==0)) {rownames(df)=unlist(ProteinGeneName[sel, selCol])
  } else {cat("gene lables not updated",sum(is.na(sel)), sum(is.na(selCol)), "\n")}
  
  return(list("df"=df, "annotation"=annotation))
})




pheatmap2_out <- eventReactive(input$plot_heatmap, {
  withProgress(message = 'Making static heatmap 1:', value = 0, {
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
    # cat("pheatmap ", dim(data.in), date(), "\n") #debug
    if (cexRow  == 0 | nrow(data.in) > 100) {
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
    
    #clean up SD=0 rows and columns
    if (input$scale=="row") {
      row_SD=apply(data.in, 1, function(x) sd(x,na.rm=T))
      data.in=data.in[row_SD!=0, ]
    }
    if (input$scale=="column") {
      col_SD=apply(data.in, 2, function(x) sd(x,na.rm=T))
      data.in=data.in[, col_SD!=0]
    }	
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
    output$pheatmap2 <- renderPlot({
      grid.draw(p$gtable)})
    return(p)
  })
})


observeEvent(input$plot_heatmap, {  
  output$pheatmap2 <- renderPlot({
    grid.draw(pheatmap2_out()$gtable)})
})

observeEvent(input$pheatmap2, {
  saved_plots$pheatmap2 <- pheatmap2_out()$gtable
}
)

staticheatmap_out <- reactive({
  withProgress(message = 'Making static heatmap 2:', value = 0, {
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
    
    #  cat(date(), dim(data.in), "layout 2\n") #debug
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