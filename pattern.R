###########################################################################################################
## Proteomics Visualization R Shiny App
##
##This software belongs to Biogen Inc. All right reserved.
##
##@file: pattern.R
##@Developer : Benbo Gao (benbo.gao@Biogen.com)
##@Date : 5/16/2018
##@version 1.0
###########################################################################################################

observeEvent(DataQCReactive(), {
  req(DataQCReactive())
	DataIn = DataQCReactive()
	tests = c("ALL",test_order())
	MetaData=DataIn$MetaData
	attributes=sort(setdiff(colnames(MetaData), c("sampleid", "Order", "ComparePairs")))
	updateSelectInput(session, "pattern_attr", choices=attributes, selected="group")  
	updateSelectizeInput(session,'pattern_test',choices=tests, selected=tests[1])
})

observeEvent(list(input$pattern_attr, DataQCReactive()), {
  req(input$pattern_attr != "")
  req(DataQCReactive())
  DataIn = DataQCReactive()
  MetaData=DataIn$MetaData
  allgroups = unique(MetaData[,input$pattern_attr])
  if (input$pattern_attr %in% names(DataIn$tmp_group)) {
    selected_groups <- DataIn$tmp_group[[input$pattern_attr]]
  } else {
    selected_groups <- allgroups
  }
  removed_groups <- allgroups[!(allgroups %in% selected_groups)]

  # If a bookmarked selection exists for this exact attribute, prefer it
  # (validated against the current data) over the "all active groups"
  # default above -- see restored_pattern_group_snapshot's definition in
  # app.R. Only applies once pattern_attr has itself been restored back to
  # the matching attribute, since group_source/group_dest are regenerated
  # with these same fixed IDs for whichever attribute is currently selected.
  snap <- restored_pattern_group_snapshot()
  if (!is.null(snap) && identical(snap$attr, input$pattern_attr)) {
    restored_source <- intersect(snap$source, allgroups)
    if (length(restored_source) > 0) {
      selected_groups <- restored_source
      removed_groups  <- allgroups[!(allgroups %in% selected_groups)]
    }
    restored_pattern_group_snapshot(NULL)
  }

  output$ui_sel_order_group <- renderUI({
    tags$div(fluidRow(
      column(7, shinyjqui::orderInput(
        inputId = "group_source",
        label = "Groups to Plot:",
        items = selected_groups,
        width = "100%",
        item_class = "success",
        connect = "group_dest"
      )),
      column(3, shinyjqui::orderInput(
        inputId = "group_dest",
        label = "Drag Here to Remove:",
        items = removed_groups,
        width = "100%",
        placeholder = "Drop items here...",
        item_class = "danger",
        connect = "group_source"
      ))
    )
    )
  })
  
  output$reset_group <- renderUI({
    actionButton("reset_group","Reset")
  })
})

observeEvent(input$reset_group, {
  req(DataQCReactive())
  DataIn = DataQCReactive()
  MetaData=DataIn$MetaData
  allgroups = unique(MetaData[,input$pattern_attr])
  output$ui_sel_order_group <- renderUI({
    tags$div(fluidRow(
      column(7, shinyjqui::orderInput(
        inputId = "group_source",
        label = "Groups to Plot:",
        items = allgroups,
        width = "100%",
        item_class = "success",
        connect = "group_dest"
      )),
      column(3, shinyjqui::orderInput(
        inputId = "group_dest",
        label = "Drag Here to Remove:",
        items = NULL,
        width = "100%",
        placeholder = "Drop items here...",
        item_class = "danger",
        connect = "group_source"
      ))
    )
    )
  })
})

observe({
	DataIn = DataQCReactive()
	ProteinGeneName = DataIn$ProteinGeneName
	results_long = DataIn$tmp_results_long
	if (!is.null(results_long)) {
	  pattern_test = input$pattern_test
	  # "ALL" uses DEGs from every test; otherwise only the selected test
	  if (!is.null(pattern_test) && pattern_test != "ALL") results_long = results_long %>% dplyr::filter(test == pattern_test)
	  pattern_fccut = log2(as.numeric(input$pattern_fccut))
	  pattern_pvalcut = as.numeric(input$pattern_pvalcut)
	  
	  if (input$pattern_psel == "Padj") {
	    filteredgene1 = results_long %>%
	      dplyr::filter(abs(logFC) > pattern_fccut & Adj.P.Value < pattern_pvalcut) %>%
	      dplyr::select(UniqueID) %>% collect %>%	.[["UniqueID"]] %>%	as.character() %>% unique()
	  } else {
	    filteredgene1 = results_long %>%
	      dplyr::filter(abs(logFC) > pattern_fccut & P.Value < pattern_pvalcut) %>%
	      dplyr::select(UniqueID) %>% collect %>%	.[["UniqueID"]] %>%	as.character() %>% unique()
	  }
	  output$patternfilteredgene <- renderText({paste("Selected Genes:",length(filteredgene1),sep="")})
	}
})

filteredgeneReactive <- reactive({
  DataIn <- DataQCReactive()
  ProteinGeneName <- DataIn$ProteinGeneName
  results_long <- DataIn$tmp_results_long
  # "ALL" uses DEGs from every test; otherwise only the selected test
  if (!is.null(results_long) && !is.null(input$pattern_test) && input$pattern_test != "ALL") {
    results_long <- results_long %>% dplyr::filter(test == input$pattern_test)
  }

  pattern_fccut <- log2(as.numeric(input$pattern_fccut))
  pattern_pvalcut <- as.numeric(input$pattern_pvalcut)
  
  if (input$pattern_subset == "subset" && !is.null(results_long)) {
    if (input$pattern_psel == "Padj") {
      return(
        results_long %>%
          dplyr::filter(abs(logFC) > pattern_fccut & Adj.P.Value < pattern_pvalcut) %>%
          dplyr::pull(UniqueID) %>% as.character() %>% unique()
      )
    } else {
      return(
        results_long %>%
          dplyr::filter(abs(logFC) > pattern_fccut & P.Value < pattern_pvalcut) %>%
          dplyr::pull(UniqueID) %>% as.character() %>% unique()
      )
    }
  }
  
  if (input$pattern_subset == "upload genes") {
    pattern_list <- input$pattern_list
    
    if (grepl("\n", pattern_list)) {
      pattern_list <- stringr::str_split(pattern_list, "\n")[[1]]
    } else if (grepl(",", pattern_list)) {
      pattern_list <- stringr::str_split(pattern_list, ",")[[1]]
    }
    
    pattern_list <- gsub(" ", "", pattern_list, fixed = TRUE)
    pattern_list <- unique(pattern_list[pattern_list != ""])

    if (length(pattern_list) <= 2) log_warn("Pattern Clustering: uploaded gene list has ", length(pattern_list), " genes (need > 2).")
    validate(need(length(pattern_list) > 2, "input gene list"))
    
    return(
      ProteinGeneName %>%
        dplyr::filter(
          UniqueID %in% pattern_list |
            Protein.ID %in% pattern_list |
            Gene.Name %in% pattern_list
        ) %>%
        dplyr::pull(UniqueID) %>% as.character()
    )
  }
  character(0)
})

output$selectGroupSamplePattern <- renderUI(shared_header_content())

# Row transform applied to each gene's group means before clustering.
# zscore: centre and scale (trend only); center: subtract the gene's own mean
# (log2 fold change vs the gene's average, keeps the size of the change);
# none: mean log expression as is.
pattern_transform_rows <- function(mat, transform) {
  if (transform == "none") return(mat)
  mat <- sweep(mat, 1, rowMeans(mat), "-")
  if (transform == "zscore") mat <- mat / apply(mat, 1, sd)
  mat
}

pattern_value_label <- function(transform) {
  switch(transform,
         zscore = "Z-score of group mean",
         center = "Centred log2 expression (log2FC vs gene mean)",
         "Mean expression")
}

DatapatternReactive <- eventReactive(input$pattern_plot, {
  req(input$pattern_attr)
  req(input$group_source)
	DataIn = DataQCReactive()
	ProteinGeneName = DataIn$ProteinGeneName
	sel_attr = input$pattern_attr
	sel_group = input$group_source
	transform <- if (is.null(input$pattern_transform)) "zscore" else input$pattern_transform
	data_long <- DataIn$tmp_data_long
	filteredgene <- filteredgeneReactive()
	
	subdatlong <- dplyr::filter(data_long, (.data[[sel_attr]] %in% sel_group) & (UniqueID %in% filteredgene)) %>%
	  group_by(., .data[[sel_attr]], UniqueID)  %>%
	  dplyr::summarise(mean=mean(expr, na.rm = TRUE), .groups = "drop") %>%
	  dplyr::mutate(mean = ifelse(is.nan(mean), NA, mean))

	rawwide <- subdatlong %>%
	  tidyr::pivot_wider(
	    names_from  = all_of(sel_attr),
	    values_from = mean
	  ) %>%
	  as.data.frame() %>%
	  tibble::column_to_rownames("UniqueID") %>%
	  dplyr::select(any_of(as.character(sel_group)))
	rawwide <- as.matrix(rawwide)

	# Genes missing a group mean cannot be clustered (a filled-in 0 would
	# create a fake trend); genes with the same mean in every group have no
	# trend and cannot be scaled.
	n_in <- nrow(rawwide)
	keep_complete <- rowSums(is.na(rawwide)) == 0
	rawwide <- rawwide[keep_complete, , drop = FALSE]
	keep_var <- if (transform == "zscore") apply(rawwide, 1, sd) > 0 else rep(TRUE, nrow(rawwide))
	rawwide <- rawwide[keep_var, , drop = FALSE]
	validate(need(nrow(rawwide) > input$k, sprintf("Only %d genes can be clustered; need more genes than the Cluster Number (%d).", nrow(rawwide), input$k)))

	subdatwide <- as.data.frame(pattern_transform_rows(rawwide, transform))
	subdatlong <- subdatlong %>%
	  dplyr::filter(UniqueID %in% rownames(subdatwide)) %>%
	  dplyr::left_join(
	    subdatwide %>% tibble::rownames_to_column("UniqueID") %>%
	      tidyr::pivot_longer(-UniqueID, names_to = sel_attr, values_to = "value"),
	    by = c("UniqueID", sel_attr)
	  )

	msg <- sprintf("%d genes clustered (%s).", nrow(subdatwide), pattern_value_label(transform))
	if (sum(!keep_complete) > 0) msg <- paste(msg, sprintf("%d genes left out: no expression in at least one group.", sum(!keep_complete)))
	if (sum(!keep_var) > 0) msg <- paste(msg, sprintf("%d genes left out: same mean in every group.", sum(!keep_var)))

	return(list("subdatlong"= subdatlong,"subdatwide"= subdatwide, "filteredgene" = filteredgene,
	            "transform" = transform, "msg" = msg))
})

output$pattern_msg <- renderText({
  req(input$pattern_plot)
  DatapatternReactive()$msg
})

# Cluster once and share the result between the plot and the data table.
pattern_cluster <- eventReactive(input$pattern_plot, {
  Datapattern <- DatapatternReactive()
  subdatwide <- Datapattern$subdatwide
  k <- input$k
  log_info("Pattern Clustering: running ", input$ClusterMethod, " with k=", k, " on ", nrow(subdatwide),
           " genes (", Datapattern$transform, ").")
  set.seed(123)
  if (input$ClusterMethod == "kmeans") {
    cl <- kmeans(subdatwide, k)
    list(method = "kmeans", cl = cl, cluster = cl$cluster)
  } else {
    tmp_expr <- new('ExpressionSet', exprs = as.matrix(subdatwide))
    m1 <- mestimate(tmp_expr)
    cl <- mfuzz(tmp_expr, c = k, m = m1, iter.max = 200)
    list(method = "mfuzz", cl = cl, cluster = cl$cluster, eset = tmp_expr)
  }
})

pattern_out <- eventReactive(input$pattern_plot, {
  withProgress(message = "Processing...", value = 0, {
    req(DatapatternReactive())
    Datapattern <-DatapatternReactive ()
    subdatwide <- Datapattern$subdatwide
    subdatlong <- Datapattern$subdatlong
    sel_attr = input$pattern_attr
    sel_group = input$group_source
    transform <- Datapattern$transform
    clres <- pattern_cluster()
    k=input$k

    if (clres$method == "kmeans") {
      cluster<-clres$cluster
      cluster.df <- data.frame(UniqueID=names(cluster), cluster=paste('Cluster',cluster,sep=' '), row.names=NULL)
      
      subdatlong <- subdatlong  %>%
        left_join(., cluster.df, by="UniqueID")
      
      # Identify the grouping column dynamically
      grp_col <- sel_attr
      
      # Check if the column is categorical
      if (is.character(subdatlong[[sel_attr]]) || is.factor(subdatlong[[sel_attr]])) {
        subdatlong[[sel_attr]] <- factor(subdatlong[[sel_attr]], levels = sel_group)
      }
      
      # Centred and z-scored profiles share one y scale so clusters can be
      # compared; raw means keep a free scale per cluster.
      p <- ggplot(subdatlong, aes(x=.data[[grp_col]], y=value))
      if (transform != "none") p <- p + geom_hline(yintercept = 0, linetype = "dashed", color = "grey60")
      p <- p +
        facet_wrap(~ cluster, scales = if (transform == "none") "free" else "free_x", ncol = input$pattern_ncol) +
        geom_line(aes(group=UniqueID, color="UniqueID")) +
        stat_summary(aes(color="red", group=1), fun=mean, geom="line", size=1.2, group=1)
      p <- p +	theme_bw(base_size = input$pattern_font) + ylab(pattern_value_label(transform)) + xlab(" ") +
        theme (plot.margin = unit(c(1,1,1,1), "cm"), axis.text.x = element_text(angle = input$pattern_Xangle),legend.position="none")
      return(p)
    } else if (clres$method == "mfuzz") {
      # Same content as Mfuzz::mfuzz.plot (genes with membership >= 0.4 in a
      # cluster, coloured by membership with the Mfuzz colour scale), drawn
      # with ggplot so axis labels, font size and x angle work as for K-means.
      min_mem <- 0.4
      mem <- clres$cl$membership
      mem_long <- as.data.frame(mem) %>%
        tibble::rownames_to_column("UniqueID") %>%
        tidyr::pivot_longer(-UniqueID, names_to = "cl", values_to = "membership") %>%
        dplyr::filter(membership >= min_mem) %>%
        dplyr::mutate(cl = as.integer(cl))
      n_per <- mem_long %>% dplyr::count(cl)
      mem_long$cluster <- factor(paste0("Cluster ", mem_long$cl, " (n=", n_per$n[match(mem_long$cl, n_per$cl)], ")"),
                                 levels = paste0("Cluster ", n_per$cl, " (n=", n_per$n, ")"))
      prof <- subdatwide %>% tibble::rownames_to_column("UniqueID") %>%
        tidyr::pivot_longer(-UniqueID, names_to = "grp", values_to = "value") %>%
        dplyr::mutate(grp = factor(grp, levels = colnames(subdatwide)))
      plotdf <- dplyr::inner_join(mem_long, prof, by = "UniqueID", relationship = "many-to-many") %>%
        dplyr::mutate(line_id = paste(UniqueID, cl)) %>%
        dplyr::arrange(membership)
      p <- ggplot(plotdf, aes(x = grp, y = value, group = line_id, color = membership))
      if (transform != "none") p <- p + geom_hline(yintercept = 0, linetype = "dashed", color = "grey60")
      p <- p + geom_line() +
        facet_wrap(~ cluster, scales = if (transform == "none") "free" else "free_x", ncol = input$pattern_ncol, drop = FALSE) +
        scale_color_gradientn(colours = c("#FFFF00", "#00FF00", "#00FFFF", "#0000FF", "#FF00FF", "#FF0000"),
                              limits = c(min_mem, 1), name = "Membership") +
        theme_bw(base_size = input$pattern_font) + ylab(pattern_value_label(transform)) + xlab(" ") +
        theme(plot.margin = unit(c(1,1,1,1), "cm"), axis.text.x = element_text(angle = input$pattern_Xangle))
      return(p)
    }
  })
})

output$pattern<- renderPlot({
	if (input$ClusterMethod == "kmeans") {
		pattern_out()
	} else if (input$ClusterMethod == "pam") {
		pattern_out()
	} else if (input$ClusterMethod == "mfuzz") {
		pattern_out()
	}
})


observeEvent(input$pattern, {
	log_info("Pattern Clustering (", input$ClusterMethod, ") plot saved to output.")
	if (input$ClusterMethod == "kmeans") {
		saved_plots$patternkmeans <- pattern_out()
	} else if (input$ClusterMethod == "pam") {
		saved_plots$patternpam <- pattern_out()
	} else if (input$ClusterMethod == "mfuzz") {
		saved_plots$patternmfuzz <- pattern_out()
	}
})


output$dat_pattern<- DT::renderDataTable({withProgress(message = 'Processing...', value = 0, {
  req(input$pattern_plot)
  Datapattern <-DatapatternReactive()
  subdatwide <- Datapattern$subdatwide
  subdatlong <- Datapattern$subdatlong
  DataIn = DataQCReactive()
  ProteinGeneName = DataIn$ProteinGeneName

  cluster <- pattern_cluster()$cluster
  
  cluster.df <- data.frame(UniqueID = names(cluster), cluster = cluster, row.names = NULL)
  
  if (input$DataFormat == "long") {
    subdatlong <- subdatlong %>% 
      left_join(ProteinGeneName[, c("UniqueID", "Gene.Name")], by = "UniqueID") %>% 
      dplyr::relocate(Gene.Name, .after = UniqueID) %>% 
      left_join(cluster.df, by = "UniqueID") %>% 
      dplyr::mutate(mean = signif(mean, 3), value = signif(value, 3))
    
    DT::datatable(subdatlong)
    
  } else if (input$DataFormat == "wide") {
    subdatwide[, sapply(subdatwide, is.numeric)] <- signif(subdatwide[, sapply(subdatwide, is.numeric)], 3)
    
    subdatwide <- subdatwide %>% 
      rownames_to_column(., var = "UniqueID") %>% 
      left_join(cluster.df, by = "UniqueID") %>% 
      left_join(ProteinGeneName[, c("UniqueID", "Gene.Name")], by = "UniqueID") %>% 
      relocate(Gene.Name, .after = UniqueID) %>% 
      tidyr::separate_wider_delim(
        cols = UniqueID, delim = "_", names = c("Gene", "ID"), 
        too_few = "align_start", cols_remove = TRUE
      ) %>% 
      dplyr::select(where(~!all(is.na(.))))
    
    DT::datatable(subdatwide, extensions = "Buttons", options = list(
      dom = "Blfrtip", 
      buttons = list("copy", list(extend = "collection", buttons = c("csv", "excel", "pdf"), text = "Download")), 
      lengthMenu = list(c(10, 20, -1), c(10, 20, "All")), 
      pageLength = 10
    ), filter = "top")
  }
  })
})



#output$nbclust <- renderPlot({withProgress(message = 'Processing...', value = 0, {
#	Datapattern <-DatapatternReactive ()
#	subdatwide <- Datapattern$subdatwide
#	fviz_nbclust(subdatwide, kmeans, method = "silhouette", k.max = 12) + theme_minimal()
#})
#})
#
