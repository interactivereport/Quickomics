###########################################################################################################
## Proteomics Visualization R Shiny App
##
##This software belongs to Biogen Inc. All right reserved.
##
##@file: network.R
##@Developer : Benbo Gao (benbo.gao@Biogen.com)
##@Date : 5/16/2018
##@version 1.0
###########################################################################################################
current_samples <- reactiveVal()
current_nw_genes <- reactiveVal()
current_network <- reactiveVal()

compute_network <- function(data_wide) {
  # Trim to 10k genes
  if (nrow(data_wide) > 10000) {
    dataSD <- apply(data_wide, 1, sd, na.rm = TRUE)
    dataM  <- rowMeans(data_wide)
    diff   <- dataSD / (dataM + median(dataM))
    data_wide <- data_wide[order(diff, decreasing = TRUE)[1:10000], ]
  }
  
  cor_res <- Hmisc::rcorr(as.matrix(t(data_wide)))
  cormat  <- cor_res$r
  pmat    <- cor_res$P
  ut      <- upper.tri(cormat)
  
  network <- tibble(
    from = rownames(cormat)[row(cormat)[ut]],
    to   = rownames(cormat)[col(cormat)[ut]],
    cor  = signif(cormat[ut], 2),
    p    = signif(pmat[ut], 2),
    direction = as.integer(sign(cormat[ut]))
  ) %>%
    mutate_if(is.factor, as.character) %>%
    filter(!is.na(cor) & abs(cor) > 0.7 & p < 0.05)
  
  # Additional filtering if too large
  if (nrow(network) > 2e6) {
    network <- network %>% filter(abs(cor) > 0.8 & p < 0.005)
  }
  if (nrow(network) > 2e6) {
    network <- network %>% filter(abs(cor) > 0.85 & p < 0.005)
  }
  return(network)
}

edges_UniqueID_GeneName_mapping <- function(edges, ProteinGeneName) {
  ProteinGeneName <- data.table::as.data.table(ProteinGeneName)
  edges <- data.table::as.data.table(edges)
  data.table::setkey(ProteinGeneName, UniqueID)
  
  edges[, from := ProteinGeneName[edges$from, Gene.Name]]
  edges[, to   := ProteinGeneName[edges$to, Gene.Name]]
  ProteinGeneName <- as.data.frame(ProteinGeneName)
  edges <- as.data.frame(edges)
  return(edges)
}

# Last "Select Gene Label" value seen by the observer below, so a real
# label switch (Gene.Name <-> UniqueID) can be told apart from the first
# run or a recomputation, and reset the selection and results.
last_network_label <- reactiveVal(NULL)

observeEvent(list(NetworkReactive(), input$network_label, restored_sel_net_gene()), {
  # restored_sel_net_gene() is included as a trigger (not just read via
  # isolate() below) specifically because onRestored() is not guaranteed to
  # run BEFORE this observer's first natural firing (e.g. when Save Session
  # also restores input$menu to "Correlation_Network" directly, this
  # observer's first run can beat onRestored() to the punch and find nothing
  # pending yet). Reacting to it too means whichever happens last -- the tab
  # visit or the restore -- still triggers a correct re-run.
  all_genes <- current_nw_genes()
  if (input$network_label == "UniqueID") {
    DataIngenes <- all_genes %>%	as.character()
  } else {
    DataIngenes <- DataQCReactive()$ProteinGeneName %>%
      filter(UniqueID %in% all_genes) %>%
      pull(Gene.Name) %>%
      na.omit() %>%
      discard(~ .x == "") %>%
      unique() %>%
      as.character()
  }
  selected_genes <- intersect(input$sel_net_gene, DataIngenes)
  prev_label <- isolate(last_network_label())
  last_network_label(input$network_label)
  # If a Save Session restore is waiting on this tab's choices to exist
  # (see restored_sel_net_gene in app.R), apply it now instead of whatever
  # input$sel_net_gene currently holds. It stays pending (cleared by the
  # observer below once the browser reports a selection) because this
  # observer can re-run during restore before input$sel_net_gene has come
  # back, and consuming it on the first run let that re-run wipe it.
  pending_restore <- isolate(restored_sel_net_gene())
  if (!is.null(pending_restore)) {
    selected_genes <- intersect(pending_restore, DataIngenes)
  }
  if (!is.null(prev_label) && !identical(prev_label, input$network_label)) {
    selected_genes <- character(0)
    restored_sel_net_gene(NULL)
    output$networkstat <- renderText(NULL)
    output$myTabUI <- renderUI(NULL)
  }
  updateSelectizeInput(session,'sel_net_gene', choices= DataIngenes, selected = selected_genes, server=TRUE)
  output$visnetwork <- renderVisNetwork({
    NULL  
  })
  output$networkD3 <- renderForceNetwork({
    NULL
  })
  # output$dat_network <- DT::renderDT(server=FALSE,{
  #   NULL
  # })
})

# Consume a pending Save Session restore of sel_net_gene once the browser
# reports a selection (the restored one, or the user's own).
observeEvent(input$sel_net_gene, {
  if (!is.null(restored_sel_net_gene())) restored_sel_net_gene(NULL)
})

NetworkReactive <- reactive({
  req(input$menu == "Correlation_Network")
  # Evaluate DataQCReactive() outside sort(): sort() is an S4 generic here
  # (BiocGenerics), and S4 dispatch turns a req() failure inside its argument
  # into a real error, which ends the session. That happened when restoring a
  # session saved on this tab, before Groups and Samples was ready.
  DataQC <- DataQCReactive()
  req(DataQC)
  if (!identical(sort(as.character(DataQC$tmp_sampleid)), current_samples())) {
    isolate({
      n_sample_all <- length(all_samples())
      n_sample_current <- length(DataQCReactive()$tmp_sampleid)
      ProteinGeneName <- DataQCReactive()$ProteinGeneName
      
      Pinfo <- ProjectInfo
      CorResFile <- Pinfo$file2
      if (is.null(CorResFile)) {
        CorResFile <- paste0("networkdata/", Pinfo$ProjectID, "_network.RData")
      }
      # API-sourced projects (ProjectInfo$source set to "published"/
      # "unpublished"/"datalake") have no local network storage to write a
      # cache file back to, so the write-back below is intentionally skipped
      # for them. A 404 from the API (no precomputed network for this project)
      # is treated the same as a missing local file -- network stays NULL
      # and falls through to the compute-fresh branch below, no error.
      is_api_sourced <- !is.null(Pinfo$source)
      network <- NULL

      if (n_sample_current == n_sample_all) {
        if (is_api_sourced) {
          log_info("Correlation Network: loading cached network through the Quickomics API (source=", Pinfo$source, ")")
          network <- api_load_network(Pinfo$ProjectID, Pinfo$source)
        } else if (file.exists(CorResFile)) {
          log_info("Correlation Network: loading cached network from ", CorResFile)
          load(CorResFile)
          stopifnot(exists("network"))
        }
      }

      if (is.null(network)) {
        log_info("Correlation Network: computing correlation network for ", n_sample_current, " samples (no cache available).")
        data_wide = DataQCReactive()$tmp_data_wide
        withProgress(message = "Compute correlation network data.",
                     detail = "This may take a few minutes...",
                     value = 0, {
                       network <- compute_network(data_wide)
                     })
        if (n_sample_current == n_sample_all && !is_api_sourced) {
          save(network, file = CorResFile)
          ProjectInfo$file2 <- CorResFile
        }
      }
      log_info("Correlation Network: ", nrow(network), " edges available.")
      current_samples(sort(as.character(DataQCReactive()$tmp_sampleid)))
      current_network(network)
      current_nw_genes(unique(c(network$from, network$to)))
    })
  }
})
output$selectGroupSampleNetwork <- renderUI(shared_header_content())

DataNetworkReactive <- eventReactive(list(input$sel_net_gene, input$network_rcut, input$network_pcut, sort(as.character(current_samples()))), {
  req(input$sel_net_gene)
  network <- current_network()
  # DataIn  <- DataReactive()
  # ProteinGeneName <- DataIn$ProteinGeneName
  DataIn <- DataQCReactive()
  ProteinGeneName <- DataIn$ProteinGeneName
  req(network, ProteinGeneName)
  sel_gene <- input$sel_net_gene
  
  # Map selected gene to UniqueID
  tmpids <- ProteinGeneName[
    unique(na.omit(c(
      apply(ProteinGeneName, 2, function(k) match(sel_gene, k))
    ))),
  ]
  
  # Filter edges involving selected gene
  edges.sel <- network %>%
    filter(from %in% tmpids$UniqueID | to %in% tmpids$UniqueID)
  
  # Apply user cutoffs
  rcutoff     <- as.numeric(input$network_rcut)
  pvalcutoff  <- as.numeric(input$network_pcut)
  
  edges <- edges.sel %>%
    filter(abs(cor) > rcutoff & p < pvalcutoff)
  
  # Build node list
  networks_ids <- unique(c(edges$from, edges$to))
  
  nodes <- ProteinGeneName %>%
    dplyr::filter(UniqueID %in% networks_ids) %>%
    dplyr::select(UniqueID, Gene.Name) %>%
    dplyr::rename(id = UniqueID, label = Gene.Name)
  
  list(nodes = nodes, edges = edges)
})


observe({
  net <-	DataNetworkReactive()
  req(net)
  output$networkstat <- renderText({
    sprintf("\nNodes:%d  Edges:%d",
            nrow(net$nodes), nrow(net$edges))
	})

	if (nrow(net$nodes) > 0 & nrow(net$nodes) <= 200){
		output$myTabUI <- renderUI({
			actionButton("gennet","Generate")
		})
	} else if (nrow(net$nodes) == 0) {
		log_warn("Correlation Network: current cutoffs/selection produced zero nodes.")
		output$myTabUI <- renderUI({
			"Zero node. Try lower r cutoff, higher P Value cutoff or select other genes."
		})
	} else {
	  log_warn("Correlation Network: current cutoffs/selection produced ", nrow(net$nodes), " nodes (over the 200 limit).")
	  output$myTabUI <- renderUI({
	    "Too many nodes (1 ~ 200 are acceptable). Try higher r cutoffs, , lower P Value cutoff or select fewer genes."
	  })
	}
})


observeEvent(input$gennet,{
	log_info("Correlation Network: Generate clicked (visNetwork).")
	output$visnetwork <- renderVisNetwork({
		withProgress(message = 'Making Network:', value = 0, {
			isolate({
				net <-	DataNetworkReactive()
				visNetwork(net$nodes,net$edges,  height = "800px", width = "100%") %>%
				visLayout(randomSeed = 123) %>%
				visInteraction(navigationButtons = TRUE)

			})
		})
	})
})


#
observeEvent(input$gennet,{
	output$networkD3 <- renderForceNetwork({
		withProgress(message = 'Making Network:', value = 0, {
			isolate({
				net <-	DataNetworkReactive()
				net$nodes$group = 1
				net$nodes$size = 10
				edgelist <- as.data.frame(net$edges)
				nodes <- as.data.frame(net$nodes)
				sources <- edgelist$from
				targets <- edgelist$to
				node_names <- factor(sort(unique(c(as.character(sources),  as.character(targets)))))
				links <- data.frame(source = match(sources, node_names) - 1,target = match(targets, node_names) - 1, value = edgelist$cor)
				nodes <- nodes[match(node_names, nodes$id),]
				forceNetwork(Links = links, Nodes = nodes, Source = "source",
					Target = "target", Value = "value", NodeID = "label",fontSize=7,
				Group = "group", opacity = 0.9,zoom = TRUE, opacityNoHover = 1)
			})
		})
	})
})


output$dat_network <- DT::renderDT(server=FALSE,{
	# DataIn <- DataReactive()
	net <-	DataNetworkReactive()
	results <- as.data.frame(net$edges)
	if (input$network_label == "Gene.Name") {
	  results <- edges_UniqueID_GeneName_mapping(results, DataQCReactive()$ProteinGeneName)
	}
	results[,sapply(results,is.numeric)] <- signif(results[,sapply(results,is.numeric)],3)
	DT::datatable(results, extensions = 'Buttons', options = list(dom = 'lBfrtip', pageLength = 15,
	     buttons = list(
	     list(extend = "csv", text = "Download Page", filename = "Page_results",
	          exportOptions = list(modifier = list(page = "current"))),
	     list(extend = "csv", text = "Download All", filename = "All_Results",
	          exportOptions = list(modifier = list(page = "all")))
	   )
	))
})

