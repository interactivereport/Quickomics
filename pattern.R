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
  
  pattern_fccut <- log2(as.numeric(input$pattern_fccut))
  pattern_pvalcut <- as.numeric(input$pattern_pvalcut)
  
  if (input$pattern_subset == "subset" && !is.null(results_long)) {
    if (input$pattern_psel == "Padj") {
      return(
        results_long %>%
          dplyr::filter(abs(logFC) > pattern_fccut & Adj.P.Value < pattern_pvalcut) %>%
          dplyr::pull(UniqueID) %>% as.character()
      )
    } else {
      return(
        results_long %>%
          dplyr::filter(abs(logFC) > pattern_fccut & P.Value < pattern_pvalcut) %>%
          dplyr::pull(UniqueID) %>% as.character()
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


DatapatternReactive <- eventReactive(input$pattern_plot, {
  req(input$pattern_attr)
  req(input$group_source)
	DataIn = DataQCReactive()
	ProteinGeneName = DataIn$ProteinGeneName
	pattern_fccut = log2(as.numeric(input$pattern_fccut))
	pattern_pvalcut = as.numeric(input$pattern_pvalcut)
	sel_attr = input$pattern_attr
	sel_group = input$group_source
	filteredgene=NULL
	results_long <- DataIn$tmp_results_long
	data_long <- DataIn$tmp_data_long
	filteredgene <- filteredgeneReactive()

	subdatlong <- dplyr::filter(data_long, (.data[[sel_attr]] %in% sel_group) & (UniqueID %in% filteredgene)) %>%
	group_by(., .data[[sel_attr]], UniqueID) %>%
	dplyr::summarise(mean=mean(expr, na.rm = TRUE))

	subdatwide <- subdatlong %>%
	  tidyr::pivot_wider(
	    names_from  = sel_attr,
	    values_from = mean,
	    values_fill = 0
	  ) %>%
	  as.data.frame() %>%
	  tibble::column_to_rownames("UniqueID") %>%
	  dplyr::select(all_of(sel_group))
	
	return(list("subdatlong"= subdatlong,"subdatwide"= subdatwide, "filteredgene" = filteredgene))
})

pattern_out <- eventReactive(input$pattern_plot, {
  withProgress(message = "Processing...", value = 0, {
    req(DatapatternReactive())
    Datapattern <-DatapatternReactive ()
    subdatwide <- Datapattern$subdatwide
    subdatlong <- Datapattern$subdatlong
    sel_attr = input$pattern_attr
    sel_group = input$group_source
    
    k=input$k
    set.seed(123)
    if (input$ClusterMethod == "kmeans") {
      cl <- kmeans(subdatwide, k)
      cluster<-cl$cluster
      cluster.df <- data.frame(UniqueID=names(cluster), cluster=paste('Cluster',cluster,sep=' '), row.names=NULL)
      
      subdatlong <- subdatlong  %>%
        left_join(., cluster.df, by="UniqueID")
      
      # subdatlong$group = factor(subdatlong$group,levels = sel_group)
      # 
      # Identify the grouping column dynamically
      grp_col <- sel_attr
      
      # Check if the column is categorical
      if (is.character(subdatlong[[sel_attr]]) || is.factor(subdatlong[[sel_attr]])) {
        subdatlong[[sel_attr]] <- factor(subdatlong[[sel_attr]], levels = sel_group)
      }
      
      p <- ggplot(subdatlong, aes(x=.data[[grp_col]], y=mean)) +
        facet_wrap(~ cluster,scales = "free", ncol = input$pattern_ncol) +
        geom_line(aes(group=UniqueID, color="UniqueID")) +
        stat_summary(aes(color="red", group=1), fun=mean, geom="line", size=1.2, group=1)
      p <- p +	theme_bw(base_size = input$pattern_font) + ylab("Expression") + xlab(" ") +
        theme (plot.margin = unit(c(1,1,1,1), "cm"), axis.text.x = element_text(angle = input$pattern_Xangle),legend.position="none")
      return(p)
      #} else if (input$ClusterMethod == "pam") {
      #	clpam <- pam(subdatwide, k)
      #	cluster <- clpam$clustering
      #
      #	cluster.df <- data.frame(UniqueID=names(cluster), cluster=cluster, row.names=NULL)
      #	subdatlong <- subdatlong  %>%
      #	left_join(., cluster.df, by="UniqueID")
      #	subdatlong$group = factor(subdatlong$group,levels = sel_group)
      #
      #	p <- ggplot(subdatlong, aes(x=group, y=mean)) +
      #	facet_wrap(~ cluster,scales = "free", ncol = 3) +
      #	geom_line(aes(group=UniqueID, color="UniqueID")) +
      #	stat_summary(aes(color="red", group=1), fun=mean, geom="line", size=1.2, group=1)
      #	p <- p + theme_bw(base_size = 14) + ylab("expr") + xlab(" ") +
      #	theme (plot.margin = unit(c(1,1,1,1), "cm"), axis.text.x = element_text(angle = 45),legend.position="none")
      #	return(p)
      #
    } else if (input$ClusterMethod == "mfuzz") {
      tmp_expr <- new('ExpressionSet', exprs = as.matrix(subdatwide))
      m1 <- mestimate(tmp_expr)
      cl <- mfuzz(tmp_expr, c = k, m = m1, iter.max = 200)
      nrow=ceiling(k/input$pattern_ncol)
      mfuzz.plot(tmp_expr, cl = cl, mfrow = c(nrow, input$pattern_ncol), min.mem=0.4, time.labels=colnames(subdatwide),new.window = FALSE)
      p = recordPlot()
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
		replayPlot(pattern_out())
	}
})


observeEvent(input$pattern, {
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
  set.seed(123)
  Datapattern <-DatapatternReactive ()
  subdatwide <- Datapattern$subdatwide
  #subdatwide[,sapply(subdatwide,is.numeric)] <- signif(subdatwide[,sapply(subdatwide,is.numeric)],3)
  
  subdatlong <- Datapattern$subdatlong
  #subdatlong[,sapply(subdatlong,is.numeric)] <- signif(subdatlong[,sapply(subdatlong,is.numeric)],3)

  k=input$k
  
  if (input$ClusterMethod == "kmeans") {
    cl <- kmeans(subdatwide, k)
    cluster<-cl$cluster
    cluster.df <- data.frame(UniqueID=names(cluster), cluster=cluster, row.names=NULL)
    
    if (input$DataFormat == "long") {
      subdatlong <- subdatlong  %>%
        left_join(., cluster.df, by="UniqueID")
      DT::datatable(subdatlong)
    } else if (input$DataFormat == "wide") {
    	subdatwide[,sapply(subdatwide,is.numeric)] <- signif(subdatwide[,sapply(subdatwide,is.numeric)],3)
      subdatwide  <- subdatwide %>%
        rownames_to_column(.,var="UniqueID") %>%
        left_join(., cluster.df, by="UniqueID")%>%
        separate(UniqueID, c("Gene", "ID"), sep = "_")
      
      DT::datatable(subdatwide, 
                    extensions = 'Buttons', 
                    options = list(dom = "Blfrtip", 
                                   buttons = list("copy", list(extend = "collection", 
                                                               buttons = c("csv", "excel", "pdf"), 
                                                               text = "Download")), 
                                   lengthMenu = list( c(10, 20, -1), c(10, 20, "All")), pageLength = 10),
                    filter = 'top')
    }
    
    #} else if (input$ClusterMethod == "pam") {
    #	clpam <- pam(subdatwide, k)
    #	cluster <- clpam$clustering
    #	cluster.df <- data.frame(UniqueID=names(cluster), cluster=cluster, row.names=NULL)
    #
    #	if (input$DataFormat == "long") {
    #		subdatlong <- subdatlong  %>%
    #		left_join(., cluster.df, by="UniqueID")
    #		DT::datatable(subdatlong)
    #
    #	} else if (input$DataFormat == "wide") {
    #		subdatwide  <- subdatwide %>%
    #		rownames_to_column(.,var="UniqueID") %>%
    #		left_join(., cluster.df, by="UniqueID")
    #		DT::datatable(subdatwide, filter = 'top')
    #	}
    #
  } else if (input$ClusterMethod == "mfuzz") {
    tmp_expr <- new('ExpressionSet', exprs = as.matrix(subdatwide))
    m1 <- mestimate(tmp_expr)
    cl <- mfuzz(tmp_expr, c = k, m = m1, iter.max = 200)
    cluster <- cl$cluster
    cluster.df <- data.frame(UniqueID=names(cluster), cluster=cluster, row.names=NULL)
    if (input$DataFormat == "long") {
      #subdatlong[,sapply(subdatlong,is.numeric)] <- signif(subdatlong[,sapply(subdatlong,is.numeric)],3)
      subdatlong <- subdatlong  %>%
        left_join(., cluster.df, by="UniqueID")
      DT::datatable(subdatlong)
      
    } else if (input$DataFormat == "wide") {
      subdatwide[,sapply(subdatwide,is.numeric)] <- signif(subdatwide[,sapply(subdatwide,is.numeric)],3)
      subdatwide  <- subdatwide %>%
        rownames_to_column(.,var="UniqueID") %>%
        left_join(., cluster.df, by="UniqueID") %>%
        separate(UniqueID, c("Gene", "ID"), sep = "_")
      
      DT::datatable(subdatwide, 
                    extensions = 'Buttons', 
                    options = list(dom = "Blfrtip", 
                                   buttons = list("copy", list(extend = "collection", 
                                                               buttons = c("csv", "excel", "pdf"), 
                                                               text = "Download")), 
                                   lengthMenu = list( c(10, 20, -1), c(10, 20, "All")), pageLength = 10),
      filter = 'top')
    }
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
