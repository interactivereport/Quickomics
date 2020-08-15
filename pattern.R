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

observe({
	DataIn = DataReactive()
	tests = c("ALL",DataIn$tests)
	allgroups = DataIn$groups
	groups = group_order()
	updateSelectizeInput(session,'pattern_group', choices=allgroups, selected=groups)
	updateSelectizeInput(session,'pattern_test',choices=tests, selected=tests[1])
})

observe({
	DataIn = DataReactive()
	results_long = DataIn$results_long
	pattern_test = input$pattern_test
	pattern_fccut = as.numeric(input$pattern_fccut)
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
})

DatapatternReactive <- reactive({
	DataIn = DataReactive()
	sample_group <- DataIn$sample_group
	pattern_fccut = as.numeric(input$pattern_fccut)
	pattern_pvalcut = as.numeric(input$pattern_pvalcut)
	sel_group = input$pattern_group

	group_order(input$pattern_group)

	results_long <- DataIn$results_long
	data_long <- DataIn$data_long

	if (input$pattern_psel == "Padj") {
		filteredgene = results_long %>%
		dplyr::filter(abs(logFC) > pattern_fccut & Adj.P.Value < pattern_pvalcut) %>%
		dplyr::select(UniqueID) %>% collect %>%	.[["UniqueID"]] %>%	as.character()

	} else {
		filteredgene = results_long %>%
		dplyr::filter(abs(logFC) > pattern_fccut & P.Value < pattern_pvalcut) %>%
		dplyr::select(UniqueID) %>% collect %>%	.[["UniqueID"]] %>%	as.character()

	}

	subdatlong <- dplyr::filter(data_long, (group %in% sel_group) & (UniqueID %in% filteredgene)) %>%
	group_by(., group, UniqueID) %>%
	dplyr::summarise(mean=mean(expr, na.rm = TRUE))

	subdatwide <- subdatlong  %>%
	tidyr::spread(.,group, mean, fill = 0) %>% as.data.frame() %>%
	remove_rownames(.) %>%
	column_to_rownames(.,var="UniqueID") %>%
	dplyr::select(sel_group)

	return(list("subdatlong"= subdatlong,"subdatwide"= subdatwide, "filteredgene" = filteredgene))
})

pattern_out <- reactive({withProgress(message = 'Processing...', value = 0, {

	Datapattern <-DatapatternReactive ()
	subdatwide <- Datapattern$subdatwide
	subdatlong <- Datapattern$subdatlong
	sel_group=input$pattern_group

	k=input$k
	set.seed(123)
	if (input$ClusterMehtod == "kmeans") {

		cl <- kmeans(subdatwide, k)
		cluster<-cl$cluster

		cluster.df <- data.frame(UniqueID=names(cluster), cluster=paste('Cluster',cluster,sep=' '), row.names=NULL)

		subdatlong <- subdatlong  %>%
		left_join(., cluster.df, by="UniqueID")

		subdatlong$group = factor(subdatlong$group,levels = sel_group)

		p <- ggplot(subdatlong, aes(x=group, y=mean)) +
		facet_wrap(~ cluster,scales = "free", ncol = 3) +
		geom_line(aes(group=UniqueID, color="UniqueID")) +
		stat_summary(aes(color="red", group=1), fun.y=mean, geom="line", size=1.2, group=1)
		p <- p +	theme_bw(base_size = input$pattern_font) + ylab("expr") + xlab(" ") +
		theme (plot.margin = unit(c(1,1,1,1), "cm"), axis.text.x = element_text(angle = input$pattern_Xangle),legend.position="none")
		return(p)
	#} else if (input$ClusterMehtod == "pam") {
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
	#	stat_summary(aes(color="red", group=1), fun.y=mean, geom="line", size=1.2, group=1)
	#	p <- p + theme_bw(base_size = 14) + ylab("expr") + xlab(" ") +
	#	theme (plot.margin = unit(c(1,1,1,1), "cm"), axis.text.x = element_text(angle = 45),legend.position="none")
	#	return(p)
  #
	} else if (input$ClusterMehtod == "mfuzz") {
		tmp_expr <- new('ExpressionSet', exprs = as.matrix(subdatwide))
		m1 <- mestimate(tmp_expr)
		cl <- mfuzz(tmp_expr, c = k, m = m1, iter.max = 200)
		if (k <= 6) {
			nrow = 2
		} else if (k <= 9) {
			nrow = 3
		} else if (k <= 12) {
			nrow = 4
		}
		mfuzz.plot(tmp_expr, cl = cl, mfrow = c(nrow, 3), min.mem=0.4, time.labels=colnames(subdatwide),new.window = FALSE)
		p = recordPlot()
		return(p)
	}

})
})

output$pattern<- renderPlot({
	if (input$ClusterMehtod == "kmeans") {
		pattern_out()
	} else if (input$ClusterMehtod == "pam") {
		pattern_out()
	} else if (input$ClusterMehtod == "mfuzz") {
		replayPlot(pattern_out())
	}
})


observeEvent(input$pattern, {
	if (input$ClusterMehtod == "kmeans") {
		saved_plots$patternkmeans <- pattern_out()
	} else if (input$ClusterMehtod == "pam") {
		saved_plots$patternpam <- pattern_out()
	} else if (input$ClusterMehtod == "mfuzz") {
		saved_plots$patternmfuzz <- pattern_out()
	}
})


output$dat_pattern<- DT::renderDataTable({withProgress(message = 'Processing...', value = 0, {
  set.seed(123)
  
	
  Datapattern <-DatapatternReactive ()
  subdatwide <- Datapattern$subdatwide
  #subdatwide[,sapply(subdatwide,is.numeric)] <- signif(subdatwide[,sapply(subdatwide,is.numeric)],3)
  
  subdatlong <- Datapattern$subdatlong
  #subdatlong[,sapply(subdatlong,is.numeric)] <- signif(subdatlong[,sapply(subdatlong,is.numeric)],3)
  
  sel_group=input$pattern_group
  
  k=input$k
  
  if (input$ClusterMehtod == "kmeans") {
    cl <- kmeans(subdatwide, k)
    cluster<-cl$cluster
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
    
    #} else if (input$ClusterMehtod == "pam") {
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
  } else if (input$ClusterMehtod == "mfuzz") {
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