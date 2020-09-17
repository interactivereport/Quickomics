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

observe({
  DataIn = DataReactive()
  ProteinGeneName = DataIn$ProteinGeneName
  if (input$network_label=="UniqueID") {
    DataIngenes <- ProteinGeneName %>% dplyr::select(UniqueID) %>% collect %>% .[["UniqueID"]] %>%	as.character()
  } else 
  {DataIngenes <- ProteinGeneName %>% dplyr::select(Gene.Name) %>% collect %>% .[["Gene.Name"]] %>%	as.character()}
  updateSelectizeInput(session,'sel_net_gene', choices= DataIngenes, server=TRUE)
})


observe({
net <-	DataNetworkReactive()

output$networkstat <- renderText({
	sprintf("\nNodes:%d  Edges:%d",
	nrow(net$nodes), nrow(net$edges))
	})

	if (nrow(net$nodes) > 0 & nrow(net$nodes) < 200){
		output$myTabUI <- renderUI({
			actionButton("gennet","Generate")
		})
	} else {
		output$myTabUI <- renderUI({
			"too many nodes or zero node"
		})
	}
})


observeEvent(input$gennet,{
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
#observeEvent(input$gennet,{
#	output$networkD3 <- renderForceNetwork({
#		withProgress(message = 'Making Network:', value = 0, {
#			isolate({
#				net <-	DataNetworkReactive()
#				net$nodes$group = 1
#				net$nodes$size = 10
#				edgelist <- as.data.frame(net$edges)
#				nodes <- as.data.frame(net$nodes)
#				sources <- edgelist$from
#				targets <- edgelist$to
#				node_names <- factor(sort(unique(c(as.character(sources),  as.character(targets)))))
#				links <- data.frame(source = match(sources, node_names) - 1,target = match(targets, node_names) - 1, value = edgelist$cor)
#				nodes <- nodes[match(node_names, nodes$id),]
#				forceNetwork(Links = links, Nodes = nodes, Source = "source",
#					Target = "target", Value = "value", NodeID = "label",
#				Group = "group", opacity = 0.9,zoom = TRUE, opacityNoHover = 1)
#			})
#		})
#	})
#})
#

output$dat_network <- DT::renderDataTable({
	DataIn <- DataReactive()
	net <-	DataNetworkReactive()
	results <- as.data.frame(net$edges)
	results[,sapply(results,is.numeric)] <- signif(results[,sapply(results,is.numeric)],3)
	DT::datatable(results, options = list(pageLength = 15))
})

