###########################################################################################################
## Proteomics Visualization R Shiny App
##
##This software belongs to Biogen Inc. All right reserved.
##
##@file: valcano.R
##@Developer : Benbo Gao (benbo.gao@Biogen.com)
##@Date : 5/16/2018
##@version 1.0
###########################################################################################################

observe({
	DataIn = DataReactive()
	tests = DataIn$tests
	ProteinGeneName = DataIn$ProteinGeneName
	updateRadioButtons(session,'valcano_label', inline = TRUE, choices=colnames(ProteinGeneName)[-1])
	updateSelectizeInput(session,'valcano_test',choices=tests, selected=tests[1])
})


observe({
	DataIn = DataReactive()
	results_long = DataIn$results_long
	test_sel = input$valcano_test
	FCcut = as.numeric(input$valcano_FCcut)
	pvalcut = as.numeric(input$valcano_pvalcut)
	if (input$valcano_psel == "Padj") {
	tmpdat = results_long %>% filter(test==test_sel & Adj.P.Value < pvalcut & abs(logFC) > FCcut) 
	} else {
		tmpdat = results_long %>% filter(test==test_sel & P.Value < pvalcut & abs(logFC) > FCcut) 
	}
	output$valcano_filteredgene <- renderText({paste("Selected Genes:",nrow(tmpdat),sep="")})
})


DataValcanoReactive <- reactive({
	DataIn = DataReactive()
	results_long = DataIn$results_long

	test_sel = input$valcano_test
	FCcut = as.numeric(input$valcano_FCcut)
	pvalcut = as.numeric(input$valcano_pvalcut)
	valcano_label = input$valcano_label
	
	res = results_long %>% filter(test==test_sel) %>%
	filter(!is.na(P.Value)) %>%
	dplyr::mutate (color="Not Significant") %>% as.data.frame() 


	res$labelgeneid = res[,match(valcano_label,colnames(res))]

	if (input$valcano_psel == "Padj") {
	res$color[which((abs(res$logFC)>FCcut)*(res$Adj.P.Value<pvalcut)==1)] = paste0("Padj","<",pvalcut," & abs(logfc)>",FCcut)
	res$color[which((abs(res$logFC)<FCcut)*(res$Adj.P.Value<pvalcut)==1)] =  paste0("Padj","<",pvalcut, " & abs(logfc)<",FCcut)
	res$color = factor(res$color,levels = unique(c("Not Significant",	paste0("Padj","<",pvalcut, " & abs(logfc)<",FCcut),	paste0("Padj","<",pvalcut, " & abs(logfc)>",FCcut))))
	} else { 
	res$color[which((abs(res$logFC)>FCcut)*(res$P.Value<pvalcut)==1)] = paste0("pval","<",pvalcut," & abs(logfc)>",FCcut)
	res$color[which((abs(res$logFC)<FCcut)*(res$P.Value<pvalcut)==1)] =  paste0("pval","<",pvalcut, " & abs(logfc)<",FCcut)
	res$color = factor(res$color,levels = unique(c("Not Significant",	paste0("pval","<",pvalcut, " & abs(logfc)<",FCcut),	paste0("pval","<",pvalcut, " & abs(logfc)>",FCcut))))
	}
	return(res)
})

output$volcanoplot <- renderPlotly({
	res = DataValcanoReactive()
	test_sel = input$valcano_test
	FCcut = as.numeric(input$valcano_FCcut)
	pvalcut = as.numeric(input$valcano_pvalcut)
	if (input$valcano_psel == "Padj") {
		p <- ggplot(res, aes(x = logFC, y =-log10(Adj.P.Value), text=UniqueID))
		ylab <- "-log10(Padj.Value)"
	} else {
		p <- ggplot(res, aes(x = logFC, y =-log10(P.Value), text=UniqueID))
		ylab <- "-log10(P.Value)"
	}
	p <- p	+
	scale_color_manual(values = c("grey", "green2","red2")) +
	geom_point(aes(color = color)) +
	theme_bw(base_size = 20) +
	geom_hline(yintercept = -log10(pvalcut), colour="grey") +
	geom_vline(xintercept = c(-FCcut,0,FCcut), colour="grey") +
	ylab(ylab) + xlab("log2 Fold Change") +
	ggtitle(test_sel)
	p$elementId <- NULL
	p

})

volcanoplotstatic_out <- reactive({
	res = DataValcanoReactive()
	test_sel = input$valcano_test
	FCcut = as.numeric(input$valcano_FCcut)
	pvalcut = as.numeric(input$valcano_pvalcut)

	if (input$valcano_psel == "Padj") {
		p <- ggplot(res, aes(x = logFC, y = -log10(Adj.P.Value)))
		ylab <- "-log10(Padj.Value)"

		filterSig <- paste0("Padj","<",pvalcut," & abs(logfc)>",FCcut)
		data.label <- filter(res, color == filterSig)
		if (nrow(data.label) > 50)
		data.label <- top_n(data.label,50,abs(logFC))

	} else {
		filterSig <- paste0("pval","<",pvalcut," & abs(logfc)>",FCcut)
		data.label <- filter(res, color == filterSig)
		if (nrow(data.label) > 50)
		data.label <- top_n(data.label,50,abs(logFC))

		p <- ggplot(res, aes(x = logFC, y = -log10(P.Value)))
		ylab <- "-log10(P.Value)"
	}
	p <- p	+
	scale_color_manual(values = c("grey", "green2","red2")) +
	geom_point(aes(color = color)) +
	theme_bw(base_size = 20) +
	geom_hline(yintercept = -log10(pvalcut), colour="grey") +
	geom_vline(xintercept = c(-FCcut,0,FCcut), colour="grey") +
	ylab(ylab) + xlab("log2 Fold Change") +
	ggtitle(test_sel) +
	theme(legend.position = "bottom", legend.text=element_text(size=input$yfontsize)) +
	geom_text_repel(data = data.label,  aes(label=labelgeneid),	size = input$lfontsize,	box.padding = unit(0.35, "lines"),	point.padding = unit(0.3, "lines"))
	return(p)
})

output$volcanoplotstatic <- renderPlot({
	volcanoplotstatic_out()
})

observeEvent(input$valcano, {
	test_sel = input$valcano_test
	saved_plots$volcano[[test_sel]] <- volcanoplotstatic_out()
})

output$valcanoData <- DT::renderDataTable({
	DataIn = DataReactive()
	results_long = DataIn$results_long
	test_sel = input$valcano_test
	FCcut = as.numeric(input$valcano_FCcut)
	pvalcut = as.numeric(input$valcano_pvalcut)
	tmpdat = results_long %>% filter(test==test_sel & P.Value < pvalcut & abs(logFC) > FCcut)
	tmpdat[,sapply(tmpdat,is.numeric)] <- signif(tmpdat[,sapply(tmpdat,is.numeric)],3)
	DT::datatable(tmpdat)
})
