###########################################################################################################
## Proteomics Visualization R Shiny App
##
##This software belongs to Biogen Inc. All right reserved.
##
##@file: scurve.R
##@Developer : Benbo Gao (benbo.gao@Biogen.com)
##@Date : 2/15/2021
##@version 1.0
###########################################################################################################
observe({
	DataIn = DataReactive()
	ProteinGeneName = DataIn$data_results
	DataIngenes <- ProteinGeneName %>% dplyr::select(UniqueID) %>% collect %>% .[["UniqueID"]] %>%	as.character()
	updateSelectizeInput(session,'sel_gene_scurve', choices= DataIngenes, server=TRUE)
})


Scurve_out <- reactive({
	DataIn = DataReactive()
	data_results <- DataIn$data_results

	gene_list <- input$gene_list_scurve

	if (grepl("\n", gene_list)) {
		gene_list <-  stringr::str_split(gene_list, "\n")[[1]]
	} else if (grepl(",", gene_list)) {
		gene_list <-  stringr::str_split(genen_list, ",")[[1]]
	}
	gene_list <- gsub(" ", "", gene_list, fixed = TRUE)
	gene_list <- unique(gene_list[gene_list != ""])

	if (length(input$sel_gene_scurve) > 0) {
		gene_list <- input$sel_gene_scurve
	}

	validate(need(length(gene_list)>0,"Please select a gene or input gene."))

	scurve.data <- data_results %>%
	dplyr::select(one_of(c("UniqueID","Gene.Name","Intensity","Protein.ID"))) %>%
	dplyr::filter((!is.na(Intensity)) & Intensity > 0) %>%
	dplyr::arrange(desc(Intensity)) %>%
	dplyr::mutate(RANK = row_number())

	scurve.label <- scurve.data %>%
	dplyr::filter((UniqueID %in% gene_list) | (Gene.Name %in% gene_list) | (Protein.ID %in% gene_list))

	if (nrow(scurve.label) > 0) {
		p <- ggplot(scurve.data, aes(x = RANK, y = Intensity)) +
		geom_point(colour = "blue")  +
		geom_point(data=scurve.label, aes(x = RANK, y = Intensity, colour='red')) +
		scale_y_continuous(trans='log10')  +
		geom_label_repel(
			data = scurve.label,aes(x = RANK, y = Intensity,  label = Gene.Name),
			fontface = 'bold', color = 'red', size = input$scurve_labelfontsize,
			box.padding = unit(0.35, "lines"),
			point.padding = unit(0.5, "lines"),
			segment.color = 'red'
		) +
		theme_bw(base_size = 14) + 
		ylab(input$scurveYlab) + 
		xlab(input$scurveXlab) +
		theme (plot.margin = unit(c(1,1,1,1), "cm"),
			text = element_text(size=input$scurve_axisfontsize),
			axis.text.x = element_text(angle = input$scurveXangle),
		legend.position="none")

		p <- ggExtra::ggMarginal(p, type = input$scurveright, margins = "y",colour="blue", fill = "gray")

	}
	return(p)
})

output$SCurve <- renderPlot({
	Scurve_out()
})





