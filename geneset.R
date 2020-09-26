###########################################################################################################
## Proteomics Visualization R Shiny App
##
##This software belongs to Biogen Inc. All right reserved.
##
##@file: geneset.R
##@Developer : Benbo Gao (benbo.gao@Biogen.com)
##@Date : 5/16/2018
##@version 1.0
###########################################################################################################


observe({
	DataIn = DataReactive()
	tests = DataIn$tests
	updateSelectizeInput(session,'geneset_test',choices=tests, selected=tests[1])
})

DataGenesetReactive <- reactive({
	DataIn = DataReactive()
	results_long = DataIn$results_long
	ProteinGeneName = DataIn$ProteinGeneName

	test_sel = input$geneset_test
	absFCcut = as.numeric(input$geneset_FCcut)
	pvalcut = as.numeric(input$geneset_pvalcut)


	if (input$geneset_psel == "Padj") {
		filteredgene = results_long %>%
		dplyr::filter(abs(logFC) >= absFCcut & Adj.P.Value < pvalcut) %>%
		dplyr::filter(test == test_sel) %>%
		dplyr::filter(!is.na(`Gene.Name`)) %>%
		dplyr::select(one_of(c("Gene.Name","logFC"))) %>%
		dplyr::mutate_at(.vars = vars(Gene.Name), .funs = toupper) %>%
		dplyr::distinct(., Gene.Name,.keep_all = TRUE)
	} else {
		filteredgene = results_long %>%
		dplyr::filter(abs(logFC) >= absFCcut & P.Value < pvalcut) %>%
		dplyr::filter(test == test_sel) %>%
		dplyr::filter(!is.na(`Gene.Name`)) %>%
		dplyr::select(one_of(c("Gene.Name","logFC"))) %>%
		dplyr::mutate_at(.vars = vars(Gene.Name), .funs = toupper) %>%
		dplyr::distinct(., Gene.Name,.keep_all = TRUE)

	}

	terminals.df <- dplyr::inner_join(hgnc,filteredgene, by=c("symbol"="Gene.Name"))

	all_genes <- dplyr::filter(ProteinGeneName, !is.na(`Gene.Name`)) %>%
	dplyr::mutate_at(.vars = vars(Gene.Name), funs(toupper)) %>%
	dplyr::select(one_of(c("Gene.Name"))) %>%
	dplyr::inner_join(hgnc,., by=c("symbol"="Gene.Name")) %>%
	dplyr::select(one_of(c("entrez_id"))) %>% collect %>%
	.[["entrez_id"]] %>% as.character() %>% unique()

	sig_genes <- as.numeric(as.data.frame(terminals.df)[,3])
	names(sig_genes) <- as.data.frame(terminals.df)[,2]

	output$geneset_filteredgene <- renderText({ paste("Selected Genes:",nrow(filteredgene), sep="")})
	return(list("sig_genes" = sig_genes,"all_genes" = all_genes, "terminals.df"=terminals.df ))
})

output$MSigDB <- DT::renderDataTable({ withProgress(message = 'Processing...', value = 0, {
	getresults <- DataGenesetReactive()
	sig_genes <- 	getresults$sig_genes
	all_genes <- 	getresults$all_genes

	if (input$MSigDB == "KEGG") {
		gsets = kegg.pathways$human$kg.sets
	} else {
		gsets <- gmtlist[[input$MSigDB]]
	}
	gsa <- ORAEnrichment (deGenes=names(sig_genes),universe=all_genes, gsets, logFC =sig_genes  )
	res <- 	gsa %>%
	#add_rownames(., var = "ID") %>%
	rownames_to_column(var="ID") %>%
	dplyr::filter( p.value < 0.05)

	res[,sapply(res,is.numeric)] <- signif(res[,sapply(res,is.numeric)],3)
	DT::datatable(
		res, selection = 'none', class = 'cell-border strip hover', options = list(pageLength = 15)
	)  %>% formatStyle(1, cursor = 'pointer',color='blue')
	})
})

observeEvent(input$MSigDB_cell_clicked, {
	info = input$MSigDB_cell_clicked
	if (is.null(info$value) || info$col != 1) return()
	updateTabsetPanel(session, 'geneset_tabset', selected = 'Gene Expression')
	updateTextInput(session, 'x1', value = info$value)
	updateTextInput(session, 'x2', value = info$value)
	updateTextInput(session, 'x3', value = info$value)
})

keggView_out <- reactive({withProgress(message = 'Making KEGG Pathway View...', value = 0, {
	validate(need(input$MSigDB == "KEGG", message = "Only works on KEGG."))
	ID = input$x2
	validate(need(ID!="", message = "Select one geneset by clicking geneset name from 'Gene Set Enrichment' tab."))
	getresults <- DataGenesetReactive()
	sig_genes <- 	getresults$sig_genes
	pid <- strsplit(ID," ")[[1]][1]
	img.file <- paste(pid,"pathview","png",sep=".")

	has_img <- file.exists(img.file )
	if (has_img) {
	file.remove(img.file)
	}
	tmp <- pathview(gene.data=sig_genes, pathway.id=pid, kegg.dir="./kegg", kegg.native = T, species="hsa",low = "green", mid = "yellow", high = "red")
	return(img.file)
})
})

observeEvent(input$keggView, {
	saved.num <- length(saved_plots$keggView) +1
	saved_plots$keggView[[saved.num]] <-  keggView_out()
})

output$keggView = renderImage({
	img.file <- keggView_out()
	has_img <- file.exists(img.file )
	if (has_img) {
		list(src = img.file, contentType = 'image/png',	alt = "This is alternate text")
	}
}, deleteFile = FALSE)


output$SetHeatMap = renderPlot({
	ID = input$x3
	validate(need(ID!="", message = "Select one geneset by clicking geneset name from 'Gene Set Enrichment' tab."))

	getresults <- DataGenesetReactive()
	sig_genes <- 	getresults$sig_genes
	all_genes <- 	getresults$all_genes
	terminals.df <- getresults$terminals.df

	DataIn = DataReactive()
	data_long = DataIn$data_long
	ProteinGeneName = DataIn$ProteinGeneName

	if (input$MSigDB == "KEGG") {
		GenesetSig = kegg.pathways$human$kg.sets[[ID]]
	} else {
		GenesetSig <- gmtlist[[input$MSigDB]][[ID]]
	}

	terminalsdf.set <- dplyr::filter(terminals.df, entrez_id %in% GenesetSig)

	terminals_id <- dplyr::filter(ProteinGeneName, Gene.Name %in% terminalsdf.set$symbol)  %>%
	dplyr::select(UniqueID) %>% collect %>% .[["UniqueID"]] %>%	as.character()

	subdatlong <- dplyr::filter(data_long, UniqueID %in% terminals_id ) %>%
	group_by(., group, UniqueID) %>%
	dplyr::summarise(mean=mean(expr, na.rm = TRUE))

	subdatwide <- subdatlong  %>%
	tidyr::spread(.,group, mean, fill = 0) %>%
	as.data.frame() %>%
	remove_rownames(.) %>%
	column_to_rownames(.,var="UniqueID") %>%
	dplyr::select(one_of(as.character(group_order())))

	pheatmap(as.matrix(subdatwide),	scale = "row", color = colorpanel (64, low = "blue",mid = "white", high = "red"),filename=NA)
}, height=800) ## need to change


output$Expression <- DT::renderDataTable({
	ID = input$x1
	validate(need(ID!="", message = "Select one geneset by clicking geneset name from 'Gene Set Enrichment' tab."))

	getresults <- DataGenesetReactive()
	terminals.df <- getresults$terminals.df

	if (input$MSigDB == "KEGG") {
		GenesetSig = kegg.pathways$human$kg.sets[[ID]]
	} else {
		GenesetSig <- gmtlist[[input$MSigDB]][[ID]]
	}

	terminalsdf.set <- dplyr::filter(terminals.df, entrez_id %in% GenesetSig)
	terminalsdf.set[,sapply(terminalsdf.set,is.numeric)] <- signif(terminalsdf.set[,sapply(terminalsdf.set,is.numeric)],3)

	DT::datatable(terminalsdf.set, options = list(pageLength = 15))
})
