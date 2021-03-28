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
	tests_more=c("None", tests)
	updateSelectizeInput(session,'geneset_test2',choices=tests_more, selected="None")
	updateSelectizeInput(session,'geneset_test3',choices=tests_more, selected="None")
	updateSelectizeInput(session,'geneset_test4',choices=tests_more, selected="None")
	updateSelectizeInput(session,'geneset_test5',choices=tests_more, selected="None")
})

DataGenesetReactive <- reactive({
	DataIn = DataReactive()
	results_long = DataIn$results_long
	ProteinGeneName = DataIn$ProteinGeneName

	test_sel = input$geneset_test
	absFCcut = log2(as.numeric(input$geneset_FCcut))
	pvalcut = as.numeric(input$geneset_pvalcut)
	filteredgene <- results_long %>%  mutate(P.stat=ifelse(rep(input$geneset_psel == "Padj", nrow(results_long)),  Adj.P.Value,  P.Value)) %>%
	dplyr::filter(test == test_sel) %>%dplyr::filter(!is.na(`Gene.Name`))%>%
	dplyr::filter(abs(logFC) >= absFCcut & P.stat < pvalcut) %>%
	dplyr::select(one_of(c("Gene.Name","logFC"))) %>%
	dplyr::distinct(., Gene.Name,.keep_all = TRUE)
	#browser()#bebug
	if (ProjectInfo$Species=="human") {
		filteredgene<-filteredgene%>%dplyr::mutate_at(.vars = vars(Gene.Name), .funs = toupper)
		terminals.df <- dplyr::inner_join(hgnc,filteredgene, by=c("symbol"="Gene.Name"))
		all_genes <- dplyr::filter(ProteinGeneName, !is.na(`Gene.Name`)) %>%
		dplyr::mutate(Gene.Name=toupper(Gene.Name)) %>%
		dplyr::select(one_of(c("Gene.Name"))) %>%
		dplyr::inner_join(hgnc,., by=c("symbol"="Gene.Name")) %>%
		dplyr::select(one_of(c("entrez_id"))) %>% collect %>%
		.[["entrez_id"]] %>% as.character() %>% unique()
	} else if  (ProjectInfo$Species=="mouse") {
		terminals.df<-filteredgene%>%transmute(mouse_Gene=Gene.Name, logFC)%>%left_join(M_match%>%transmute(mouse_Gene=mouse_symbol, symbol=human_symbol))%>%
		mutate(symbol=ifelse(is.na(symbol), toupper(mouse_Gene), symbol))%>%left_join(hgnc)%>%filter(!is.na(entrez_id))%>%
		dplyr::select(symbol, entrez_id, logFC, mouse_Gene)
		all_genes<-ProteinGeneName%>%left_join(M_match%>%transmute(Gene.Name=mouse_symbol, symbol=human_symbol))%>%
		mutate(symbol=ifelse(is.na(symbol), toupper(Gene.Name), symbol))%>%left_join(hgnc)%>%filter(!is.na(entrez_id))%>%
		dplyr::select(one_of(c("entrez_id"))) %>% collect %>%.[["entrez_id"]] %>% as.character() %>% unique()
	} else if  (ProjectInfo$Species=="rat") {
		terminals.df<-filteredgene%>%transmute(rat_Gene=Gene.Name, logFC)%>%left_join(R_match%>%transmute(rat_Gene=rat_symbol, symbol=human_symbol))%>%
		mutate(symbol=ifelse(is.na(symbol), toupper(rat_Gene), symbol))%>%left_join(hgnc)%>%filter(!is.na(entrez_id))%>%
		dplyr::select(symbol, entrez_id, logFC, rat_Gene)
		all_genes<-ProteinGeneName%>%left_join(R_match%>%transmute(Gene.Name=rat_symbol, symbol=human_symbol))%>%
		mutate(symbol=ifelse(is.na(symbol), toupper(Gene.Name), symbol))%>%left_join(hgnc)%>%filter(!is.na(entrez_id))%>%
		dplyr::select(one_of(c("entrez_id"))) %>% collect %>%.[["entrez_id"]] %>% as.character() %>% unique()
	}


	sig_genes <- as.numeric(as.data.frame(terminals.df)[,3])
	names(sig_genes) <- as.data.frame(terminals.df)[,2]

	sig_genes_Dir=sig_genes
	if (input$geneset_direction=="Up") {
		sig_genes_Dir=sig_genes[terminals.df$logFC>0]
	} else if (input$geneset_direction=="Down") {
		sig_genes_Dir=sig_genes[terminals.df$logFC<0]
	}
	#browser() #debug

	match_info=ifelse(ProjectInfo$Species=="mouse", "(mouse genes mapped to human)", ifelse (ProjectInfo$Species=="rat", "(rat genes mapped to human)", ""))
	output$geneset_filteredgene <- renderText({ paste("Selected Genes", match_info, ":",length(sig_genes_Dir), " (",input$geneset_direction, " Direction)",  sep="")})
	#	browser()#debug
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
	gsa <- ORAEnrichment (deGenes=names(sig_genes),universe=all_genes, gsets, logFC =sig_genes, Dir=input$geneset_direction  )
	res <- 	gsa %>%
	#add_rownames(., var = "ID") %>%
	rownames_to_column(var="ID") %>%
	dplyr::filter( p.value < 0.05)

	res[,sapply(res,is.numeric)] <- signif(res[,sapply(res,is.numeric)],3)
	DT::datatable(
		res,  extensions = 'Buttons', selection = 'none', class = 'cell-border strip hover', 
		options = list( dom = 'lBfrtip', buttons = c('csv', 'excel', 'print'), pageLength = 15)
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
	dataIn=DataReactive()
	data_results=dataIn$data_results
	
	if (ProjectInfo$Species=="human") {
		data_results<-data_results%>%left_join(hgnc, by=c("Gene.Name"="symbol"))
	} else if  (ProjectInfo$Species=="mouse") {
		data_results<-data_results%>%mutate(mouse_Gene=Gene.Name)%>%left_join(M_match%>%transmute(mouse_Gene=mouse_symbol, symbol=human_symbol))%>%
		mutate(symbol=ifelse(is.na(symbol), toupper(mouse_Gene), symbol))%>%left_join(hgnc)
	} else if  (ProjectInfo$Species=="rat") {
		data_results<-data_results%>%mutate(rat_Gene=Gene.Name)%>%left_join(R_match%>%transmute(rat_Gene=rat_symbol, symbol=human_symbol))%>%
		mutate(symbol=ifelse(is.na(symbol), toupper(rat_Gene), symbol))%>%left_join(hgnc)
	}

	if (file.exists(img.file)) {
		file.remove(img.file)
	}
	#get logFC data from data_results
	tests=input$geneset_test
	if (input$kegg_more_tests=="Yes") {
		if (input$geneset_test2!="None") {tests=c(tests, input$geneset_test2)}
		if (input$geneset_test3!="None") {tests=c(tests, input$geneset_test3)}
		if (input$geneset_test4!="None") {tests=c(tests, input$geneset_test4)}
		if (input$geneset_test5!="None") {tests=c(tests, input$geneset_test5)}
	}
	selCol=rep(NA, length(tests))
	all_names=names(data_results)
	for (i in 1:length(tests)) {
		sel_i=which(str_detect(all_names, regex(str_c("^", tests[i]), ignore_case=T)) & str_detect(all_names, regex("logFC$", ignore_case=T)) )
		if (length(sel_i)==1) {selCol[i]=sel_i}
	}
	#selCol=match(str_c(tests, ".logFC"), names(data_results) )
	if (sum(is.na(selCol))==0) {
		sel_gene=which(!is.na(data_results$entrez_id))
		FCdata=data.matrix(data_results[sel_gene, selCol])
		rownames(FCdata)=data_results$entrez_id[sel_gene]
		tmp <- pathview(gene.data=FCdata, pathway.id=pid, kegg.dir="./kegg", kegg.native = T, species="hsa",low = "green", mid = "yellow", high = "red",
		same.layer = F, map.symbol=as.logical(input$kegg_mapsample), limit=list(gene=as.numeric(input$kegg_logFC), cpd=1) )
			if (ncol(FCdata)>1) {	img.file <- paste(pid,"pathview.multi.png",sep=".")}
		} else {
			cat("looking for additional comparisons failed", tests, selCol, "\ngo back to the first comparison\n")
			#browser()#debug
			tmp <- pathview(gene.data=sig_genes, pathway.id=pid, kegg.dir="./kegg", kegg.native = T, species="hsa",low = "green", mid = "yellow", high = "red",
			same.layer = F, map.symbol=as.logical(input$kegg_mapsample), limit=list(gene=as.numeric(input$kegg_logFC), cpd=1) )
				img.file <- paste(pid,"pathview","png",sep=".")
			}
			#browser() #debug
			return(img.file)
		})
	})


observeEvent(input$keggSave, {
	ID = input$x3
	img.file <- keggView_out()
	if (file.exists(img.file)) {
		img <- readPNG(img.file)
		saved_plots$keggSave[[ID]] <- img
	}
})

output$keggView = renderImage({
	img.file <- keggView_out()
	if (file.exists(img.file)) {
		list(src = img.file, contentType = 'image/png',	alt = "This is alternate text")
	}
}, deleteFile = FALSE)

genesetheatmap_out <- reactive({withProgress(message = 'Making heatmap...', value = 0, {
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

	if (ProjectInfo$Species=="human") {
		terminals_id <- dplyr::filter(ProteinGeneName, toupper(Gene.Name) %in% terminalsdf.set$symbol)  %>%
		dplyr::select(UniqueID) %>% collect %>% .[["UniqueID"]] %>%	as.character()
	} else if (ProjectInfo$Species=="mouse"){
		terminals_id <- dplyr::filter(ProteinGeneName, Gene.Name %in% terminalsdf.set$mouse_Gene)  %>%
		dplyr::select(UniqueID) %>% collect %>% .[["UniqueID"]] %>%	as.character()
	}	else if (ProjectInfo$Species=="rat"){
		terminals_id <- dplyr::filter(ProteinGeneName, Gene.Name %in% terminalsdf.set$rat_Gene)  %>%
		dplyr::select(UniqueID) %>% collect %>% .[["UniqueID"]] %>%	as.character()
	}

	subdatlong <- dplyr::filter(data_long, UniqueID %in% terminals_id ) %>%
	group_by(., group, UniqueID) %>%
	dplyr::summarise(mean=mean(expr, na.rm = TRUE))
	subdatlong<-subdatlong%>%left_join(ProteinGeneName)
	subdatwide <- subdatlong  %>%
	tidyr::spread(.,group, mean, fill = 0) %>%
	as.data.frame() %>%
	remove_rownames(.) %>%
	column_to_rownames(.,var="UniqueID") %>%
	dplyr::select(one_of(as.character(group_order())))
	subdatwide=data.matrix(subdatwide)
	if (input$gs_heatmap_label=="Gene.Name") {
		sel_col=match(rownames(subdatwide), ProteinGeneName$UniqueID)
		rownames(subdatwide)=ProteinGeneName$Gene.Name[sel_col]
	}

	#remove rows with same values across all samples, which can cause hcluster error
	row_SD=apply(subdatwide, 1, function(x) sd(x,na.rm=T))
	subdatwide=subdatwide[row_SD!=0, ]
	#p <- pheatmap::pheatmap(as.matrix(subdatwide),	scale = "row", color = colorpanel (64, low = "blue",mid = "white", high = "red"),filename=NA)
	scaled_data=t(scale(t(subdatwide))); scaled_data=pmin(scaled_data, 3); scaled_data=pmax(scaled_data, -3)
	p <- ComplexHeatmap::Heatmap(scaled_data,
      column_names_gp = gpar(fontsize = as.numeric(as.character(input$hxfontsize_gsh))),
      row_names_gp = gpar(fontsize = as.numeric(as.character(input$hyfontsize_gsh))),
      column_title_gp = gpar(fontsize = as.numeric(as.character(input$htfontsize_gsh))),
      heatmap_legend_param = list(
        title = "Z Score",
        color_bar = "continuous",
        title_gp = gpar(fontsize = as.numeric(as.character(input$hlfontsize_gsh))),
        labels_gp = gpar(fontsize = as.numeric(as.character(input$hlfontsize_gsh))-1)
      ),
      column_title = ID, cluster_columns =F)
	return(p)
})
})

output$SetHeatMap = renderPlot({
	#grid.draw(genesetheatmap_out()$gtable)
	draw(genesetheatmap_out(), merge_legend=T,  auto_adjust = FALSE)
}, height=800) 


observeEvent(input$genesetheatmap, {
	#ID = input$x3
 # saved_plots$genesetheatmap[[ID]] <- genesetheatmap_out() #this only works on R4.0
  saved_plots$genesetheatmap<- genesetheatmap_out() #this works on R3.5 - 3.6
  
  #saved_plots$genesetheatmap <- genesetheatmap_out()$gtable
}
)

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

	DT::datatable(terminalsdf.set,  extensions = 'Buttons', options = list( dom = 'lBfrtip', buttons = c('csv', 'excel', 'print'), pageLength = 15))
})
