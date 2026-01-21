###########################################################################################################
##This software belongs to Biogen Inc. All right reserved.
##
##@file: correlationfmodule.R
##@Developer : Wei Li (wei.li@Biogen.com)
##@Date : 11/06/2025
##@version 1.0
###########################################################################################################
## correlation analysis
##########################################################################################################
#pkgs:  

library(shiny)
library(DT)
library(ggplot2)

get_regression_stats <- function(x, y, corr_method = "pearson") {
  if (corr_method == "spearman") {
    x_fit <- rank(x, ties.method = "average")
    y_fit <- rank(y, ties.method = "average")
  } else {
    x_fit <- x
    y_fit <- y
  }
  fit <- lm(x_fit ~ y_fit)
  s   <- summary(fit)
  res_corr <- cor.test(x, y, method = corr_method, exact = FALSE)
  data.frame(
    intercept   = coef(fit)[1],
    slope       = coef(fit)[2],
    r_squared   = s$r.squared,
    correlation = unname(res_corr$estimate),
    p.value     = res_corr$p.value
  )
}

get_gene_counts <- function(ProteinGeneName, gene_list, UniqueIDonly = F) {
  ProteinGeneName_sel <- data.frame()
  gene_label_reset <- F
  ProteinGeneName_sel <- dplyr::filter(ProteinGeneName, (UniqueID %in% gene_list) | (Protein.ID %in% gene_list) | (toupper(Gene.Name) %in% toupper(gene_list)))
  if (UniqueIDonly) {
    ProteinGeneName_sel <- dplyr::filter(ProteinGeneName, UniqueID %in% gene_list)
  }

  if (nrow(ProteinGeneName_sel) < 2) {
    msg_filter1 <- "Please input at least 2 matched gene."
    msg_filter <- ""
  } else {
    msg_filter1 <- paste("Selected Genes:",length(ProteinGeneName_sel %>% dplyr::pull(Gene.Name) %>% unique()),sep="")
    msg_filter <- paste("Selected IDs:", length(ProteinGeneName_sel %>% dplyr::pull(UniqueID) %>% unique()),sep="")
    if (length(ProteinGeneName_sel %>% dplyr::pull(UniqueID) %>% unique()) > length(ProteinGeneName_sel %>% dplyr::pull(Gene.Name) %>% unique())) {
      msg_filter <- paste(msg_filter, " (Gene(s) match multiple IDs, use UniqueID for Gene Label)", sep="")
      gene_label_reset <- T
    }
  }
  return(list("msg_filter1" = msg_filter1, "msg_filter" = msg_filter, "ProteinGeneName_sel" = ProteinGeneName_sel, "gene_label_reset" = gene_label_reset))
}

# ---- UI Module ----
correlation_ui <- function(id) {
  ns <- shiny::NS(id)
  fluidRow(
    column(3,
           wellPanel(
             column(width=12,uiOutput(ns("selectGroupSample"))),
             radioButtons(ns("correlation_type"), "Correlation Type:", choices = c("Gene-Gene" = "gene", "Sample-Sample" = "sample", "Group-Group" = "group"), inline = TRUE, selected = "gene"),
             radioButtons(ns("gene_subset"),label="Genes Used in Correlation Analysis", choices=c(""),inline = TRUE),
             conditionalPanel(ns=ns, "input.gene_subset=='Select'",
                              selectizeInput(ns("sel_gene"),	label="Gene Name (Select 2 or more)",	choices = NULL,	multiple=TRUE, options = list(placeholder =	'Type to search')),
                              span(textOutput(ns("filteredgene_Select")), style = "color:red; font-size:15px; font-family:arial; font-style:italic"),
                              span(textOutput(ns("filteredUniqueID_Select")), style = "color:red; font-size:15px; font-family:arial; font-style:italic"),
                              tags$hr(style="border-color: black;")
             ),
             conditionalPanel(ns=ns, "input.gene_subset=='Upload Genes'",
                              textAreaInput(ns("gene_list"), label="Enter Gene List", "", cols = 5, rows=6),
                              span(textOutput(ns("filteredgene_Upload")), style = "color:red; font-size:15px; font-family:arial; font-style:italic"),
                              span(textOutput(ns("filteredUniqueID_Upload")), style = "color:red; font-size:15px; font-family:arial; font-style:italic"),
                              tags$hr(style="border-color: black;")
             ),
             conditionalPanel(ns=ns, "input.gene_subset=='Geneset'",
                              selectizeInput(ns("sel_geneset"), label="Available GeneSet", choices = NULL, multiple = FALSE),
                              textAreaInput(ns("geneset_genes"), "Genes in Geneset", "", cols = 5, rows=6),
                              span(textOutput(ns("filteredgene_Geneset")), style = "color:red; font-size:15px; font-family:arial; font-style:italic"),
                              span(textOutput(ns("filteredUniqueID_Geneset")), style = "color:red; font-size:15px; font-family:arial; font-style:italic"),
                              tags$hr(style="border-color: black;")
             ),
             conditionalPanel(ns=ns, "input.gene_subset=='Browsing'",
                              selectInput(ns("sel_test"), label="Select Test", choices=NULL),
                              fluidRow(
                                column(width=6, radioButtons(ns("psel"), label= "P value or P.adj Value?", choices= c("Pval"="Pval","Padj"="Padj"), inline = TRUE)),
                                column(width=6, radioButtons(ns("updown"), label= "All, Up or Down?", choices= c("All"="All","Up"="Up","Down"="Down"), inline = TRUE))
                              ),
                              fluidRow(
                                column(width=6, numericInput(ns("fccut"), label= "Fold Change Threshold", value = 1.2, min=1, step=0.1)),
                                column(width=6, numericInput(ns("pvalcut"), label= "P-value Threshold", value=0.01, min=0, step=0.01))
                              ),
                              span(textOutput(ns("filteredgene_Browsing")), style = "color:red; font-size:15px; font-family:arial; font-style:italic"),
                              span(textOutput(ns("filteredUniqueID_Browsing")), style = "color:red; font-size:15px; font-family:arial; font-style:italic"),
                              tags$hr(style="border-color: black;")
             ),
             conditionalPanel(ns = ns, "input.correlation_type=='gene'",
                              radioButtons(ns("gene_label"),label="Select Gene Label",inline = TRUE, choices=c("UniqueID", "Gene.Name"), selected="Gene.Name")
             ),
             conditionalPanel(ns = ns, "input.correlation_type=='sample'",
                              selectizeInput(ns("sel_sample"),	label="Select Samples",	choices = NULL,	multiple=TRUE),
                              span(textOutput(ns("Selected_samples")), style = "color:red; font-size:15px; font-family:arial; font-style:italic"),
                              tags$hr(style="border-color: black;")
             ),             
             conditionalPanel(ns = ns, "input.correlation_type=='group'",
                              selectizeInput(ns("sel_attribute"), label="Select a Attribute", choices = NULL, multiple = FALSE),
                              selectizeInput(ns("sel_group"),	label="Select Groups",	choices = NULL,	multiple=TRUE),
                              span(textOutput(ns("Selected_groups")), style = "color:red; font-size:15px; font-family:arial; font-style:italic"),
                              tags$hr(style="border-color: black;")
             ),
             radioButtons(ns("correlation_method"), "Correlation Method:", choices = c("Pearson", "Spearman"), inline = TRUE, selected = "Pearson"),
           )
    ),
    column(9,
           tabsetPanel(id = ns("correlation_tabset"),
                       tabPanel(title = "Result Table", 
                                actionButton(ns("compute_corr"), "Compute/Refresh",style = "color: #0961E3; background-color: #F6E98C; border-color: #2e6da4"),
                                tags$br(), 
                                tags$hr(),
                                DT::dataTableOutput(ns("corr_table"))),
                       tabPanel(title = "Correlation Plot", 
                                actionButton(ns("CorrPlot"), "Save to output"),
                                tags$br(), 
                                tags$hr(),
                                plotOutput(ns("corr_plot"), width = '800px', height = "800px")),
                       tabPanel(title = "Data Table", 
                                DT::dataTableOutput(ns("data_table")))
           )         
    )
  )}

# ---- Server Module ----
correlation_server <- function(id) {
  moduleServer(id, 
               function(input, output, session) {
                 ProteinGeneName_sel <- reactiveVal()
                 CorrPlot <- reactiveVal()
                 CorrPlot_ID <- reactiveVal()
                 genelabel <- reactiveVal()
                 CorrMethod <- reactiveVal()
                 observe({
                   DataIn = DataReactive()
                   ProteinGeneName = DataIn$ProteinGeneName
                   MetaData = DataIn$MetaData
                   
                   if (input$gene_label=="UniqueID") {
                     DataIngenes <- ProteinGeneName %>% dplyr::select(UniqueID) %>% collect %>% .[["UniqueID"]] %>%	as.character()
                   } else {
                     DataIngenes <- ProteinGeneName %>% dplyr::select(Gene.Name) %>% collect %>% .[["Gene.Name"]] %>%	as.character()
                   }
                   updateSelectizeInput(session,'sel_gene',choices=DataIngenes, server = TRUE)
                 })
                 
                 observe({
                   tests = all_tests()
                   ProteinGeneName_Header = ProteinGeneNameHeader()
                   updateSelectizeInput(session,'sel_test',choices=tests, selected=tests[1])
                 })
                 
                 observe({
                   preset_samples=sample_order()
                   updateSelectizeInput(session,'sel_sample',choices=preset_samples)
                 })
                 
                 observe({
                   DataIn = DataReactive()
                   MetaData = DataIn$MetaData
                   preset_group=group_order() 
                   attributes=sort(setdiff(colnames(MetaData), c("sampleid", "Order", "ComparePairs") ))
                   updateSelectInput(session, "sel_attribute", choices=attributes, selected="group")
                   updateSelectizeInput(session,'sel_group',choices=preset_group)
                 })
                 
                 observe({
                   CorrMethod(input$correlation_method)
                 })
                 
                 observeEvent(input$correlation_type, {
                   if (input$correlation_type == 'gene') {
                     updateRadioButtons(
                       session,
                       "gene_subset",
                       choices = c("Select", "Browsing", "Upload Genes", "Geneset"),
                       inline = TRUE,  
                       selected = "Select"
                     )
                   } else {
                     updateRadioButtons(
                       session,
                       "gene_subset",
                       choices = c("All", "Browsing", "Upload Genes", "Geneset"),
                       inline = TRUE,
                       selected = "All"
                     )
                   }
                 })
                 
                 observeEvent(input$sel_gene, {
                   req(input$gene_subset == "Select")
                   req(input$sel_gene!="")
                   DataIn = DataReactive()
                   ProteinGeneName = DataIn$ProteinGeneName
                   gene_list = input$sel_gene
                   res_gene_count <- get_gene_counts(ProteinGeneName, gene_list)
                   ProteinGeneName_sel(res_gene_count$ProteinGeneName_sel)
                   genelabel(input$gene_label)
                   if (res_gene_count$gene_label_reset) genelabel("UniqueID")
                   output$filteredgene_Select <- renderText({res_gene_count$msg_filter1})
                   output$filteredUniqueID_Select <- renderText({res_gene_count$msg_filter})
                 })
                 
                 observe({
                   p_sel   <- input$psel
                   test_sel <- input$sel_test
                   FCcut <- log2(as.numeric(input$fccut))
                   pvalcut <- as.numeric(input$pvalcut)
                   Updown <- input$updown
                   
                   DataIn = DataReactive()
                   results_long = DataIn$results_long
                   ProteinGeneName = DataIn$ProteinGeneName

                   tmpdat <- GeneFilter(results_long, test_sel, p_sel, Updown, pvalcut, FCcut,'UniqueID')
                   tmpids <- tmpdat %>% as.data.frame() %>% dplyr::pull(UniqueID)
                   
                   res_gene_count <- get_gene_counts(ProteinGeneName, tmpids, UniqueIDonly = T)
                   ProteinGeneName_sel(res_gene_count$ProteinGeneName_sel)
                   genelabel(input$gene_label)
                   if (res_gene_count$gene_label_reset) genelabel("UniqueID")
                   output$filteredgene_Browsing <- renderText({res_gene_count$msg_filter1})
                   output$filteredUniqueID_Browsing <- renderText({res_gene_count$msg_filter})
                 })
                 
                 observeEvent(input$gene_list, {
                   req(input$gene_subset == "Upload Genes")
                   req(input$gene_list)
                   gene_list <- input$gene_list
                   gene_list <- ProcessUploadGeneList(gene_list)

                   DataIn = DataReactive()
                   ProteinGeneName = DataIn$ProteinGeneName

                   res_gene_count <- get_gene_counts(ProteinGeneName, gene_list)
                   ProteinGeneName_sel(res_gene_count$ProteinGeneName_sel)
                   genelabel(input$gene_label)
                   if (res_gene_count$gene_label_reset) genelabel("UniqueID")
                   output$filteredgene_Upload <- renderText({res_gene_count$msg_filter1})
                   output$filteredUniqueID_Upload <- renderText({res_gene_count$msg_filter})
                 })
                 
                 observeEvent(input$gene_subset , {
                   req(input$gene_subset == "Geneset")
                   genesetnames <- GetGeneSetNames()
                   updateSelectizeInput(session, "sel_geneset", choices =  c('Type to Search' = '', genesetnames), server = TRUE)
                 })
                 
                 observeEvent(input$sel_geneset, {
                   req(input$gene_subset == "Geneset")
                   req(input$sel_geneset!="")
                   sel_geneset <- input$sel_geneset
                   gene_list <- GetGenesFromGeneSet(sel_geneset)
                   updateTextAreaInput(session, "geneset_genes", value=paste(gene_list, collapse=","))
                   
                   DataIn = DataReactive()
                   ProteinGeneName = DataIn$ProteinGeneName

                   res_gene_count <- get_gene_counts(ProteinGeneName, gene_list)
                   genelabel(input$gene_label)
                   if (res_gene_count$gene_label_reset) genelabel("UniqueID")
                   ProteinGeneName_sel(res_gene_count$ProteinGeneName_sel)
                   output$filteredgene_Geneset <- renderText({res_gene_count$msg_filter1})
                   output$filteredUniqueID_Geneset <- renderText({res_gene_count$msg_filter})
                 })
                 
                 observeEvent(input$sel_sample, {
                   sample_list = input$sel_sample
                   
                   msg_sample <- paste("Selected Samples:",length(sample_list),sep="")
                   output$Selected_samples <- renderText({msg_sample})
                 })
                 
                 observeEvent(input$sel_attribute, {
                   req(input$sel_attribute)
                   DataIn = DataReactive()
                   MetaData = DataIn$MetaData %>% filter(sampleid %in% sample_order())
                   groups = unique(MetaData[, input$sel_attribute])
                   updateSelectizeInput(session,'sel_group',choices=groups)
                   output$Selected_groups <- renderText({""})
                 })
                   
                 observeEvent(input$sel_group, {
                   group_list = input$sel_group
                 
                   msg_group <- paste("Selected Groups:",length(group_list),sep="")
                   output$Selected_groups <- renderText({msg_group})
                 })

                 CorrResult <- eventReactive(input$compute_corr, { 
                   DataIn <- DataReactive()
                   data_long <- DataIn$data_long
                   preset_group=group_order() 
                   preset_samples=sample_order()
                   if (input$gene_subset == "All") {
                     tmpids <- DataIn$ProteinGeneName %>% dplyr::pull(UniqueID) %>% unique()
                   } else {
                     ProteinGeneName_sel <- ProteinGeneName_sel() 
                     genelabel <- genelabel() 
                     validate(need(nrow(ProteinGeneName_sel) > 1, message = "Please input at least 2 matched gene."))
                     tmpids <- ProteinGeneName_sel %>% dplyr::pull(UniqueID) %>% unique()
                   }
                   if (input$correlation_type=='gene') {
                     exp_tmp = data_long %>% 
                       dplyr::filter(UniqueID %in% tmpids, group %in% preset_group, sampleid %in% preset_samples) %>%
                       dplyr::select(gene = !!sym(genelabel), sampleid, expr) %>%
                       dplyr::filter(!is.na(expr)) %>% 
                       tidyr::pivot_wider(
                         id_cols = gene,          
                         names_from = sampleid,       
                         values_from = expr           
                       ) %>% 
                       tibble::column_to_rownames("gene") %>%
                       as.matrix()
                   } else if (input$correlation_type=='sample') {
                     exp_tmp = data_long %>% 
                       dplyr::filter(UniqueID %in% tmpids, group %in% preset_group, sampleid %in% input$sel_sample) %>%
                       dplyr::select(gene = UniqueID, sampleid, expr) %>%
                       dplyr::filter(!is.na(expr)) %>% 
                       tidyr::pivot_wider(
                         id_cols = sampleid,          
                         names_from = gene,       
                         values_from = expr           
                       ) %>% 
                       tibble::column_to_rownames("sampleid") %>%
                       as.matrix()
                   } else if (input$correlation_type=='group') {
                     adding_number <- ifelse(exp_unit() == "Expression Level", 0, as.numeric(exp_unit()))
                     data_long <- DataIn$data_long
                     sel_attribute <- input$sel_attribute
                     sel_group <- input$sel_group
                     MetaData = DataIn$MetaData %>% filter(sampleid %in% sample_order(), !!sym(sel_attribute) %in% sel_group)
                     
                     exp_tmp = data_long %>% 
                       dplyr::left_join(MetaData %>% dplyr::select(sampleid, !!sym(sel_attribute)),
                                 by = "sampleid") %>%
                       dplyr::filter(UniqueID %in% tmpids, sampleid %in% MetaData$sampleid) %>%
                       dplyr::select(gene = UniqueID, !!sym(sel_attribute), expr) %>%
                       dplyr::filter(!is.na(expr)) %>% 
                       dplyr::mutate(TPM = 2^expr-adding_number) %>%
                       dplyr::group_by(!!sym(sel_attribute), gene) %>%
                       dplyr::summarise(mean_TPM = mean(TPM), .groups = "drop") %>%
                       dplyr::mutate(log2_mean_TPM = log2(mean_TPM + adding_number)) %>%
                       tidyr::pivot_wider(
                         id_cols = !!sym(sel_attribute),          
                         names_from = gene,       
                         values_from = log2_mean_TPM           
                       ) %>% 
                       tibble::column_to_rownames(sel_attribute) %>%
                       as.matrix()
                   }
                   t_exp_tmp <- t(exp_tmp)   # genes as columns
                   toCheck_list <- colnames(t_exp_tmp)
                   
                   corr_method = tolower(CorrMethod())

                   res_corr <- do.call(rbind, combn(toCheck_list, 2, FUN = function(pair) {
                     stats <- get_regression_stats(t_exp_tmp[, pair[1]], t_exp_tmp[, pair[2]], corr_method)
                     cbind(gene1 = pair[1], gene2 = pair[2], stats)
                   }, simplify = FALSE))
                   
                   if (input$correlation_type=='sample') {
                     names(res_corr)[1:2] = c('sample1', 'sample2')
                   } else if (input$correlation_type=='group') {
                     names(res_corr)[1:2] = c('group1', 'group2')
                   }
                   
                   res_corr_sorted <- res_corr[order(-res_corr$correlation), ]
                   rownames(res_corr_sorted) <- seq(1:nrow(res_corr_sorted))
                   res_corr_sorted$rank <- rownames(res_corr_sorted)
                   return(list("corr_table" = res_corr_sorted, "exp_filtered" = exp_tmp ))
                 })
                 
                 output$corr_table <-  DT::renderDT({
                   res <- CorrResult()
                   corr_table <- res$corr_table
                   
                   DT::datatable(
                     corr_table,  extensions = 'Buttons', escape = FALSE, selection = 'none', class = 'cell-border strip hover',
                     options = list(    dom = 'lBfrtip', pageLength = 15,
                                        buttons = list(
                                          list(extend = "csv", text = "Download Page", filename = "Page_results",
                                               exportOptions = list(modifier = list(page = "current"))),
                                          list(extend = "csv", text = "Download All", filename = "All_Results",
                                               exportOptions = list(modifier = list(page = "all")))
                                        )
                     )) %>% 
                     formatSignif(columns=c('intercept', 'slope', 'r_squared', 'correlation', 'p.value'), digits=3) %>% 
                     formatStyle(8, cursor = 'pointer',color='blue')
                 })
                 
                 observeEvent(input$corr_table_cell_clicked, {
                   clicked <- input$corr_table_cell_clicked
                   res <- CorrResult()
                   corr_table <- res$corr_table
                   exp_filtered <- res$exp_filtered
                   corr_method = CorrMethod()
                   
                   if (length(clicked) > 0) {
                     row <- clicked$row
                     gene1 <- corr_table[row, 1]
                     gene2 <- corr_table[row, 2]
                     
                     eqn <- sprintf("Correlation: %.3f\ny = %.3f %+.3fx\nR² = %.3f, p.value = %.3f",
                                    corr_table$correlation[row],
                                    corr_table$intercept[row],
                                    corr_table$slope[row],
                                    corr_table$r_squared[row],
                                    corr_table$p.value[row])
                     
                     if (corr_method == "Spearman") {
                       df_plot <- data.frame(
                         x = rank(exp_filtered[gene1, ], ties.method = "average"),
                         y = rank(exp_filtered[gene2, ], ties.method = "average")
                       )
                     } else {
                       df_plot <- data.frame(
                         x = exp_filtered[gene1, ],
                         y = exp_filtered[gene2, ]
                       )
                     }
                     
                     p <- ggplot(df_plot, aes(x, y)) +
                       geom_point() +
                       geom_smooth(method = "lm") +
                       labs(
                         x = gene1,
                         y = gene2,
                         title = paste(corr_method, " Correlation: ", paste(gsub(" ", "_", gene1), gsub(" ", "_", gene2), sep = "-"))
                       ) +
                       annotate("text",
                                x = min(df_plot$x, na.rm = TRUE),
                                y = max(df_plot$y, na.rm = TRUE),
                                label = eqn,
                                # hjust = 0, vjust = 1, size = 5, color = "blue") +
                                hjust = 0, vjust = 0.5, size = 5, color = "blue") +
                       theme(
                         axis.title.x = element_text(size = 16),   
                         axis.title.y = element_text(size = 16),   
                         axis.text.x  = element_text(size = 14),   
                         axis.text.y  = element_text(size = 14),   
                         plot.title   = element_text(size = 18, face = "bold") 
                       )
                     CorrPlot(p)
                     CorrPlot_ID(paste(gsub(" ", "_", gene1), gsub(" ", "_", gene2), sep = "-"))
                     output$corr_plot <- renderPlot({CorrPlot()})
                   }
                 })
                 
                 observeEvent(input$CorrPlot, {
                   saved_plots$CorrPlot[[CorrPlot_ID()]] <- CorrPlot()
                 })
                 
                 observeEvent(input$corr_table_cell_clicked, {
                   if (length(input$corr_table_cell_clicked) > 0) {
                     updateTabsetPanel(session, "correlation_tabset", selected = "Correlation Plot")
                   }
                 })
                 
                 output$data_table <- DT::renderDataTable({
                   req(public_dataset)
                   res <- CorrResult()
                   if (input$correlation_type == 'gene') {
                     exp_filtered <- res$exp_filtered
                   } else {
                     exp_filtered <- t(res$exp_filtered)
                   }
                   DT::datatable(exp_filtered, extensions = c('FixedColumns', 'Buttons'),
                                 options = list(
                                   pageLength = 15,
                                   dom = 'lBfrtip', 
                                   buttons = list(
                                     list(extend = "csv", text = "Download Page", filename = "Page_results",
                                          exportOptions = list(modifier = list(page = "current"))),
                                     list(extend = "csv", text = "Download All", filename = "All_Results",
                                          exportOptions = list(modifier = list(page = "all")))
                                   ),
                                   scrollX = TRUE,
                                   fixedColumns = list(leftColumns = 1)
                                 )) %>% 
                     formatSignif(columns=1:ncol(exp_filtered), digits=3)
                 })
                 
                 observe({
                   if (public_dataset) {
                     showTab(inputId = "correlation_tabset", target = "data_table")
                   } else {
                     hideTab(inputId = "correlation_tabset", target = "data_table")
                   }
                 })
                 
               })
}



