###########################################################################################################
##This software belongs to Biogen Inc. All right reserved.
##
##@file: TimeSeriesmodule.R
##@Developer : Wei Li (wei.li@Biogen.com)
##@Date : 11/06/2025
##@version 1.0
###########################################################################################################
## correlation analysis
##########################################################################################################
#pkgs:  

library(DESeq2)
library(ggplot2)
library(DEGreport)

get_count_mtx <- function(data_long, linear_base, linear_small_value) {
  if (is.na(linear_base)) linear_base = 2
  if (is.na(linear_small_value)) linear_small_value = 1
  data_long_linear <- data_long %>% mutate(expr=linear_base^(expr-linear_small_value))
  
  data_wide_linear <- data_long_linear %>%
    dplyr::select(UniqueID, sampleid, expr) %>%
    tidyr::pivot_wider(
      names_from = sampleid,
      values_from = expr
    )
  data_wide_linear <- as.data.frame(data_wide_linear)
  rownames(data_wide_linear) <- data_wide_linear$UniqueID
  data_wide_linear$UniqueID <- NULL
  data_wide_linear[is.na(data_wide_linear)] <- 0
  
  data_wide_linear[] <- lapply(data_wide_linear, function(col) {
    if (is.numeric(col)) {
      as.integer(round(col))
    } else {
      col
    }
  })
  return(data_wide_linear)
}

parse_time <- function(x) {
  # Extract numeric part
  value <- as.numeric(sub("([0-9.]+).*", "\\1", x))
  
  # Extract unit part
  unit <- tolower(sub("[0-9.]+\\s*(.*)", "\\1", x))
  unit <- trimws(unit)
  
  # Unit conversion map (to hours)
  unit_map <- c(
    # seconds
    "s" = 1/3600, "sec" = 1/3600, "second" = 1/3600, "seconds" = 1/3600,
    # minutes
    "min" = 1/60, "mins" = 1/60, "minute" = 1/60, "minutes" = 1/60,
    # hours
    "h" = 1, "hr" = 1, "hrs" = 1, "hour" = 1, "hours" = 1,
    # days
    "d" = 24, "day" = 24, "days" = 24,
    # weeks
    "wk" = 24*7, "week" = 24*7, "weeks" = 24*7,
    # months
    "m" = 24*30, "mo" = 24*30, "mon" = 24*30, "month" = 24*30, "months" = 24*30,
    # years
    "yr" = 24*365, "year" = 24*365, "years" = 24*365
  )
  
  # Match units
  multiplier <- unit_map[unit]
  if (any(is.na(multiplier))) {
    stop("Unknown time units: ", paste(unique(unit[is.na(multiplier)]), collapse = ", "))
  }
  
  # Convert to hours
  hours <- value * multiplier
  
  # Canonical unit mapping
  canonical_unit <- c(
    "s"="s","sec"="s","second"="s","seconds"="s",
    "m"="min","min"="min","mins"="min","minute"="min","minutes"="min",
    "h"="h","hr"="h","hrs"="h","hour"="h","hours"="h",
    "d"="d","day"="d","days"="d",
    "wk"="wk","week"="wk","weeks"="wk",
    "mo"="mo","mon"="mo","month"="mo","months"="mo",
    "yr"="yr","year"="yr","years"="yr"
  )
  
  # Determine majority unit
  majority_unit <- names(sort(table(canonical_unit[unit]), decreasing = TRUE))[1]
  
  # Conversion factors (hours → target unit)
  unit_to_hours <- c(
    "s" = 1/3600,
    "min" = 1/60,
    "h" = 1,
    "d" = 24,
    "wk" = 24*7,
    "mo" = 24*30,
    "yr" = 24*365
  )
  
  # Convert hours → majority unit
  converted <- hours / unit_to_hours[majority_unit]
  
  # ---- NEW: format to max 3 decimals, no trailing zeros ----
  converted_fmt <- sub("\\.?0+$", "", format(round(converted, 3), nsmall = 3))
  
  # Build canonical labels
  canonical <- paste0(converted_fmt, majority_unit)
  
  # Order levels by numeric time
  ord <- order(hours)
  
  # Return UNORDERED factor with ordered levels
  factor(canonical, levels = unique(canonical[ord]), ordered = FALSE)
}

run_timecourse_DEG <- function(MetaData, count_mtx,time_variable,Condition) {
  library(DESeq2)
  library(dplyr)
  library(tibble)
  
  # -----------------------------
  # 1. Prepare metadata
  # -----------------------------
  sample_meta <- MetaData %>%
    filter(sampleid %in% colnames(count_mtx)) %>%
    arrange(match(sampleid, colnames(count_mtx)))
  # Convert all factor columns to character (DEGreport safety)
  sample_meta[] <- lapply(sample_meta, function(x) {
    if (is.factor(x)) as.character(x) else x
  })
  # # Create required columns
  # sample_meta$sample <- sample_meta$sampleid
  # 
  # Create ordered time factor
  if (!is.numeric(sample_meta[[time_variable]])) {
    sample_meta$nominaltimef <- parse_time(sample_meta[[time_variable]])
  } else {
    sample_meta$nominaltimef <- factor(sample_meta[[time_variable]], levels = sort(unique(sample_meta[[time_variable]])))
    
  }
  # -----------------------------
  # 1b. Factor all Condition columns
  # -----------------------------
  for (col in Condition) {
    sample_meta[[col]] <- factor(sample_meta[[col]])
  }
  
  # Last covariate is the one to test interaction
  interaction_cov <- tail(Condition, 1)
  
  # Set rownames
  rownames(sample_meta) <- sample_meta$sampleid
  
  # -----------------------------
  # 2. Build DESeq2 model formulas
  # -----------------------------
  # Main effects: all covariates
  main_effects <- paste(Condition, collapse = " + ")
  
  # Full model: all covariates + time + lastCov:time
  full_formula <- as.formula(
    paste("~", main_effects, "+ nominaltimef +", 
          paste0(interaction_cov, ":nominaltimef"))
  )
  
  # Reduced model: drop interaction
  reduced_formula <- as.formula(
    paste("~", main_effects, "+ nominaltimef")
  )
  
  # -----------------------------
  # 3. Build DESeq2 object
  # -----------------------------
  dds <- tryCatch({
    DESeqDataSetFromMatrix(
      countData = count_mtx,
      colData   = sample_meta,
      design    = full_formula
    )
  }, error = function(e) {
    if (grepl("full rank", e$message)) {
      return("RANK_ERROR")
    }
    stop(e$message)
  })
  
  # Use validate to show the message in the UI
  validate(
    need(!is.character(dds) || dds != "RANK_ERROR", "Error: The model is not full rank. Check your metadata for confounded variables.")
  )
  
  # Filter low-count genes
  dds <- dds[rowSums(counts(dds)) > 1, ]
  # LRT test
  dds <- DESeq(dds, test = "LRT", reduced = reduced_formula)
  
  # # -----------------------------
  # # 4. Extract LRT results
  # # -----------------------------
  # res_lrt <- results(dds)
  return(dds)
}

run_DEG_cluster <- function(dds, sig_gene, Condition) {
  res_lrt <- results(dds)
  sample_meta <- as.data.frame(colData(dds))
  
  # clustering_sig_genes <- res_lrt %>%
  #   data.frame() %>%
  #   rownames_to_column("gene") %>%
  #   as_tibble() %>%
  #   filter(padj < padj_cutoff) %>%
  #   head(n = top_n_cutoff)
  
  rld <- vst(dds, blind = TRUE)
  rld_mat <- assay(rld)
  cluster_rlog <- rld_mat[sig_gene, ]
  
  clusters <- degPatterns(
    cluster_rlog,
    metadata = sample_meta,
    time     = "nominaltimef",
    col      = tail(Condition, 1),
    plot     = FALSE
  )
  
  return(list(
    rlog = rld_mat,
    cluster_rlog = cluster_rlog,
    clusters = clusters
  ))
}

# ---- UI Module ----
TimeSeries_ui <- function(id) {
  ns <- shiny::NS(id)
  fluidRow(
    column(3,
           wellPanel(
             column(width=12,uiOutput(ns("selectGroupSampleTimeSeries"))), 
             radioButtons(ns("ts_genelable"),label="Select Gene Label",inline = TRUE, choices=c("Gene.Name","UniqueID"), selected="Gene.Name"),
             conditionalPanel(ns = ns, "input.TimeSeries_tabset =='Time-Series DE Analysis'",
                              selectInput(ns("sel_time_var"), label="Select Time Variable", choices=NULL),
                              selectizeInput(ns("sel_condition_var"), label="Select Condition Variable", choices=NULL,multiple=TRUE),
                              actionButton(ns("compute_DE"), "Run Time-Series DE Analysis", style="color: #0961E3; background-color: #F6E98C ; border-color: #2e6da4"),
                              tags$hr()
             ),
             conditionalPanel(ns = ns, "input.TimeSeries_tabset =='Time-Series DE Analysis' || input.TimeSeries_tabset =='DEG Clustering Analysis'",
                              sliderInput(ns("padjCutoff"), label="Set the Adjusted p-value Cutoff:", min = 0, max = 1, step = 0.01, value = 0.05),
                              sliderInput(ns("top_n"), label="Select Number of Top DEGs for Clustering", min = 300, max = 1500, step = 100, value = 1000),
                              span(textOutput(ns("filtered_DEG")), style = "color:blue; font-size:15px; font-family:arial; font-style:italic")
             ),
             conditionalPanel(ns = ns, "input.TimeSeries_tabset =='DEG Clustering Analysis'",
                              actionButton(ns("compute_cluster"), "Run Clustering", style="color: #0961E3; background-color: #F6E98C ; border-color: #2e6da4"),
                              conditionalPanel(ns = ns, "input.ts_cluster_tabset =='Time Series Plot'",
                                               tags$hr(),
                                               tags$p(
                                                 "Customize Plot:",
                                                 style = "font-size: 20px; font-weight: bold;"
                                               ),
                                               sliderInput(ns("plot_width"), label="Plot Width", min = 600, max = 2000, step = 100, value = 1200),
                                               sliderInput(ns("plot_height"), label="Plot Height", min = 600, max = 2000, step = 100, value = 800),
                                               sliderInput(ns("strip_text"), "Facet Title Size", min = 12, max = 25, value = 18),
                                               column(width=6,sliderInput(ns("axis_title"), "Axis Title Size", min = 12, max = 25, value = 20)),
                                               column(width=6,sliderInput(ns("axis_text"), "Axis Text Size", min = 12, max = 25, value = 18)),
                                               column(width=6,sliderInput(ns("legend_title"), "Legend Title Size", min = 12, max = 25, value = 20)),
                                               column(width=6,sliderInput(ns("legend_text"), "Legend Text Size", min = 12, max = 25, value = 18))
                              )
             )
           )
    ),
    column(9,
           tabsetPanel(id = ns("TimeSeries_tabset"),
                       tabPanel(title="Time-Series DE Analysis",
                                br(),
                                tags$details(
                                  tags$summary(
                                    tags$strong("Click to view DESeq2 LRT Methodology & Example", 
                                                style = "color: #007bff; cursor: pointer;")
                                  ),
                                  tags$div(
                                    style = "padding: 15px; background-color: #f0f7fb; border-left: 5px solid #007bff; margin-top: 10px;",
                                    tags$p("DESeq2 offers the Likelihood ratio test (LRT) test which is used to identify any genes that show change in expression across the different time points."),
                                    tags$p(
                                      "The LRT is comparing the full model to the reduced model to identify significant genes. ", 
                                      tags$strong("The p-values are determined solely by the difference in deviance between the ‘full’ and ‘reduced’ model formula (not log2 fold changes)."), 
                                      "Essentially the LRT test is testing whether the term(s) removed in the ‘reduced’ model explains a significant amount of variation in the data"
                                    ),
                                    tags$p(tags$strong("Example:")),
                                    tags$p(
                                      "User select ", 
                                      tags$strong("Age"), 
                                      " as the Time Variable and ",
                                      tags$strong("Tissue, Treatment"),     
                                      " as the Condition Variables."
                                    ),
                                    tags$p("Then the full model will be:"),
                                    tags$p(tags$strong("Tissue + Treatment + Treatment:Age")),
                                    tags$p(" and the reduced module will be:"),
                                    tags$p(tags$strong("Tissue + Treatment")),
                                    tags$p(
                                      "It always test the effect of the last variable in the Select Condition Variables input box with time. ",
                                      tags$strong("Any DEGs identified will show different time-course trends between the Condition groups.")
                                    )
                                  )
                                ),
                                br(),
                                actionButton(ns("save_ts_deg"), "Save to output"),
                                br(),br(),
                                DT::DTOutput(ns("ts_deg"))),
                       tabPanel(title="DEG Clustering Analysis",
                                tabsetPanel(id=ns("ts_cluster_tabset"),
                                            tabPanel(title = "Time Series Plot",
                                                     br(),
                                                     tags$details(
                                                       # This is the clickable header that is always visible
                                                       tags$summary(
                                                         tags$strong("Click to view Clustering Methodology (degPatterns & DIANA)", 
                                                                     style = "color: #007bff; cursor: pointer;")
                                                       ),
                                                       # Everything inside this div is hidden until expanded
                                                       tags$div(
                                                         style = "padding: 15px; background-color: #f8f9fa; border-radius: 5px; margin-top: 10px;",
                                                         tags$p(tags$strong("The degPatterns() R function in DEGreport library is used to do the DEG clustering.")),
                                                         tags$p(
                                                           "It can work with one or more groups with 2 or more several time points. Before calculating the genes similarity among samples, all samples inside the same time point (time parameter) and group (col parameter) are collapsed together, and the mean value is the representation of the group for the gene abundance. Then, all pair-wise gene expression is calculated using cor.test R function using kendall as the statistical method. A distance matrix is created from those values. After that, cluster::diana() is used for the clustering of gene-gene distance matrix and cut the tree using the divisive coefficient of the clustering, giving as well by diana."
                                                         ),
                                                         tags$p(tags$strong("How cluster::diana() Calculates Clusters")),
                                                         tags$p(
                                                           'The cluster::diana() function, which stands for DIvisive ANAlysis, is a "top-down" hierarchical clustering algorithm. Unlike common "bottom-up" (agglomerative) methods like hclust, it starts with all your data in one single cluster and systematically splits them. The algorithm follows a divisive process to build a hierarchy:'
                                                         ), 
                                                         tags$p(
                                                           tags$strong("Start at the Top:"), 
                                                           " Initially, all observations (genes, in your case) are in one large cluster."
                                                         ),
                                                         tags$p(
                                                           tags$strong('Find the "Splinter" Group:'), 
                                                           " In each step, the algorithm identifies the cluster with the largest diameter (the maximum distance between any two elements)."
                                                         ),
                                                         tags$p(
                                                           tags$strong("The Seed:"), 
                                                           'Within that cluster, it finds the most "disparate" element—the one with the highest average distance to all other elements in the cluster. This element starts a new "splinter group."'
                                                         ),
                                                         tags$p(
                                                           tags$strong("Reassignment:"), 
                                                           " The algorithm then looks at all other elements in the original cluster. If an element is closer to the splinter group than to the remainder of the original group, it is moved to the splinter group."
                                                         ),
                                                         tags$p(
                                                           tags$strong("Iterate:"), 
                                                           " This process repeats until every observation is its own individual cluster, creating a dendrogram (tree)."
                                                         ),
                                                         tags$p(
                                                           tags$strong("Distance Metric:"), 
                                                           " By default, it uses Euclidean distance, though the degPatterns function specifically calculates a distance matrix using Kendall correlation before passing it to diana()."
                                                         )
                                                       )
                                                     ),
                                                     br(),
                                                     actionButton(ns("ts_cluster_plot"), "Save to output"),
                                                     br(),br(),
                                                     plotOutput(ns("ts_plot"), width = '1200px', height = "1200px")),
                                            tabPanel(title="DEG Cluster Result",
                                                     br(),
                                                     DT::DTOutput(ns("ts_cluster_table")),
                                                     tags$script(HTML(sprintf("
                                                     $(document).on('click', '[id^=plotExp_]', function() {
                                                     var genes = $(this).attr('data-genes');
                                                     Shiny.setInputValue('%s', genes, {priority: 'event'});
                                                     });", ns("plotExp_trigger"))))
                                            ),
                                            tabPanel(title="Time‑series Heatmap",
                                                     br(),
                                                     plotOutput(ns("ts_heatmap"), height = 800))
                                )
                       )
           )         
    )
)}

# ---- Server Module ----
TimeSeries_server <- function(id, parent_session) {
  moduleServer(id, 
               function(input, output, session) {
                 linear_small_value <- reactiveVal()
                 linear_base <- reactiveVal()
                 
                 output$selectGroupSampleTimeSeries <- renderUI(shared_header_content())

                 observe( {
                   expU = exp_unit()
                   linear_base(as.numeric(str_replace(str_extract(expU,"log\\d+"),"log","")))
                   linear_small_value(as.numeric(str_replace(str_split_fixed(expU, "\\+", 2)[2], "\\)", "")))
                   unit=str_replace_all(str_extract(expU, "\\(.+\\+"), "(\\(|\\+)", "")
                 })
                 
                 observeEvent(DataQCReactive(), {
                   req(DataQCReactive())
                   DataIn = DataQCReactive()
                   MetaData=DataIn$MetaData
                   MetaData_clean <- MetaData[, sapply(MetaData, function(x) length(unique(x)) > 1)]
                   attributes=sort(setdiff(colnames(MetaData_clean), c("sampleid", "Order", "ComparePairs")))
                   updateSelectInput(session, "sel_time_var", choices=attributes, selected=NULL)  
                   updateSelectizeInput(session,'sel_condition_var',choices=attributes, selected=NULL)
                 })
                 
                 DEGReactive <- eventReactive(input$compute_DE, {
                   withProgress(message = 'Processing...', value = 0, {
                     req(DataQCReactive())
                     req(input$sel_time_var)
                     req(input$sel_condition_var)
                     DataIn = DataQCReactive()
                     data_long <- DataIn$tmp_data_long
                     MetaData <- DataIn$MetaData
                     ProteinGeneName <- DataIn$ProteinGeneName
                     
                     count_mtx <- get_count_mtx(data_long, linear_base(), linear_small_value())
                     
                     run_timecourse_DEG(MetaData,
                                        count_mtx,
                                        time_variable = input$sel_time_var,
                                        Condition = input$sel_condition_var
                     )
                   })
                 })
                 
                 # filtered_DE <- eventReactive(list(input$ts_genelable, input$padjCutoff, input$top_n), {
                 filtered_DE <- reactive({
                   req(DEGReactive())
                   dds <- DEGReactive()
                   
                   df<- results(dds) %>%
                     data.frame() %>%
                     rownames_to_column("gene") %>%
                     as_tibble() %>%
                     filter(padj < input$padjCutoff) %>%
                     arrange(padj) %>%
                     head(n = input$top_n)
                   
                   if (input$ts_genelable == 'Gene.Name') {
                     req(DataQCReactive())
                     ProteinGeneName <- DataQCReactive()$ProteinGeneName
                     df <- df %>% 
                       dplyr::left_join(ProteinGeneName %>% dplyr::select(UniqueID, Gene.Name), by = c("gene" = "UniqueID")) %>%
                       dplyr:: rename(UniqueID = gene) %>%
                       dplyr:: filter(!is.na(Gene.Name), Gene.Name != "") %>%
                       dplyr:: group_by(Gene.Name) %>%
                       dplyr:: slice_max(baseMean, n = 1, with_ties = FALSE) %>%
                       dplyr:: ungroup() %>%
                       dplyr:: rename(gene = Gene.Name) %>%
                       dplyr:: relocate(gene, .before = 1)
                   }
                   
                   nrow = nrow(df)
                   max_padj = max(df$padj)
                   list(
                     df = df,
                     msg = sprintf("Top %s DEGs are selected with maximum padj %.3f.", nrow, max_padj)
                   )
                 })
                     
                 output$filtered_DEG <- renderText({filtered_DE()$msg})
                     
                 observeEvent(input$save_ts_deg, {
                   df <- filtered_DE()$df
                   if (input$ts_genelable == 'Gene.Name') {
                     df <- df %>%
                       dplyr::select(-UniqueID)
                   }
                   num_cols <- names(df)[sapply(df, is.numeric)]
                   
                   saved_table$DEG_data <- df
                 })
                 
                 output$ts_deg<- DT::renderDT({
                   req(filtered_DE())
                   df <- filtered_DE()$df
                   if (input$ts_genelable == 'Gene.Name') {
                     df <- df %>%
                       dplyr::select(-UniqueID)
                   }
                   num_cols <- names(df)[sapply(df, is.numeric)]
                   DT::datatable(df, extensions = 'Buttons', options = list(
                     dom = "Blfrtip", buttons = c("csv", "excel", "print"), pageLength = 20), rownames= FALSE) %>% 
                     formatSignif(columns= num_cols, digits=3)
                 })
                 
                 DataClusterReactive <- eventReactive(input$compute_cluster, {
                   withProgress(message = 'Processing...', value = 0, {
                     req(DEGReactive())
                     req(filtered_DE())
                     req(input$sel_condition_var)
                     
                     dds <- DEGReactive()
                     filtered_dds <- filtered_DE()$df
                     sig_gene <- filtered_dds$gene
                     if (input$ts_genelable == 'Gene.Name') {
                       sig_gene <- filtered_dds$UniqueID
                     }
                     run_DEG_cluster(dds, sig_gene, input$sel_condition_var)
                   })
                 })
                 
                 cluster_table <- eventReactive(input$compute_cluster, {
                   withProgress(message = "Processing...", value = 0, {
                     req(DataClusterReactive())
                     req(DataQCReactive())
                     ProteinGeneName <- DataQCReactive()$ProteinGeneName
                     DataCluster <-DataClusterReactive()
                     clusters <- DataCluster$clusters
                     df = clusters$df

                     if (input$ts_genelable == 'Gene.Name') {
                       df = df %>% 
                         dplyr::left_join(ProteinGeneName %>% dplyr::select(UniqueID, Gene.Name), by = c("genes" = "UniqueID")) %>%
                         dplyr:: select(-genes) %>%
                         dplyr:: filter(!is.na(Gene.Name), Gene.Name != "") %>%
                         dplyr:: rename(genes = Gene.Name) %>%
                         dplyr:: relocate(genes, .before = 1)
                     }
                     df = df %>%
                       dplyr:: group_by(cluster) %>%
                       dplyr:: summarise(
                         n_gene = n(),
                         genes = paste(genes, collapse = ","), .groups = "drop"
                       ) %>%   
                       dplyr::mutate(
                         copy = purrr::map_chr(seq_len(n()), ~ paste0(
                           as.character(
                             rclipButton(
                               paste0("clipbtn_", .x),
                               label = "Copy all genes",
                               clipText = genes[.x],
                               icon = icon("copy", lib = "glyphicon"),
                               class = "btn-primary btn-sm"
                             )
                           ),
                           " ",
                           sprintf(
                             '<button id="plotExp_%s" class="btn btn-info btn-sm" data-genes="%s">Create Expression Plot</button>',
                             .x,
                             genes[.x]
                           )
                         ))
                       )  %>%
                       dplyr::select(cluster, n_gene, copy, genes)
                   })
                 })
                 
                 output$ts_cluster_table<- DT::renderDT({
                   req(cluster_table())
                   df <- cluster_table()
                   DT::datatable(df, 
                                 rownames = FALSE,
                                 escape = FALSE,
                                 selection = "none",
                                 colnames=c("Cluster", "Number of genes", "Action","Genes in cluster"),
                                 extensions = 'Buttons', 
                                 options = list(dom = "Blfrtip", 
                                                buttons = c("csv", "excel", "print")
                                                )
                                 ) 
                 })
                 
                 observeEvent(input$plotExp_trigger, {
                   gene_list <- input$plotExp_trigger
                   updateNavbarPage(parent_session, inputId = "menu", selected = "Exp_Plot")
                   updateTabsetPanel(parent_session, inputId = "expression_tabset", selected = "Searched Expression Data")
                   parent_session$onFlushed(function() {
                     updateSelectInput(parent_session, inputId = "exp_subset", selected = "Upload Genes")
                     updateTextAreaInput(parent_session, inputId = "exp_list", value = gene_list)
                   }, once = TRUE)
                 })
                 
                 cluster_plot <- eventReactive(input$compute_cluster, {
                   withProgress(message = "Processing...", value = 0, {
                     req(DataClusterReactive())
                     DataCluster <-DataClusterReactive()
                     clusters <- DataCluster$clusters
                     clusters$plot
                   })
                 })
                 
                 output$ts_plot<- renderPlot({
                   cluster_plot() +
                     theme(
                       text = element_text(size = input$base_text),        # all text
                       axis.text = element_text(size = input$axis_text),   # tick labels
                       axis.title = element_text(size = input$axis_title),  # axis titles
                       strip.text = element_text(size = input$strip_text),  # facet titles
                       legend.text = element_text(size = input$legend_text),
                       legend.title = element_text(size = input$legend_title)
                     )
                 },
                 width  = function() input$plot_width,
                 height = function() input$plot_height
                 )
                 
                 observeEvent(input$ts_cluster_plot, {
                   saved.num <- length(saved_plots$ts_cluster_plot) + 1
                   p <- cluster_plot() +
                     theme(
                       text = element_text(size = input$base_text),        # all text
                       axis.text = element_text(size = input$axis_text),   # tick labels
                       axis.title = element_text(size = input$axis_title),  # axis titles
                       strip.text = element_text(size = input$strip_text),  # facet titles
                       legend.text = element_text(size = input$legend_text),
                       legend.title = element_text(size = input$legend_title)
                     )
                   
                   saved_plots$ts_cluster_plot[[saved.num]] <- p
                 })
                 
                 heatmap_plot <- eventReactive(input$compute_cluster, {
                   withProgress(message = "Processing...", value = 0, {
                     req(DataClusterReactive())
                     DataCluster <-DataClusterReactive()
                     clusters <- DataCluster$clusters
                     df <- clusters$df %>% dplyr::select(genes, cluster) %>% arrange(cluster)
                     gene_order <- df  %>% pull(genes)
                     mat_ordered <- DataCluster$cluster_rlog[gene_order, ]
                     pheatmap(
                       mat_ordered,
                       scale = "row",
                       cluster_rows = FALSE,
                       cluster_cols = FALSE,
                       annotation_row = df["cluster"],
                       show_rownames = FALSE,
                       main = "Time‑series Heatmap Ordered by Cluster"
                     )
                   })
                 })
                 
                 output$ts_heatmap<- renderPlot({
                   heatmap_plot()
                 })
                 
               })
}



