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
  fit <- lm(y_fit ~ x_fit)  # y on x, matching the plot (x = first item, y = second)
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

# ---- Gene-Target helpers ----
# Pull the first number out of each label (e.g. "4wk" -> 4, "Day 12" -> 12).
# NA when a label has no number.
ct_parse_numbers <- function(x) {
  x <- as.character(x)
  out <- rep(NA_real_, length(x))
  ok <- !is.na(x)
  pos <- regexpr("[0-9]*\\.?[0-9]+", x[ok])
  hit <- pos > 0
  vals <- rep(NA_real_, sum(ok))
  vals[hit] <- as.numeric(regmatches(x[ok], pos))
  neg <- hit & grepl("^\\s*-[0-9.]", x[ok])
  vals[neg] <- -vals[neg]
  out[ok] <- vals
  out
}

# How a metadata column can be used as a target:
#   numeric = numeric column; labels = every label carries a number;
#   text    = only usable as ordered levels.
# levels comes back in a sensible default order (by number when there is one).
ct_attr_info <- function(v) {
  if (is.numeric(v)) {
    lv <- as.character(sort(unique(v[!is.na(v)])))
    return(list(kind = "numeric", levels = lv, msg = "Numeric attribute."))
  }
  vc <- as.character(v)
  vc[vc %in% c("", "NA")] <- NA
  lv <- if (is.factor(v)) intersect(levels(v), vc) else unique(vc[!is.na(vc)])
  nums <- ct_parse_numbers(lv)
  # Labels whose numbers carry different units (e.g. 2mo, 1yr) can't be compared by number:
  # treat them as ordered levels in the metadata's level order.
  units <- unique(tolower(trimws(sub("-?[0-9]*\\.?[0-9]+", "", lv))))
  if (length(lv) > 0 && all(!is.na(nums)) && length(units) > 1) {
    return(list(kind = "text", levels = lv,
                msg = paste0("Text labels with mixed units (", paste(units, collapse = ", "),
                             "): using 'Ordered levels' in the metadata order; reorder below if needed.")))
  }
  if (length(lv) > 0 && all(!is.na(nums))) {
    o <- order(nums)
    list(kind = "labels", levels = lv[o],
         msg = paste0("Text labels with numbers: 'Numeric values' uses the number in each label (e.g. ",
                      lv[o][1], " = ", nums[o][1], ")."))
  } else {
    list(kind = "text", levels = lv,
         msg = "Text levels without numbers: use 'Ordered levels' and set the order below.")
  }
}

# Correlation of every row of X (genes x samples) with y.
# Genes without missing values are done in one matrix step; genes with
# missing values are done one by one on their non-missing samples.
ct_cor_rows <- function(X, y, method = "pearson") {
  r <- rep(NA_real_, nrow(X))
  n <- rep(0L, nrow(X))
  if (length(y) < 3 || length(unique(y)) < 2) return(list(r = r, n = n))
  complete <- rowSums(is.na(X)) == 0
  if (any(complete)) {
    Xc <- X[complete, , drop = FALSE]
    yc <- y
    if (method == "spearman") {
      if (requireNamespace("matrixStats", quietly = TRUE)) {
        Xc <- matrixStats::rowRanks(Xc, ties.method = "average")
      } else {
        Xc <- t(apply(Xc, 1, rank, ties.method = "average"))
      }
      yc <- rank(yc, ties.method = "average")
    }
    Xc <- Xc - rowMeans(Xc)
    yc <- yc - mean(yc)
    den <- sqrt(rowSums(Xc^2) * sum(yc^2))
    rr <- as.vector(Xc %*% yc) / den
    rr[!is.finite(rr)] <- NA
    r[complete] <- pmax(pmin(rr, 1), -1)
    n[complete] <- length(y)
  }
  for (i in which(!complete)) {
    ok <- !is.na(X[i, ])
    n[i] <- sum(ok)
    if (n[i] >= 3) {
      r[i] <- suppressWarnings(stats::cor(X[i, ok], y[ok], method = method))
    }
  }
  list(r = r, n = n)
}

# Two-sided p-value for r from n samples (t distribution, n - 2 df).
ct_cor_p <- function(r, n) {
  df <- n - 2
  t <- r * sqrt(df / pmax(1 - r^2, .Machine$double.eps))
  p <- 2 * stats::pt(-abs(t), df)
  p[is.na(r) | df < 1] <- NA
  p
}

# ---- UI Module ----
correlation_ui <- function(id) {
  ns <- shiny::NS(id)
  fluidRow(
    column(3,
           wellPanel(
             column(width=12,uiOutput(ns("selectGroupSampleCorr"))), 
             radioButtons(ns("correlation_type"), "Correlation Type:", choices = c("Gene-Gene" = "gene", "Sample-Sample" = "sample", "Group-Group" = "group", "Gene-Target" = "target"), inline = TRUE, selected = "gene"),
             conditionalPanel(ns = ns, "input.correlation_type!='target'",
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
             ),
             conditionalPanel(ns = ns, "input.correlation_type=='gene' || input.correlation_type=='target'",
                              radioButtons(ns("gene_label"),label="Select Gene Label",inline = TRUE, choices=c("UniqueID", "Gene.Name"), selected="Gene.Name")
             ),
             conditionalPanel(ns = ns, "input.correlation_type=='target'",
                              radioButtons(ns("target_kind"), "Correlate all genes with:", choices = c("A gene" = "gene", "A metadata attribute" = "attr"), inline = TRUE, selected = "gene"),
                              conditionalPanel(ns = ns, "input.target_kind=='gene'",
                                               selectizeInput(ns("target_gene"), label = "Target Gene", choices = NULL, multiple = FALSE, options = list(placeholder = 'Type to search'))
                              ),
                              conditionalPanel(ns = ns, "input.target_kind=='attr'",
                                               selectInput(ns("target_attr"), label = "Target Attribute", choices = NULL),
                                               radioButtons(ns("target_scale"), "Use the attribute as:", choices = c("Numeric values" = "numeric", "Ordered levels" = "ordered"), inline = TRUE, selected = "numeric"),
                                               span(textOutput(ns("target_attr_info")), style = "color:#555; font-size:13px; font-style:italic"),
                                               conditionalPanel(ns = ns, "input.target_scale=='ordered'",
                                                                selectizeInput(ns("target_levels"), label = "Levels, lowest to highest", choices = NULL, multiple = TRUE, options = list(plugins = list('remove_button'))),
                                                                helpText("To move a level to the end, remove it and add it again. Samples in levels not listed are left out. Ordered levels always use Spearman.")
                                               )
                              ),
                              selectInput(ns("restrict_attr"), label = "Restrict samples to (optional)", choices = c("None" = "none")),
                              conditionalPanel(ns = ns, "input.restrict_attr!='none'",
                                               selectizeInput(ns("restrict_levels"), label = "Keep samples in", choices = NULL, multiple = TRUE, options = list(plugins = list('remove_button')))
                              ),
                              selectInput(ns("strata_attr"), label = "Compute separately for each level of (optional)", choices = c("None" = "none")),
                              tags$hr(style = "border-color: black;")
             ),
             conditionalPanel(ns = ns, "input.correlation_type=='sample'",
                              selectizeInput(ns("sel_sample"),	label="Select At Least 2 Samples",	choices = NULL,	multiple=TRUE),
                              span(textOutput(ns("Selected_samples")), style = "color:red; font-size:15px; font-family:arial; font-style:italic"),
                              tags$hr(style="border-color: black;")
             ),             
             conditionalPanel(ns = ns, "input.correlation_type=='group'",
                              selectizeInput(ns("sel_attribute"), label="Select a Attribute", choices = NULL, multiple = FALSE),
                              selectizeInput(ns("sel_group"),	label="Select At Least 2 Groups",	choices = NULL,	multiple=TRUE),
                              span(textOutput(ns("Selected_groups")), style = "color:red; font-size:15px; font-family:arial; font-style:italic"),
                              tags$hr(style="border-color: black;")
             ),
             radioButtons(ns("correlation_method"), "Correlation Method:", choices = c("Pearson", "Spearman"), inline = TRUE, selected = "Pearson"),
           )
    ),
    column(9,
           conditionalPanel(ns = ns, "input.correlation_type=='target'",
           tabsetPanel(id = ns("target_tabset"),
                       tabPanel(title = "Gene-Target Results",
                                br(),
                                tags$details(
                                  tags$summary(
                                    tags$strong("Click to view Gene-Target Correlation Instructions",
                                                style = "color: #28a745; cursor: pointer;")
                                  ),
                                  tags$div(
                                    style = "padding: 15px; background-color: #f4faf6; border-left: 5px solid #28a745; margin-top: 10px;",
                                    tags$p("Gene-Target correlates every gene with one target across individual samples (not group means). The target is either one gene or one metadata attribute."),
                                    tags$p(tags$strong("Setup Steps:")),
                                    tags$ul(
                                      tags$li(tags$strong("Target gene:"), " finds genes whose expression goes up and down with the target gene's expression."),
                                      tags$li(tags$strong("Target attribute, Numeric values:"), " uses numeric attributes (e.g. Age, dose, a score) as they are. For text labels that carry a number (e.g. 4wk, 8wk, Day 12) the number in each label is used."),
                                      tags$li(tags$strong("Target attribute, Ordered levels:"), " for text levels with an order (e.g. Low, Medium, High) list the levels from lowest to highest; they are coded 1, 2, 3, ... and Spearman is used. A two-level attribute (e.g. WT, KO) works the same way; it then ranks genes much like a two-group test. Unordered categories have no correlation; use the DE modules for those."),
                                      tags$li(tags$strong("Restrict samples to:"), " only uses samples in the chosen levels of an attribute (e.g. Condition = Treated)."),
                                      tags$li(tags$strong("Compute separately for each level of:"), " computes one correlation per level (e.g. one per time point), with its own p-values and adjusted p-values."),
                                      tags$li(tags$strong("Method:"), " 'Pearson' for linear relationships, 'Spearman' for monotonic trends that need not be linear.")
                                    ),
                                    tags$p(tags$strong("Column Definitions:")),
                                    tags$ul(
                                      tags$li(tags$strong("r:"), " correlation coefficient from -1 to 1."),
                                      tags$li(tags$strong("p.value:"), " two-sided test of r = 0 (t distribution with n - 2 degrees of freedom)."),
                                      tags$li(tags$strong("padj:"), " Benjamini-Hochberg adjusted p-value, computed separately within each level."),
                                      tags$li(tags$strong("n:"), " samples used for the gene (samples with missing values are left out). A gene needs at least 4 samples to get an r; levels with fewer than 5 samples get a warning because r is unreliable there."),
                                      tags$li(tags$strong("max_abs_r:"), " when computed per level, the largest absolute r across levels; the table is sorted by it.")
                                    )
                                  )
                                ),
                                br(),
                                actionButton(ns("compute_target"), "Compute/Refresh", style = "color: #0961E3; background-color: #F6E98C; border-color: #2e6da4"),
                                actionButton(ns("save_target_table"), "Save table to output"),
                                br(), br(),
                                span(htmlOutput(ns("target_msg")), style = "color:#333; font-size:14px"),
                                tags$hr(),
                                DT::dataTableOutput(ns("target_table"))),
                       tabPanel(title = "Gene-Target Plot",
                                actionButton(ns("save_target_plot"), "Save to output"),
                                tags$br(),
                                tags$hr(),
                                plotly::plotlyOutput(ns("target_plot"), width = "900px", height = "650px")),
                       tabPanel(title = "Correlation Heatmap",
                                fluidRow(
                                  column(width = 4, sliderInput(ns("target_hm_n"), "Top genes to show", min = 10, max = 100, value = 30, step = 5)),
                                  column(width = 4, br(), actionButton(ns("save_target_heatmap"), "Save to output"))
                                ),
                                tags$hr(),
                                morpheus::morpheusOutput(ns("target_heatmap"), width = "800px", height = "900px"))
           )),
           conditionalPanel(ns = ns, "input.correlation_type!='target'",
           tabsetPanel(id = ns("correlation_tabset"),
                       tabPanel(title = "Result Table", 
                                br(),
                                tags$details(
                                  tags$summary(
                                    tags$strong("Click to view Correlation Analysis Instructions", 
                                                style = "color: #28a745; cursor: pointer;")
                                  ),
                                  tags$div(
                                    style = "padding: 15px; background-color: #f4faf6; border-left: 5px solid #28a745; margin-top: 10px;",
                                    tags$p("This tool calculates the statistical relationship between different biological variables (Gene_Gene, Sample_Sample, Group-Group). Use it to identify genes that function together or respond similarly across samples, or to identify samples or groups that function together or respond similarly across genes."),
                                    tags$p(tags$strong("Setup Steps:")),
                                    tags$ul(
                                      tags$li(tags$strong("Correlation Type:"), " Choose  one from 'Gene-Gene', 'Sample-Sample' and 'Group-Group'. To correlate all genes with one gene or one metadata attribute, choose 'Gene-Target' (it has its own instructions)."),
                                      tags$li(tags$strong("Gene Selection:"), " Type or select at least 2 genes no matter what method is used to input the gene list"),
                                      tags$li(tags$strong("Sample or Group Selection:"), " Select at least 2 samples or groups."),
                                      tags$li(tags$strong("Method:"), " Use 'Pearson' for linear relationships or 'Spearman' if you expect non-linear but monotonic trends.")
                                    ),
                                    tags$p(tags$strong("Interpreting Results:"), " After clicking 'Compute/Refresh', the ranked ", tags$strong("Correlation Result Table"), " will display. It displays pairwise statistical relationships between the selected genes, samples or groups. Each row represents a unique comparison."),
                                    tags$p(tags$strong("Column Definitions:")),
                                    tags$ul(
                                      tags$li(tags$strong("Correlation:"), " The coefficient (r) ranging from -1 to 1. Values near 1 indicate a strong positive relationship, while values near -1 indicate an inverse relationship."),
                                      tags$li(tags$strong("R-squared:"), " Represents the proportion of variance shared between the two variables. Higher values suggest a more predictable relationship."),
                                      tags$li(tags$strong("P-value:"), " Indicates statistical significance. A value below 0.05 typically suggests the correlation is not due to random chance.")
                                    ),
                                    tags$p(
                                      tags$strong("Action Column:"), 
                                      " Click ", tags$em("View Correlation Plot"), " in any row to generate a scatter plot and regression line for that specific pair."
                                    )
                                  )
                                ),
                                tags$div(
                                  style = "padding: 15px; margin-top: 20px; margin-bottom: 20px; border-left: 5px solid #28a745; background-color: #f4faf6; border-radius: 4px;",
                                  tags$strong("Note: "),
                                  "If you want to check the global sample-sample correlation for all the samples, please go to the ",
                                  actionLink(ns("go_to_qc"), "sample correlation Heatmap"),
                                  " in QC Plots module."
                                ),
                                br(),
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
           ))
    )
  )}

# ---- Server Module ----
correlation_server <- function(id, parent_session) {
  moduleServer(id,
               function(input, output, session) {
                 # actionButton click counts get bookmarked like any other input --
                 # restoring a non-zero count silently re-triggers that button's own
                 # observeEvent right after restore (Shiny can't tell "value changed
                 # because it was restored" from "value changed because it was
                 # clicked"). None of these buttons' click counts are meaningful
                 # state to restore, so exclude them all.
                 session$setBookmarkExclude(c("go_to_qc", "compute_corr", "CorrPlot",
                                               "compute_target", "save_target_table", "save_target_plot", "save_target_heatmap"))
                 ProteinGeneName_sel <- reactiveVal()
                 CorrPlot <- reactiveVal()
                 CorrPlot_ID <- reactiveVal()
                 genelabel <- reactiveVal()
                 CorrMethod <- reactiveVal()

                 # sel_test/sel_sample/sel_attribute/sel_group/sel_gene/
                 # gene_subset/sel_geneset all have their choices+selected
                 # forced back to a hardcoded default by observers that fire
                 # on every project load (see below) -- including at restore
                 # time. Reapply unconditionally for a few seconds after
                 # restore, same pattern as app.R.
                 # Attributes offered for target/restrict/strata (same rule as
                 # the attribute-choices observer below): more than one
                 # non-missing value.
                 ct_restore_attrs <- function() {
                   req(DataQCReactive())
                   MetaData <- DataQCReactive()$MetaData
                   req(MetaData)
                   attrs <- sort(setdiff(colnames(MetaData), c("sampleid", "Order", "ComparePairs")))
                   attrs[vapply(attrs, function(a) {
                     v <- as.character(MetaData[[a]])
                     length(unique(v[!is.na(v) & !(v %in% c("", "NA"))])) > 1
                   }, logical(1))]
                 }

                 session$onRestored(function(state) {
                   restored <- state$input
                   pending <- list(
                     sel_test = list(value = restored$sel_test, apply = function(v) {
                       updateSelectizeInput(session, "sel_test", choices = test_order(), selected = v)
                     }),
                     sel_sample = list(value = restored$sel_sample, apply = function(v) {
                       updateSelectizeInput(session, "sel_sample", choices = sample_order(), selected = v)
                     }),
                     sel_attribute = list(value = restored$sel_attribute, apply = function(v) {
                       req(DataQCReactive())
                       attrs <- sort(setdiff(colnames(DataQCReactive()$MetaData), c("sampleid", "Order", "ComparePairs")))
                       updateSelectInput(session, "sel_attribute", choices = attrs, selected = v)
                     }),
                     # sel_group's choices depend on whichever attribute is
                     # currently selected -- read isolate(input$sel_attribute)
                     # fresh on every reapply, same as sel_group's own
                     # observeEvent(input$sel_attribute, ...) does live.
                     sel_group = list(value = restored$sel_group, apply = function(v) {
                       req(DataQCReactive())
                       attr <- isolate(input$sel_attribute)
                       req(attr)
                       DataIn <- DataQCReactive()
                       if (attr %in% names(DataIn$tmp_group)) {
                         groups <- DataIn$tmp_group[[attr]]
                       } else {
                         groups <- all_group_list()[[attr]]
                       }
                       updateSelectizeInput(session, "sel_group", choices = groups, selected = v)
                     }),
                     # Must re-supply the full choices list on every apply --
                     # a server=TRUE selectize update with no choices
                     # registers an empty searchable dataset (see sel_gene in
                     # app.R for the same gotcha).
                     sel_gene = list(value = restored$sel_gene, apply = function(v) {
                       req(DataQCReactive())
                       req(isolate(input$gene_label))
                       ProteinGeneName <- DataQCReactive()$ProteinGeneName
                       req(ProteinGeneName)
                       if (isolate(input$gene_label) == "UniqueID") {
                         DataIngenes <- ProteinGeneName %>% dplyr::select(UniqueID) %>% collect %>% .[["UniqueID"]] %>% as.character()
                       } else {
                         DataIngenes <- ProteinGeneName %>% dplyr::select(Gene.Name) %>% collect %>% .[["Gene.Name"]] %>% as.character()
                       }
                       updateSelectizeInput(session, "sel_gene", choices = DataIngenes, selected = v, server = TRUE)
                     }),
                     # gene_subset's choices depend on correlation_type --
                     # same cascading pattern as sel_group/sel_attribute.
                     gene_subset = list(value = restored$gene_subset, apply = function(v) {
                       if (identical(isolate(input$correlation_type), "gene")) {
                         choices <- c("Select", "Browsing", "Upload Genes", "Geneset")
                       } else {
                         choices <- c("All", "Browsing", "Upload Genes", "Geneset")
                       }
                       updateRadioButtons(session, "gene_subset", choices = choices, inline = TRUE, selected = v)
                     }),
                     sel_geneset = list(value = restored$sel_geneset, apply = function(v) {
                       genesetnames <- GetGeneSetNames()
                       updateSelectizeInput(session, "sel_geneset", choices = c('Type to Search' = '', genesetnames), selected = v, server = TRUE)
                     }),
                     # Gene vs. attribute inputs: target_attr/restrict_attr/
                     # strata_attr are reset by the attribute-choices observer
                     # on every project load; target_scale/target_levels and
                     # restrict_levels are reset whenever their attribute
                     # changes (including when restored above), so they read
                     # the current attribute fresh on every reapply.
                     target_attr = list(value = restored$target_attr, apply = function(v) {
                       updateSelectInput(session, "target_attr", choices = ct_restore_attrs(), selected = v)
                     }),
                     target_scale = list(value = restored$target_scale, apply = function(v) {
                       updateRadioButtons(session, "target_scale", selected = v)
                     }),
                     target_levels = list(value = restored$target_levels, apply = function(v) {
                       req(DataQCReactive())
                       attr <- isolate(input$target_attr)
                       MetaData <- DataQCReactive()$MetaData
                       req(attr, attr %in% colnames(MetaData))
                       info <- ct_attr_info(MetaData[[attr]])
                       updateSelectizeInput(session, "target_levels", choices = info$levels, selected = v)
                     }),
                     restrict_attr = list(value = restored$restrict_attr, apply = function(v) {
                       updateSelectInput(session, "restrict_attr", choices = c("None" = "none", ct_restore_attrs()), selected = v)
                     }),
                     restrict_levels = list(value = restored$restrict_levels, apply = function(v) {
                       req(DataQCReactive())
                       attr <- isolate(input$restrict_attr)
                       MetaData <- DataQCReactive()$MetaData
                       req(attr, attr != "none", attr %in% colnames(MetaData))
                       x <- MetaData[[attr]]
                       lv <- if (is.factor(x)) levels(droplevels(x)) else sort(unique(as.character(x)))
                       lv <- lv[!is.na(lv) & !(lv %in% c("", "NA"))]
                       updateSelectizeInput(session, "restrict_levels", choices = lv, selected = v)
                     }),
                     strata_attr = list(value = restored$strata_attr, apply = function(v) {
                       updateSelectInput(session, "strata_attr", choices = c("None" = "none", ct_restore_attrs()), selected = v)
                     }),
                     # Must re-supply the full choices (server=TRUE), same as sel_gene.
                     target_gene = list(value = restored$target_gene, apply = function(v) {
                       req(DataQCReactive())
                       PGN <- DataQCReactive()$ProteinGeneName
                       req(PGN, isolate(input$gene_label))
                       uid <- as.character(PGN$UniqueID)
                       lab <- if (isolate(input$gene_label) == "UniqueID") uid else as.character(PGN$Gene.Name)
                       dup <- duplicated(lab) | duplicated(lab, fromLast = TRUE)
                       lab[dup] <- paste0(lab[dup], " (", uid[dup], ")")
                       updateSelectizeInput(session, "target_gene", choices = stats::setNames(uid, lab), selected = v, server = TRUE)
                     })
                   )
                   pending <- Filter(function(p) !is.null(p$value), pending)
                   if (length(pending) > 0) {
                     attempts_left <- 10  # ~3s at 300ms
                     restore_observer <- NULL
                     restore_observer <- observe({
                       invalidateLater(300, session)
                       isolate({
                         attempts_left <<- attempts_left - 1
                         for (nm in names(pending)) {
                           p <- pending[[nm]]
                           tryCatch(p$apply(p$value), error = function(e) NULL)
                         }
                         if (attempts_left <= 0) restore_observer$destroy()
                       })
                     })
                   }
                 })

                 observeEvent(input$go_to_qc, {
                   updateNavbarPage(parent_session, inputId = "menu", selected = "QC_Plots")
                   updateTabsetPanel(parent_session, inputId = "groupplot_tabset", selected = "Sample-sample Distance")
                   updateTabsetPanel(parent_session, inputId = "sd_tabset", selected = "Correlation Heatmap")
                 })
                 
                 
                 observe({
                   req(DataQCReactive())
                   req(input$gene_label)
                   DataIn = DataQCReactive()
                   ProteinGeneName = DataIn$ProteinGeneName
                   req(ProteinGeneName)
                   MetaData = DataIn$MetaData
                   
                   if (input$gene_label=="UniqueID") {
                     DataIngenes <- ProteinGeneName %>% dplyr::select(UniqueID) %>% collect %>% .[["UniqueID"]] %>%	as.character()
                   } else {
                     DataIngenes <- ProteinGeneName %>% dplyr::select(Gene.Name) %>% collect %>% .[["Gene.Name"]] %>%	as.character()
                   }
                   updateSelectizeInput(session,'sel_gene',choices=DataIngenes, server = TRUE)
                 })
                 
                 observe({
                   tests = test_order()
                   ProteinGeneName_Header = ProteinGeneNameHeader()
                   updateSelectizeInput(session,'sel_test',choices=tests, selected=tests[1])
                 })
                 
                 observe({
                   preset_samples=sample_order()
                   updateSelectizeInput(session,'sel_sample',choices=preset_samples)
                 })
                 
                 observe({
                   DataIn = DataQCReactive()
                   MetaData = DataIn$MetaData
                   if ('group' %in% names(DataIn$tmp_group)) {
                     preset_group <- as.factor(DataIn$tmp_group$group)
                   } else {
                     preset_group <- as.factor(all_group_list()$group)
                   }
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
                 
                 observe({
                   req(input$gene_subset == "Select")
                   req(DataQCReactive())
                   gene_list = input$sel_gene
                   if (length(gene_list) <= 1) log_warn("Correlation Analysis: fewer than 2 genes selected.")
                   validate(need(length(gene_list) > 1, message = "Please input at least 2 matched genes."))
                   DataIn = DataQCReactive()
                   ProteinGeneName = DataIn$ProteinGeneName
                   res_gene_count <- get_gene_counts(ProteinGeneName, gene_list)
                   ProteinGeneName_sel(res_gene_count$ProteinGeneName_sel)
                   genelabel(input$gene_label)
                   if (res_gene_count$gene_label_reset) genelabel("UniqueID")
                   output$filteredgene_Select <- renderText({res_gene_count$msg_filter1})
                   output$filteredUniqueID_Select <- renderText({res_gene_count$msg_filter})
                 })
                 
                 observe({
                   req(input$gene_subset == "Browsing")
                   req(DataQCReactive())
                   req(input$sel_test != "")
                   p_sel   <- input$psel
                   test_sel <- input$sel_test
                   FCcut <- log2(as.numeric(input$fccut))
                   pvalcut <- as.numeric(input$pvalcut)
                   Updown <- input$updown
                   
                   DataIn = DataQCReactive()
                   results_long = DataIn$tmp_results_long
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
                 
                 observe({
                   # Plain observe() (like the Browsing block above), not
                   # observeEvent(input$gene_list, ...) -- that only fired
                   # when the uploaded gene text itself changed, so
                   # switching "Genes Used in Correlation Analysis" from
                   # Browsing to Upload Genes without touching the text box
                   # left ProteinGeneName_sel() holding the stale Browsing
                   # selection, and Compute/Refresh kept producing the same
                   # result. This re-runs on every switch into "Upload
                   # Genes" too, since it reads input$gene_subset below.
                   req(input$gene_subset == "Upload Genes")
                   req(DataQCReactive())
                   gene_list <- input$gene_list
                   gene_list <- ProcessUploadGeneList(gene_list)
                   if (length(gene_list) <= 1) log_warn("Correlation Analysis: uploaded gene list has fewer than 2 matched genes.")
                   validate(need(length(gene_list) > 1, message = "Please input at least 2 matched genes."))
                   
                   DataIn = DataQCReactive()
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
                   
                   DataIn = DataQCReactive()
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
                   if (!is.null(sample_list) && length(sample_list) > 1) {
                     msg_sample <- paste("Selected Samples:",length(sample_list),sep="")
                   } else {
                     msg_sample <- "Please select at least 2 samples."
                   }
                   output$Selected_samples <- renderText({msg_sample})
                 })
                 
                 observeEvent(input$sel_attribute, {
                   req(input$sel_attribute)
                   DataIn = DataQCReactive()
                   MetaData = DataIn$MetaData # %>% filter(sampleid %in% sample_order())
                   # groups = unique(MetaData[, input$sel_attribute])
                   if (input$sel_attribute %in% names(DataIn$tmp_group)) {
                     groups <- DataIn$tmp_group[[input$sel_attribute]]
                   } else {
                     groups <- all_group_list()[[input$sel_attribute]]
                   }
                   
                   updateSelectizeInput(session,'sel_group',choices=groups)
                   output$Selected_groups <- renderText({""})
                 })
                 
                 observeEvent(input$sel_group, {
                   group_list = input$sel_group
                   if (!is.null(group_list) && length(group_list) > 1) {
                     msg_sample <- paste("Selected Groups:",length(group_list),sep="")
                   } else {
                     msg_sample <- "Please select at least 2 groups."
                   }
                   
                   msg_group <- paste("Selected Groups:",length(group_list),sep="")
                   output$Selected_groups <- renderText({msg_group})
                 })
                 
                 CorrResult <- eventReactive(input$compute_corr, {
                   log_info("Correlation Analysis: Compute/Refresh clicked (type=", input$correlation_type, ").")
                   DataIn <- DataQCReactive()
                   data_long <- DataIn$tmp_data_long
                   # preset_group=group_order()
                   # preset_samples=sample_order()
                   if (input$gene_subset == "All") {
                     tmpids <- DataIn$ProteinGeneName %>% dplyr::pull(UniqueID) %>% unique()
                   } else {
                     ProteinGeneName_sel <- ProteinGeneName_sel()
                     genelabel <- genelabel()
                     if (is.null(ProteinGeneName_sel) || nrow(ProteinGeneName_sel) <= 1) log_warn("Correlation Analysis: fewer than 2 matched genes selected.")
                     validate(need(nrow(ProteinGeneName_sel) > 1, message = "Please input at least 2 matched genes."))
                     tmpids <- ProteinGeneName_sel %>% dplyr::pull(UniqueID) %>% unique()
                   }
                   if (input$correlation_type=='gene') {
                     exp_tmp = data_long %>% 
                       dplyr::filter(UniqueID %in% tmpids) %>%        #, group %in% preset_group, sampleid %in% preset_samples) %>%
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
                     if (length(input$sel_sample) <= 1) log_warn("Correlation Analysis: fewer than 2 samples selected.")
                     validate(need(length(input$sel_sample) > 1, message = "Please input at least 2 samples."))
                     exp_tmp = data_long %>% 
                       dplyr::filter(UniqueID %in% tmpids, sampleid %in% input$sel_sample) %>%        # , group %in% preset_group, sampleid %in% input$sel_sample) %>%
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
                     #browser()
                     sel_group <- input$sel_group
                     if (length(sel_group) <= 1) log_warn("Correlation Analysis: fewer than 2 groups selected.")
                     validate(need(length(sel_group) > 1, message = "Please input at least 2 groups."))
                     adding_number <- ifelse(exp_unit() == "Expression Level", 0, as.numeric(str_extract(exp_unit(), "(?<=\\+)\\d*\\.?\\d+")))
                     if (is.na(adding_number)) {adding_number=0}
                     # data_long <- DataIn$data_long
                     sel_attribute <- input$sel_attribute
                     MetaData = DataIn$MetaData # %>% filter(sampleid %in% sample_order(), !!sym(sel_attribute) %in% sel_group)
                     
                     exp_tmp = data_long %>%
                       # dplyr::left_join(MetaData %>% dplyr::select(sampleid, !!sym(sel_attribute)), by = "sampleid", suffix = c("", "_meta")) %>%
                       dplyr::filter(UniqueID %in% tmpids) %>%       # , sampleid %in% MetaData$sampleid) %>%
                       dplyr::select(gene = UniqueID, !!sym(sel_attribute), expr) %>%
                       # Restrict to the groups the user actually selected --
                       # without this, every level of sel_attribute (not just
                       # the selected ones) ends up as a row after
                       # pivot_wider below, so toCheck_list/combn() further
                       # down produces every pairwise combination of ALL
                       # groups for that attribute, not just the selected
                       # ones (confirmed directly: selecting 2 of 3 Genotype
                       # levels still returned all 3 pairwise rows).
                       dplyr::filter(!!sym(sel_attribute) %in% sel_group) %>%
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
                   toCheck_list <- toCheck_list[toCheck_list != ""]
                   
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
                   
                   res_corr_sorted$Action <- "View Correlation Plot"
                   res_corr_sorted <- res_corr_sorted %>%
                     dplyr::relocate(Action, .before = intercept) %>%
                     dplyr::relocate(rank, .before = 1)

                   log_info("Correlation Analysis: computed ", nrow(res_corr_sorted), " pairwise correlation(s).")
                   return(list("corr_table" = res_corr_sorted, "exp_filtered" = exp_tmp ))
                 })
                 
                 output$corr_table <-  DT::renderDT({
                   res <- CorrResult()
                   corr_table <- res$corr_table
                   
                   DT::datatable(corr_table,  
                                 rownames = FALSE, extensions = 'Buttons', escape = FALSE, selection = 'none', class = 'cell-border strip hover',
                                 options = list(dom = 'lBfrtip', pageLength = 15,
                                                buttons = list(
                                                  list(extend = "csv", text = "Download Page", filename = "Page_results",
                                                       exportOptions = list(modifier = list(page = "current"))),
                                                  list(extend = "csv", text = "Download All", filename = "All_Results",
                                                       exportOptions = list(modifier = list(page = "all")))
                                                )
                                 )) %>% 
                     formatSignif(columns=c('intercept', 'slope', 'r_squared', 'correlation', 'p.value'), digits=3) %>% 
                     formatStyle(4, cursor = 'pointer',color='blue')
                 })
                 
                 observeEvent(input$corr_table_cell_clicked, {
                   clicked <- input$corr_table_cell_clicked
                   res <- CorrResult()
                   corr_table <- res$corr_table
                   exp_filtered <- res$exp_filtered
                   corr_method = CorrMethod()
                   
                   if (length(clicked) > 0) {
                     row <- clicked$row
                     gene1 <- corr_table[row, 2]
                     gene2 <- corr_table[row, 3]
                     
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
                   log_info("Correlation Plot (", CorrPlot_ID(), ") saved to output.")
                   saved_plots$CorrPlot[[CorrPlot_ID()]] <- CorrPlot()
                 })
                 
                 observeEvent(input$corr_table_cell_clicked, {
                   if (length(input$corr_table_cell_clicked) > 0) {
                     updateTabsetPanel(session, "correlation_tabset", selected = "Correlation Plot")
                   }
                 })
                 # ---- Gene-Target ----
                 # Attribute choices: metadata columns with more than one value.
                 observe({
                   req(DataQCReactive())
                   MetaData <- DataQCReactive()$MetaData
                   req(MetaData)
                   attrs <- sort(setdiff(colnames(MetaData), c("sampleid", "Order", "ComparePairs")))
                   attrs <- attrs[vapply(attrs, function(a) {
                     v <- as.character(MetaData[[a]])
                     length(unique(v[!is.na(v) & !(v %in% c("", "NA"))])) > 1
                   }, logical(1))]
                   cur <- isolate(input$target_attr)
                   updateSelectInput(session, "target_attr", choices = attrs,
                                     selected = if (isTruthy(cur) && cur %in% attrs) cur else attrs[1])
                   cur <- isolate(input$restrict_attr)
                   updateSelectInput(session, "restrict_attr", choices = c("None" = "none", attrs),
                                     selected = if (isTruthy(cur) && cur %in% attrs) cur else "none")
                   cur <- isolate(input$strata_attr)
                   updateSelectInput(session, "strata_attr", choices = c("None" = "none", attrs),
                                     selected = if (isTruthy(cur) && cur %in% attrs) cur else "none")
                 })

                 # Target gene choices follow the gene label; values are UniqueIDs.
                 observe({
                   req(DataQCReactive())
                   req(input$gene_label)
                   PGN <- DataQCReactive()$ProteinGeneName
                   req(PGN)
                   uid <- as.character(PGN$UniqueID)
                   lab <- if (input$gene_label == "UniqueID") uid else as.character(PGN$Gene.Name)
                   dup <- duplicated(lab) | duplicated(lab, fromLast = TRUE)
                   lab[dup] <- paste0(lab[dup], " (", uid[dup], ")")
                   updateSelectizeInput(session, "target_gene", choices = stats::setNames(uid, lab),
                                        selected = isolate(input$target_gene), server = TRUE)
                 })

                 target_attr_info <- reactive({
                   req(DataQCReactive())
                   req(input$target_attr)
                   MetaData <- DataQCReactive()$MetaData
                   req(input$target_attr %in% colnames(MetaData))
                   ct_attr_info(MetaData[[input$target_attr]])
                 })

                 observeEvent(target_attr_info(), {
                   info <- target_attr_info()
                   updateRadioButtons(session, "target_scale", selected = if (info$kind == "text") "ordered" else "numeric")
                   updateSelectizeInput(session, "target_levels", choices = info$levels, selected = info$levels)
                 })

                 output$target_attr_info <- renderText({
                   info <- target_attr_info()
                   if (info$kind == "text" && identical(input$target_scale, "numeric")) {
                     paste(info$msg, "'Numeric values' is not available for this attribute.")
                   } else {
                     info$msg
                   }
                 })

                 observeEvent(input$restrict_attr, {
                   req(input$restrict_attr != "none")
                   MetaData <- DataQCReactive()$MetaData
                   req(input$restrict_attr %in% colnames(MetaData))
                   v <- MetaData[[input$restrict_attr]]
                   lv <- if (is.factor(v)) levels(droplevels(v)) else sort(unique(as.character(v)))
                   lv <- lv[!is.na(lv) & !(lv %in% c("", "NA"))]
                   updateSelectizeInput(session, "restrict_levels", choices = lv, selected = lv)
                 })

                 gene_display_label <- function(uid, PGN) {
                   if (identical(input$gene_label, "UniqueID")) return(uid)
                   lab <- as.character(PGN$Gene.Name[match(uid, PGN$UniqueID)])
                   ifelse(is.na(lab) | lab == "", uid, lab)
                 }

                 TargetResult <- eventReactive(input$compute_target, {
                   withProgress(message = "Computing Gene-Target correlations...", value = 0, {
                   req(DataQCReactive())
                   DataIn <- DataQCReactive()
                   X <- DataIn$tmp_data_wide
                   MetaData <- DataIn$MetaData
                   PGN <- DataIn$ProteinGeneName
                   validate(need(!is.null(X) && !is.null(MetaData), "No sample data is loaded."))
                   samples <- colnames(X)
                   md <- MetaData[match(samples, as.character(MetaData$sampleid)), , drop = FALSE]
                   method <- tolower(input$correlation_method)
                   x_levels <- NULL

                   if (input$target_kind == "gene") {
                     validate(need(isTruthy(input$target_gene) && input$target_gene %in% rownames(X),
                                   "Please select a target gene."))
                     target_uid <- input$target_gene
                     y <- as.numeric(X[target_uid, ])
                     target_name <- gene_display_label(target_uid, PGN)
                     target_axis <- paste(target_name, "expression")
                   } else {
                     attr <- input$target_attr
                     validate(need(isTruthy(attr) && attr %in% colnames(md), "Please select a target attribute."))
                     info <- ct_attr_info(MetaData[[attr]])
                     v <- md[[attr]]
                     vc <- as.character(v)
                     vc[vc %in% c("", "NA")] <- NA
                     if (input$target_scale == "numeric") {
                       validate(need(info$kind != "text",
                                     paste0("'", attr, "' has no numbers in its values. Use 'Ordered levels' instead.")))
                       y <- if (is.numeric(v)) as.numeric(v) else ct_parse_numbers(vc)
                     } else {
                       lv <- input$target_levels
                       validate(need(length(lv) >= 2, "Please list at least 2 levels, lowest to highest."))
                       y <- as.numeric(match(vc, lv))
                       x_levels <- lv
                       method <- "spearman"
                     }
                     target_uid <- NULL
                     target_name <- attr
                     target_axis <- attr
                   }

                   keep <- !is.na(y)
                   restrict_note <- ""
                   if (isTruthy(input$restrict_attr) && input$restrict_attr != "none") {
                     validate(need(length(input$restrict_levels) > 0, "Please choose at least one level under 'Keep samples in'."))
                     keep <- keep & as.character(md[[input$restrict_attr]]) %in% input$restrict_levels
                     restrict_note <- paste0(" Samples restricted to ", input$restrict_attr, " = ",
                                             paste(input$restrict_levels, collapse = ", "), ".")
                   }

                   strata_attr <- input$strata_attr
                   if (isTruthy(strata_attr) && strata_attr != "none") {
                     validate(need(input$target_kind == "gene" || strata_attr != input$target_attr,
                                   "The target attribute and 'Compute separately for each level of' must be different."))
                     sv <- md[[strata_attr]]
                     strata_vec <- as.character(sv)
                     strata_vec[strata_vec %in% c("", "NA")] <- NA
                     keep <- keep & !is.na(strata_vec)
                     strata <- if (is.factor(sv)) intersect(levels(sv), strata_vec[keep]) else sort(unique(strata_vec[keep]))
                   } else {
                     strata_attr <- NULL
                     strata_vec <- rep("All", length(samples))
                     strata <- "All"
                   }

                   if (!is.null(target_uid)) X <- X[rownames(X) != target_uid, , drop = FALSE]
                   uid <- rownames(X)
                   stats_list <- list()
                   notes <- character(0)
                   for (k in seq_along(strata)) {
                     s <- strata[k]
                     idx <- which(keep & strata_vec == s)
                     lab <- if (is.null(strata_attr)) "" else paste0(strata_attr, " = ", s, ": ")
                     if (length(idx) < 4) {
                       notes <- c(notes, paste0(lab, "skipped, only ", length(idx), " sample(s) (at least 4 needed)."))
                       next
                     }
                     if (length(unique(y[idx])) < 2) {
                       notes <- c(notes, paste0(lab, "skipped, the target has the same value in all ", length(idx), " samples."))
                       next
                     }
                     if (length(idx) < 5) notes <- c(notes, paste0(lab, "only ", length(idx), " samples, r is unreliable."))
                     cr <- ct_cor_rows(X[, idx, drop = FALSE], y[idx], method)
                     cr$r[cr$n < 4] <- NA
                     p <- ct_cor_p(cr$r, cr$n)
                     stats_list[[s]] <- data.frame(r = cr$r, p.value = p, padj = stats::p.adjust(p, method = "BH"), n = cr$n)
                     incProgress(1 / length(strata))
                   }
                   validate(need(length(stats_list) > 0, paste(c("No correlations could be computed.", notes), collapse = " ")))
                   used <- names(stats_list)

                   tbl <- data.frame(UniqueID = uid, Gene.Name = as.character(PGN$Gene.Name[match(uid, PGN$UniqueID)]),
                                     check.names = FALSE, stringsAsFactors = FALSE)
                   if (length(used) == 1 && is.null(strata_attr)) {
                     tbl <- cbind(tbl, stats_list[[1]])
                     tbl <- tbl[!is.na(tbl$r), , drop = FALSE]
                     tbl <- tbl[order(tbl$p.value, -abs(tbl$r)), , drop = FALSE]
                   } else {
                     for (s in used) {
                       st <- stats_list[[s]]
                       names(st) <- paste0(s, "_", names(st))
                       tbl <- cbind(tbl, st)
                     }
                     rmat <- as.matrix(tbl[, paste0(used, "_r"), drop = FALSE])
                     tbl$max_abs_r <- suppressWarnings(apply(abs(rmat), 1, max, na.rm = TRUE))
                     tbl <- tbl[is.finite(tbl$max_abs_r), , drop = FALSE]
                     tbl <- tbl[order(-tbl$max_abs_r), , drop = FALSE]
                   }
                   rownames(tbl) <- NULL

                   n_used <- sum(keep & strata_vec %in% used)
                   method_lab <- if (method == "spearman") "Spearman" else "Pearson"
                   msg <- paste0(method_lab, " correlation of ", format(nrow(tbl), big.mark = ","), " genes with ",
                                 if (input$target_kind == "gene") paste("gene", target_name) else target_name,
                                 " across ", n_used, " samples",
                                 if (is.null(strata_attr)) "." else paste0(", computed separately for each level of ", strata_attr, "."),
                                 restrict_note)
                   if (!is.null(x_levels)) msg <- paste0(msg, " Ordered levels: ", paste(x_levels, collapse = " < "), ".")
                   msg <- paste0(msg, "<br><span style='color:#0961E3'>The table includes all genes. Click a row to show its scatter plot in the Gene-Target Plot tab.",
                                 if (is.null(strata_attr)) "" else " The Correlation Heatmap tab shows the top rows of this table as currently sorted and filtered.",
                                 "</span>")
                   if (length(notes) > 0) msg <- paste0(msg, "<br><span style='color:red'>", paste(notes, collapse = "<br>"), "</span>")
                   log_info("Correlation Analysis: Gene-Target computed (target=", target_name, ", method=", method,
                            ", strata=", if (is.null(strata_attr)) "none" else strata_attr, ", genes=", nrow(tbl), ").")

                   list(table = tbl, strata = used, strata_attr = strata_attr, strata_vec = strata_vec, keep = keep,
                        y = y, X = DataIn$tmp_data_wide, samples = samples, PGN = PGN, method = method_lab,
                        target_name = target_name, target_axis = target_axis, x_levels = x_levels, msg = msg)
                   })
                 })

                 output$target_msg <- renderUI({
                   HTML(TargetResult()$msg)
                 })

                 output$target_table <- DT::renderDT({
                   tbl <- TargetResult()$table
                   num_cols <- which(vapply(tbl, is.double, logical(1)))
                   DT::datatable(tbl, rownames = FALSE, extensions = 'Buttons', selection = 'single',
                                 class = 'cell-border strip hover',
                                 options = list(dom = 'lBfrtip', pageLength = 15, scrollX = TRUE,
                                                buttons = list(
                                                  list(extend = "csv", text = "Download Page", filename = "Page_results",
                                                       exportOptions = list(modifier = list(page = "current"))),
                                                  list(extend = "csv", text = "Download All", filename = "All_Results",
                                                       exportOptions = list(modifier = list(page = "all")))
                                                ))) %>%
                     formatSignif(columns = num_cols, digits = 3)
                 })

                 # Row shown in the plot: the clicked row, or the top row after each Compute.
                 target_row <- reactiveVal(1)
                 observeEvent(TargetResult(), target_row(1))
                 observeEvent(input$target_table_cell_clicked, {
                   clicked <- input$target_table_cell_clicked
                   if (length(clicked) > 0 && !is.null(clicked$row)) {
                     target_row(clicked$row)
                     updateTabsetPanel(session, "target_tabset", selected = "Gene-Target Plot")
                   }
                 })

                 target_selected_uid <- reactive({
                   tbl <- TargetResult()$table
                   validate(need(nrow(tbl) > 0, "No genes in the result table."))
                   row <- min(target_row(), nrow(tbl))
                   tbl$UniqueID[row]
                 })

                 target_plot <- reactive({
                   res <- TargetResult()
                   uid <- target_selected_uid()
                   gene_lab <- gene_display_label(uid, res$PGN)
                   df <- data.frame(x = res$y, expr = as.numeric(res$X[uid, res$samples]),
                                    stratum = res$strata_vec, sample = res$samples, stringsAsFactors = FALSE)
                   df <- df[res$keep & df$stratum %in% res$strata & !is.na(df$expr), , drop = FALSE]
                   # Hover text for the interactive plot.
                   x_txt <- if (!is.null(res$x_levels)) res$x_levels[df$x] else signif(df$x, 4)
                   df$hover <- paste0("Sample: ", df$sample, "<br>", res$target_axis, ": ", x_txt,
                                      "<br>", gene_lab, ": ", signif(df$expr, 4),
                                      if (is.null(res$strata_attr)) "" else paste0("<br>", res$strata_attr, ": ", df$stratum))
                   tbl <- res$table[res$table$UniqueID == uid, , drop = FALSE]
                   r_lab <- function(s) {
                     pre <- if (length(res$strata) == 1 && is.null(res$strata_attr)) "" else paste0(s, "_")
                     sprintf("r = %.2f, padj = %.2g", tbl[[paste0(pre, "r")]], tbl[[paste0(pre, "padj")]])
                   }
                   if (is.null(res$strata_attr)) {
                     df$stratum <- r_lab("All")
                     p <- ggplot(df, aes(x, expr)) + geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "steelblue")
                     colour_title <- NULL
                   } else {
                     labs_map <- stats::setNames(paste0(res$strata, " (", vapply(res$strata, r_lab, character(1)), ")"), res$strata)
                     df$stratum <- factor(labs_map[df$stratum], levels = labs_map)
                     p <- ggplot(df, aes(x, expr, colour = stratum)) + geom_smooth(method = "lm", formula = y ~ x, se = FALSE)
                     colour_title <- res$strata_attr
                   }
                   if (!is.null(res$x_levels)) {
                     p <- p + suppressWarnings(geom_jitter(aes(text = hover), width = 0.1, height = 0, size = 2.5, alpha = 0.8)) +
                       scale_x_continuous(breaks = seq_along(res$x_levels), labels = res$x_levels)
                   } else if (length(unique(df$x)) <= 10) {
                     # Few distinct target values (e.g. ages 1yr, 2yr): jitter so points don't stack.
                     p <- p + suppressWarnings(geom_jitter(aes(text = hover), width = 0.02 * diff(range(df$x)), height = 0, size = 2.5, alpha = 0.8))
                   } else {
                     p <- p + suppressWarnings(geom_point(aes(text = hover), size = 2.5, alpha = 0.8))
                   }
                   subtitle <- if (is.null(res$strata_attr)) paste(res$method, r_lab("All")) else paste(res$method, "correlation per", res$strata_attr)
                   p + labs(x = res$target_axis, y = paste(gene_lab, "expression"), colour = colour_title,
                            title = paste(gene_lab, "vs", res$target_name), subtitle = subtitle) +
                     theme_bw() +
                     theme(axis.title = element_text(size = 16), axis.text = element_text(size = 14),
                           plot.title = element_text(size = 18, face = "bold"), plot.subtitle = element_text(size = 14),
                           legend.text = element_text(size = 12))
                 })

                 # Interactive view of the same ggplot; the static ggplot is what Save to output keeps.
                 output$target_plot <- plotly::renderPlotly({
                   p <- target_plot()
                   plotly::ggplotly(p, tooltip = "text") %>%
                     plotly::layout(title = list(text = paste0("<b>", p$labels$title, "</b><br><sup>", p$labels$subtitle, "</sup>"),
                                                 x = 0.02, xanchor = "left"),
                                    margin = list(t = 70))
                 })

                 # r matrix of the top genes (rows) by level (columns), shared by the
                 # interactive (morpheus) view and the static heatmap saved to output.
                 target_hm_matrix <- reactive({
                   res <- TargetResult()
                   validate(need(length(res$strata) > 1,
                                 "The Correlation Heatmap needs 'Compute separately for each level of' with at least 2 levels that could be computed."))
                   tbl <- res$table
                   # Follow the table as the user last sorted/searched it: rows_all is the
                   # row order (original indices) after sorting and filtering. Before the
                   # table has been shown, fall back to its default order (max_abs_r).
                   ord <- input$target_table_rows_all
                   ord <- ord[!is.na(ord) & ord >= 1 & ord <= nrow(tbl)]
                   if (length(ord) == 0) ord <- seq_len(nrow(tbl))
                   top <- tbl[ord[seq_len(min(input$target_hm_n, length(ord)))], , drop = FALSE]
                   m <- as.matrix(top[, paste0(res$strata, "_r"), drop = FALSE])
                   colnames(m) <- res$strata
                   rownames(m) <- make.unique(gene_display_label(top$UniqueID, res$PGN))
                   m
                 })

                 target_heatmap <- reactive({
                   res <- TargetResult()
                   m <- target_hm_matrix()
                   Heatmap(m, name = "r",
                           col = colorRamp2(c(-1, 0, 1), c("#2166AC", "white", "#B2182B")),
                           na_col = "grey85",
                           cluster_columns = FALSE,
                           cluster_rows = nrow(m) > 1 && !anyNA(m),
                           row_names_gp = grid::gpar(fontsize = 9),
                           column_title = paste(res$method, "r with", res$target_name, "per", res$strata_attr),
                           heatmap_legend_param = list(at = c(-1, -0.5, 0, 0.5, 1)))
                 })

                 output$target_heatmap <- morpheus::renderMorpheus({
                   m <- target_hm_matrix()
                   cluster <- nrow(m) > 1 && !anyNA(m)
                   w <- morpheus::morpheus(m,
                                           Rowv = cluster, Colv = FALSE,
                                           dendrogram = if (cluster) "row" else "none",
                                           colorScheme = list(scalingMode = "fixed", min = -1, max = 1,
                                                              fractions = c(0, 0.5, 1), values = c(-1, 0, 1),
                                                              colors = c("#2166AC", "#FFFFFF", "#B2182B")))
                   # Same as the Heatmap module: morpheus() leaves a global htmlwidgets
                   # option set that breaks DT tables rendered afterwards; reset it.
                   later::later(function() options(htmlwidgets.TOJSON_ARGS = NULL), delay = 0)
                   w
                 })

                 observeEvent(input$save_target_table, {
                   if (input$compute_target == 0) {
                     showNotification("Click Compute/Refresh first.", type = "warning")
                     return()
                   }
                   saved_table$Correlation_GeneTarget <- TargetResult()$table
                   log_info("Gene-Target correlation table saved to output.")
                   showNotification("Gene-Target table saved to output.", type = "message")
                 })

                 observeEvent(input$save_target_plot, {
                   if (input$compute_target == 0) {
                     showNotification("Click Compute/Refresh first.", type = "warning")
                     return()
                   }
                   res <- TargetResult()
                   id <- paste(gsub(" ", "_", gene_display_label(target_selected_uid(), res$PGN)),
                               gsub(" ", "_", res$target_name), sep = "-")
                   saved_plots$CorrTargetPlot[[id]] <- target_plot()
                   log_info("Gene-Target plot (", id, ") saved to output.")
                   showNotification("Plot saved to output.", type = "message")
                 })

                 observeEvent(input$save_target_heatmap, {
                   if (input$compute_target == 0) {
                     showNotification("Click Compute/Refresh first.", type = "warning")
                     return()
                   }
                   if (length(TargetResult()$strata) < 2) {
                     showNotification("The Correlation Heatmap needs 'Compute separately for each level of'.", type = "warning")
                     return()
                   }
                   n <- length(saved_plots$CorrTargetHeatmap)
                   saved_plots$CorrTargetHeatmap[[n + 1]] <- target_heatmap()
                   log_info("Gene-Target Correlation Heatmap saved to output.")
                   showNotification("Heatmap saved to output.", type = "message")
                 })

                 output$selectGroupSampleCorr<- renderUI(shared_header_content())
                 
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



