###########################################################################################################
## Proteomics Visualization R Shiny App
##
##This software belongs to Biogen Inc. All right reserved.
##
##@file: network.R
##@Developer : Lin Tinchi(tinchi.lin@biogen.com); Benbo Gao (benbo.gao@Biogen.com); Kyra Griffin-Mitchell (kyra.griffinmitchell@Biogen.com)
##@Date : 02/23/2022
##@version 1.0
###########################################################################################################

##########################################################################################################
## WGCNA
##########################################################################################################
library(WGCNA)
library(plotly)

clean_expression_data <- function(df) {
  df <- as.data.frame(df)
  numeric_cols <- sapply(df, is.numeric)
  df <- df[, numeric_cols, drop = FALSE]
  is_inf <- is.infinite(as.matrix(df))
  if (any(is_inf)) {
    message("Converting ", sum(is_inf), " infinite values to NA")
    df[is_inf] <- NA
  }
  keep <- complete.cases(df)
  df <- df[keep, , drop = FALSE]
  return(df)
}
get_wgcna_netwk <-function(dataExpr, picked_power, scenario_number, set_mergeCutHeight, set_maxBlockSize, ProjectID) {
  if ( scenario_number == 1 ) {
    set_loadTOM = TRUE
    TOMFileBase = paste0("./data/wgcna_data/TOM_", ProjectID)
  } else if ( scenario_number == 2 ) {
    set_loadTOM = FALSE
    TOMFileBase = paste0("./data/wgcna_data/TOM_", ProjectID)
  } else if ( scenario_number == 3 ) {
    set_loadTOM = FALSE
    TOMFileBase = "ER"
  }
  
  print(paste0("**** scenario ", as.character(scenario_number), " ****"))
  
  t3 <- Sys.time()
  temp_cor <- cor
  cor <- WGCNA::cor         # Force it to use WGCNA cor function (fix a namespace conflict issue)
  netwk <- blockwiseModules(dataExpr,                # <= input here
                            # == Adjacency Function ==
                            power = picked_power,                # <= power here
                            networkType = "signed",
                            # == Tree and Block Options ==
                            deepSplit = 2L,
                            pamRespectsDendro = F,
                            # detectCutHeight = 0.75,
                            minModuleSize = min(20, ncol(dataExpr/2)), # al# 30, #input$minModuleSize, #30,
                            # set block size to be number of genes, so that all
                            # genes will be analyzed in a single block
                            maxBlockSize = set_maxBlockSize,
                            # == Module Adjustments ==
                            reassignThreshold = 0,
                            mergeCutHeight = set_mergeCutHeight,#,0.25,
                            # == TOM == Archive the run results in TOM file (saves time)
                            saveTOMs = F,
                            loadTOM = set_loadTOM,
                            # Note: When launching from server, the path for TOM should be
                            # paste0("/mnt/depts/dept04/compbio/projects/xOmicsShiny/data/wgcna_data/TOM_",x)
                            saveTOMFileBase = TOMFileBase,
                            # == Output Options
                            numericLabels = T,
                            verbose = 3L)
  t4 <- Sys.time()
  cat(paste0("scenario ", as.character(scenario_number), " run WGCNA: ", round(difftime(t4, t3, units='mins'),2), " min\n"))
  cor <- temp_cor
  return(netwk)
}

get_wgcna_table <-function(wgcna, ME_name_updated, ProteinGeneName, gene_label) {
  module_map <- tibble::tibble(
    ME_name = ME_name_updated,
    color = sub("^ME\\d+_", "", ME_name_updated),
  )
  
  t2 <- tibble::tibble(
    UniqueID = names(wgcna$colors),
    color = labels2colors(wgcna$colors)
  ) %>%
    dplyr::left_join(ProteinGeneName[, c("UniqueID","Gene.Name")], by = "UniqueID") %>%
    dplyr::select(color, all_of(gene_label)) %>%
    dplyr::rename(gene = gene_label) %>%
    dplyr::group_by(color) %>%
    dplyr::summarise(
      n_gene = n(),
      gene_group = paste0(gene, collapse = ","),
      .groups = "drop"
    ) %>%
    dplyr::left_join(module_map, by = "color") %>%   
    dplyr::mutate(
      copy = purrr::map_chr(seq_len(n()), ~ paste0(
        as.character(
          rclipButton(
            paste0("clipbtn_", .x),
            label = "Copy all genes",
            clipText = gene_group[.x],
            icon = icon("copy", lib = "glyphicon"),
            class = "btn-primary btn-sm"
          )
        ),
        " ",
        sprintf(
          '<button id="runORA_%s" class="btn btn-info btn-sm" data-genes="%s">Run ORA</button>',
          .x,
          gene_group[.x]
        )
      ))
  )  %>%
  dplyr::select(color = ME_name, n_gene, copy, gene_group)
  return(t2)
}

plot_gene_variance_distribution <- function(gene_variance, top_N = NULL) {
  hist(gene_variance,
       breaks = 50,
       main = "Distribution of Gene Variance",
       xlab = "Normalized Variability Score",
       col = "lightblue",
       border = "white")
  # If top_N is provided, add cutoff line + legend
  if (!is.null(top_N)) {
    sorted_var <- sort(gene_variance, decreasing = TRUE)
    cutoff_value <- sorted_var[top_N]
    abline(v = cutoff_value,
           col = "red", lty = 2, lwd = 2)
    legend("topright",
           legend = paste0("Top ", top_N, " genes cutoff"),
           col = "red", lty = 2, lwd = 2,
           bty = "n")
  }
}

plot_sample_clustering <- function(dataExpr, cutoff = NULL) {
  dataExpr <- t(dataExpr)
  sample_tree <- hclust(dist(dataExpr), method = "average")
  
  plot(sample_tree,
       main = "Sample Clustering to Detect Outliers",
       sub= "",
       xlab = "",
       cex = 0.7)
  # Optional cutoff line
  if (!is.null(cutoff)) {
    abline(h = cutoff, col = "red", lty = 2, lwd = 2)
  }
}

plot_soft_threshold_diagnose <- function(sft, powers) {
  library(ggplot2)
  library(patchwork)
  
  df <- sft$fitIndices
  df$Power <- df[, 1]
  df$SFT_R2 <- -sign(df[, 3]) * df[, 2]
  df$MeanK <- df[, 5]
  
  # Panel 1: Scale-free topology fit
  p1 <- ggplot(df, aes(x = Power, y = SFT_R2, label = powers)) +
    geom_point(color = "red") +
    geom_text(vjust = -0.5, color = "red", size = 5) +
    geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
    labs(
      x = "Soft Threshold (power)",
      y = "Scale Free Topology Model Fit, signed R^2",
      title = "Scale Independence"
    ) +
    theme_minimal(base_size = 16)
  
  # Panel 2: Mean connectivity
  p2 <- ggplot(df, aes(x = Power, y = MeanK, label = powers)) +
    geom_point(color = "red") +
    geom_text(vjust = -0.5, color = "red", size = 5) +
    labs(
      x = "Soft Threshold (power)",
      y = "Mean Connectivity",
      title = "Mean Connectivity"
    ) +
    theme_minimal(base_size = 16)
  
  # Combine panels
  p1 + p2
}

rename_MEs_blockwise <- function(MEs, moduleColors) {
  colorNames <- labels2colors(moduleColors)
  ME_colors <- labels2colors(as.numeric(gsub("ME", "", names(MEs))))
  unique_colors <- unique(ME_colors)
  colorIDs <- match(ME_colors, unique_colors)
  newNames <- paste0("ME", colorIDs, "_", ME_colors)
  names(MEs) <- newNames
  return(MEs)
}

is_categorical <- function(x) {
  is.factor(x) || 
    is.character(x) 
}

build_traits_matrix <- function(MetaData, attrs, base_levels) {
  # Identify categorical traits (factor or character)
  categorical_cols <- names(base_levels)
  # Numeric traits = selected attrs minus categorical 
  numeric_cols <- setdiff(attrs, categorical_cols)
  # Extract numeric traits 
  numeric_traits <- MetaData[, numeric_cols, drop = FALSE]
  # Convert selected categorical columns to factors with correct base level
  factor_traits <- lapply(categorical_cols, function(col) {
    vals <- MetaData[[col]]
    base <- base_levels[[col]]
    factor(vals, levels = c(base, setdiff(unique(vals), base)))
  })
  factor_traits <- as.data.frame(factor_traits)
  names(factor_traits) <- categorical_cols
  if (ncol(factor_traits) > 0) {
    dummy_traits <- model.matrix(~ . - 1, data = factor_traits)
    dummy_traits <- as.data.frame(dummy_traits)
    return(cbind(numeric_traits, dummy_traits))
  } else {
    dummy_traits <- NULL
    return(numeric_traits)
  }
}

make_plotly_heatmap <- function(cor_mat, p_mat, row_labels, col_labels) {
  hover_text <- matrix("", nrow = nrow(cor_mat), ncol = ncol(cor_mat))
  for (i in seq_len(nrow(cor_mat))) {
    for (j in seq_len(ncol(cor_mat))) {
      hover_text[i, j] <- paste0(
        "Module: ", row_labels[i], "<br>",
        "Trait: ", col_labels[j], "<br>",
        "Correlation: ", signif(cor_mat[i, j], 3), "<br>",
        "P-value: ", signif(p_mat[i, j], 3)
      )
    }
  }
  plot_ly(
    x = col_labels,
    y = row_labels,
    z = cor_mat,
    type = "heatmap",
    colors = colorRamp(c("blue", "white", "red")),
    text = hover_text,
    hoverinfo = "text",
    source = "moduleTraitHeatmap"
  )
}

compute_module_trait_heatmap <- function(MEs, trait_matrix) {
  nSamples <- nrow(trait_matrix)
  
  cor_matrix <- cor(MEs, trait_matrix, use = "pairwise.complete.obs")
  p_mat  <- corPvalueStudent(cor_matrix, nSamples)
  
  make_plotly_heatmap(
    cor_mat = cor_matrix,
    p_mat = p_mat,
    row_labels = colnames(MEs),
    col_labels = colnames(trait_matrix)
  )
}

wgcna_ui <- function(id) {
  ns <- shiny::NS(id)
  fluidRow(
    rclipboard::rclipboardSetup(),
    column(3,
           wellPanel(
             uiOutput(ns('loadedprojects')),
             radioButtons(ns("WGCNAgenelable"),label="Select Gene Label",inline = TRUE, choices=c("Gene.Name","UniqueID"), selected="Gene.Name"),
             conditionalPanel(ns = ns, "input.WGCNA_tabset=='WGCNA Result' || input.WGCNA_tabset=='Module Eigengenes' || input.WGCNA_tabset=='WGCNA QC'",
                              numericInput(ns("WGCNAtopNum"), label= "Select Top N Genes, where N is :",  value=250L, min=250L, step=25L, max = 10000L),
                              numericInput(ns("mergeCutHeight"), label= "Dendrogram Cut Height for Merging:",  value=0.25, min= 0, max = 1.0, step = 0.01),
             ),
             conditionalPanel(ns = ns, "input.WGCNA_tabset=='WGCNA Result'",
                              actionButton(ns("plotwgcna"),"Run"),
                              br(),
                              strong("Running the data could take 3-10 minutes",style="color:red"), span("depending on data size and complexity; once Run starts, please refrain from clicking the button repeatedly.",style="color:red", inline = TRUE)
             ),
             conditionalPanel(ns = ns, "input.WGCNA_tabset=='Module-Trait Relationships'",
                              conditionalPanel(ns = ns, "input.Module_Trait=='Module Trait Heatmap'",
                                               selectizeInput(ns("WGCNA_trait_var"),label="Select attributes", choices = NULL,multiple = TRUE,options = list(placeholder = "Choose one or more attributes")),
                                               uiOutput(ns("attribute_settings_ui")),
                                               actionButton(ns("plot_module_trait"),"Plot Module Trait Correlation")
                                               ),
                              conditionalPanel(ns = ns, "input.Module_Trait=='Hub Gene Identification Table' || input.Module_Trait == 'Hub Gene Identification Plot'",
                                               selectInput(ns("WGCNA_trait"),label="Select a trait", choices = NULL ,multiple = FALSE),
                                               selectInput(ns("WGCNA_module"),label="Select a module (optional)", choices = NULL ,multiple = FALSE),
                                               actionButton(ns("plot_module_hub"),"Run")
                                               )
                              )
             )
    ),
    column(9,
           tabsetPanel(id=ns("WGCNA_tabset"),
                       tabPanel(title="WGCNA Result",
                                tabsetPanel(id=ns("WGCNA_Result"),
                                            tabPanel(title="Dendrogram Plot", 
                                                     tags$details(
                                                       tags$summary(
                                                         tags$strong("Click to view Dendrogram & Network Interpretation", 
                                                                     style = "color: #28a745; cursor: pointer;")
                                                       ),
                                                       tags$div(
                                                         style = "padding: 15px; background-color: #f4faf6; border-left: 5px solid #28a745; margin-top: 10px;",
                                                         tags$p("The ", tags$strong("Cluster Dendrogram"), " visualizes how genes are grouped based on their co-expression distance. Genes on lower branches are more tightly correlated."),
                                                         tags$p(tags$strong("How to Interpret:")),
                                                         tags$p(
                                                           "The ", tags$strong("Dendrogram Cut Height"), 
                                                           " (currently 0.25) determines the sensitivity of module merging. A higher cut height results in fewer, larger modules, while a lower height preserves smaller, more distinct gene sets."
                                                         ),
                                                         tags$p(
                                                           tags$strong("Eigengene Network:"), 
                                                           " The heatmap and tree show the relationships between the modules themselves. Modules that cluster together (represented by red blocks in the heatmap) show similar expression trends across your experimental conditions."
                                                         ),
                                                         tags$p(
                                                           tags$strong("Grey Module Note:"), 
                                                           " The 'Grey' module contains genes that did not fit into any distinct cluster. These are typically considered background noise and are usually excluded from biological interpretation."
                                                         )
                                                       )
                                                     ),
                                                     plotOutput(ns("Dendrogram"), height=800)),
                                            tabPanel(title="Gene Cluster Summary Table",
                                                     br(),
                                                     tags$details(
                                                       tags$summary(
                                                         tags$strong("Click to view Module Summary Table Guide", 
                                                                     style = "color: #28a745; cursor: pointer;")
                                                       ),
                                                       tags$div(
                                                         style = "padding: 15px; background-color: #f4faf6; border-left: 5px solid #28a745; margin-top: 10px;",
                                                         tags$p("This table summarizes the gene clusters (modules) identified by WGCNA. Each module represents a group of genes that share a highly similar expression profile across your samples."),
                                                         tags$p(tags$strong("Column Definitions:")),
                                                         tags$ul(
                                                           tags$li(tags$strong("Cluster:"), " The module name (e.g., ME1_turquoise). 'ME' stands for ", tags$em("Module Eigengene"), ", which is the mathematical representative of all genes in that cluster."),
                                                           tags$li(tags$strong("Number of Genes:"), " The total count of genes assigned to that module based on the Dendrogram cut."),
                                                           tags$li(tags$strong("Action:"), " Provides tools to export the gene list or run ", tags$strong("ORA (Over-Representation Analysis)"), " to identify biological pathways.")
                                                         ),
                                                         tags$p(tags$strong("Example Interpretation:")),
                                                         tags$p(
                                                           "If you observe that ", tags$strong("ME3_blue"), " has a high correlation with a specific condition in the 'Module-Trait Relationships' tab, you can use the ", tags$strong("Run ORA"), " button here to see if those 1,018 genes are involved in processes like 'Cell Cycle' or 'Immune Response'."
                                                         ),
                                                         tags$p(
                                                           tags$span(style = "color: #6c757d; font-style: italic;", 
                                                                     "Note: The 'Grey' module contains genes that could not be assigned to any cluster; these are typically excluded from downstream analysis.")
                                                         )
                                                       )
                                                     ),
                                                     br(),
                                                     actionButton(ns("wgcna_gct"), "Save GCT data file to output"),
                                                     br(),br(),
                                                     DT::dataTableOutput(ns("gene_cluster")),
                                                     tags$script(HTML(sprintf("
                                                     $(document).on('click', '[id^=runORA_]', function() {
                                                     var genes = $(this).attr('data-genes');
                                                     Shiny.setInputValue('%s', genes, {priority: 'event'});
                                                     });", ns("runORA_trigger"))))
                                            )
                                )
                       ),
                       tabPanel(title="Module Eigengenes",
                                tabsetPanel(id=ns("Module_Eigengenes"),
                                            tabPanel(title="Eigengene table", 
                                                     br(),
                                                     tags$details(
                                                       tags$summary(
                                                         tags$strong("Click to view Module Eigengene Table Guide", 
                                                                     style = "color: #28a745; cursor: pointer;")
                                                       ),
                                                       tags$div(
                                                         style = "padding: 15px; background-color: #f4faf6; border-left: 5px solid #28a745; margin-top: 10px;",
                                                         tags$p("This table displays the ", tags$strong("Module Eigengene (ME)"), " values for each individual sample. These values represent the relative expression level of an entire module within a specific sample."),
                                                         tags$p(
                                                           "Each row corresponds to a sample (e.g., ", tags$code("2mo-WT-10F"), ") and each column represents a module identified in the co-expression network. ",
                                                           tags$strong("Positive values"), " indicate that genes in that module are upregulated in that sample relative to the average, while ", 
                                                           tags$strong("negative values"), " indicate downregulation."
                                                         ),
                                                         tags$p(tags$strong("How to use this data:")),
                                                         tags$ul(
                                                           tags$li(tags$strong("Pattern Recognition:"), " Look for columns where values change consistently between groups (e.g., all WT samples are negative and all Het samples are positive)."),
                                                           tags$li(tags$strong("Statistical Input:"), " These ME values are used as 'synthetic genes' for the Module-Trait Relationship analysis to correlate expression with clinical metadata."),
                                                           tags$li(tags$strong("Outlier Detection:"), " Samples with ME values that differ drastically from others in the same group may indicate biological outliers.")
                                                         )
                                                       )
                                                     ),
                                                     br(),
                                                     actionButton(ns("Eigengene"), "Save to output"),
                                                     br(),
                                                     DT::dataTableOutput(ns("MEs"))),
                                            tabPanel(title = "Eigengene Network",
                                                     br(),
                                                     tags$details(
                                                       tags$summary(
                                                         tags$strong("Click to view Eigengene Network Guide", 
                                                                     style = "color: #28a745; cursor: pointer;")
                                                       ),
                                                       tags$div(
                                                         style = "padding: 15px; background-color: #f4faf6; border-left: 5px solid #28a745; margin-top: 10px;",
                                                         tags$p("The Eigengene Network visualizes how different gene modules relate to one another. This helps identify 'meta-modules'—groups of clusters that follow similar biological trends."),
                                                         
                                                         tags$p(tags$strong("1. Eigengene Dendrogram (Top):")),
                                                         tags$ul(
                                                           tags$li("This tree clusters modules (not individual genes) based on their expression similarity."),
                                                           tags$li("Modules that branch together (e.g., ", tags$strong("ME5_blue"), " and ", tags$strong("ME6_brown"), ") are highly correlated and may be part of the same broad biological response.")
                                                         ),
                                                         
                                                         tags$p(tags$strong("2. Adjacency Heatmap (Bottom):")),
                                                         tags$ul(
                                                           tags$li(tags$strong("Red Squares:"), " Indicate high positive correlation. If two modules are dark red, they increase and decrease in expression together across your samples."),
                                                           tags$li(tags$strong("Blue Squares:"), " Indicate negative correlation. One module is 'turned on' while the other is 'turned off'."),
                                                           tags$li(tags$strong("Diagonal Line:"), " The dark red diagonal represents a module's correlation with itself (always 1.0).")
                                                         ),
                                                         
                                                         tags$p(tags$strong("Why this is useful:")),
                                                         tags$p(
                                                           "If multiple modules show nearly identical patterns in this network, you might consider increasing the ", 
                                                           tags$strong("Dendrogram Cut Height"), 
                                                           " in the sidebar to merge them into a single, more robust biological signature."
                                                         )
                                                       )
                                                     ),
                                                     br(),
                                                     plotOutput(ns("Eigenene_Network"), height = "1000px"))
                                )
                       ),
                       tabPanel(title="Module-Trait Relationships", 
                                tabsetPanel(id=ns("Module_Trait"),
                                            tabPanel(title="Module Trait Heatmap", 
                                                     br(),
                                                     tags$details(
                                                       tags$summary(
                                                         tags$strong("Click to view Module-Trait Relationship Guide", 
                                                                     style = "color: #28a745; cursor: pointer;")
                                                       ),
                                                       tags$div(
                                                         style = "padding: 15px; background-color: #f4faf6; border-left: 5px solid #28a745; margin-top: 10px;",
                                                         tags$p("This heatmap illustrates the correlation between ", tags$strong("Module Eigengenes"), " and your experimental design factors (Traits). It is the key step in identifying which gene clusters are driven by your biological conditions."),
                                                         tags$p(tags$strong("How to read the Heatmap:")),
                                                         tags$ul(
                                                           tags$p("Mouse over a cell to show correlation results."),
                                                           tags$li(tags$strong("Color Scale:"), " Red indicates a positive correlation (the module is upregulated in that condition), while Blue indicates a negative correlation (downregulated)."),
                                                           tags$li(tags$strong("Cell Values:"), " The top number in each cell is the Correlation Coefficient, and the bottom number (in parentheses) is the ", tags$strong("p-value"), "."),
                                                           tags$li(tags$strong("Significance:"), " Focus on cells where the p-value is < 0.05. These modules are statistically associated with that specific trait.")
                                                         ),
                                                         tags$p(tags$strong("Example Interpretation:")),
                                                         tags$p(
                                                           "If the column for ", tags$strong("group2mo_WT"), " shows a high positive correlation (Red) with ", tags$strong("ME4_red"), 
                                                           ", it suggests that the genes within the blue module are highly active specifically in the wT samples."
                                                         ),
                                                         tags$p(
                                                           "Click on a cell will launch the selected modlue and trait to ", tags$strong("Hub Gene Identification Table"), " Tab."
                                                         )
                                                       )
                                                     ),
                                                     br(),
                                                     plotlyOutput(ns("module_trait_hmap"), height=800)
                                            ),
                                            tabPanel(title="Hub Gene Identification Table",
                                                     br(),
                                                     tags$details(
                                                       tags$summary(
                                                         tags$strong("Click to view Hub Gene Identification Guide", 
                                                                     style = "color: #28a745; cursor: pointer;")
                                                       ),
                                                       tags$div(
                                                         style = "padding: 15px; background-color: #f4faf6; border-left: 5px solid #28a745; margin-top: 10px;",
                                                         tags$p("This table identifies 'Hub Genes'—highly connected central nodes within a module that are often the most biologically significant drivers of a phenotype."),
                                                         tags$p(tags$strong("Metric Definitions:")),
                                                         tags$ul(
                                                           tags$li(tags$strong("Connectivity:"), " Measures how strongly a gene is connected to all other genes in the network. Higher values indicate a more central role."),
                                                           tags$li(tags$strong("ModuleMembership (kME):"), " The correlation between an individual gene's expression and the Module Eigengene. A value near 1.0 means the gene is a perfect representative of that module's trend."),
                                                           tags$li(tags$strong("GeneSignificance (GS):"), " The absolute correlation between the gene and the selected Trait (e.g., Genotype). High GS indicates the gene is strongly associated with the biological condition of interest.")
                                                         ),
                                                         tags$p(tags$strong("What makes a 'Hub Gene'?") ),
                                                         tags$p(
                                                           "Look for genes with both ", tags$strong("high ModuleMembership"), " and ", tags$strong("high GeneSignificance")
                                                         ),
                                                         tags$p(
                                                           "You can filter the columns to find genes with the highest connectivity to prioritize candidates for knock-down or validation experiments."
                                                         )
                                                       )
                                                     ),
                                                     br(),
                                                     actionButton(ns("hub_gene"), "Save to output"),
                                                     br(),
                                                     DT::dataTableOutput(ns("hub_gene_table"))
                                            ),
                                            tabPanel(
                                              title = "Hub Gene Identification Plot",
                                              br(),
                                              tags$details(
                                                tags$summary(
                                                  tags$strong("Click to view Hub Gene Plot Interpretation", 
                                                              style = "color: #28a745; cursor: pointer;")
                                                ),
                                                tags$div(
                                                  style = "padding: 15px; background-color: #f4faf6; border-left: 5px solid #28a745; margin-top: 10px;",
                                                  tags$p("These plots provide a visual validation of the 'Hub' status of genes within a specific module by comparing their internal network importance against their biological significance."),
                                                  tags$p(tags$strong("1. MM vs. GS Scatter Plot (Left):")),
                                                  tags$ul(
                                                    tags$li("This plot correlates ", tags$strong("Module Membership (MM)"), " with ", tags$strong("Gene Significance (GS)"), "."),
                                                    tags$li("A high positive correlation (e.g., cor > 0.8) indicates that genes most central to the module are also the ones most strongly associated with your trait (e.g., Genotype or Age)."),
                                                    tags$li("Genes in the top-right corner are your primary candidates for further biological study.")
                                                  ),
                                                  tags$p(tags$strong("2. Hub Gene Identification Plot (Right):")),
                                                  tags$ul(
                                                    tags$li("This plot ranks genes by their ", tags$strong("Intramodular Connectivity"), "."),
                                                    tags$li("The ", tags$span(style = "color: red; font-weight: bold;", "Red Dots"), " highlight the top 10% most connected genes. These 'Hubs' are theoretically the most influential regulators within this specific module.")
                                                  ),
                                                  tags$p(tags$strong("Usage Tip:")),
                                                  tags$p(
                                                    "Select different modules from the sidebar to see if the correlation holds true across different clusters. Strong MM-GS correlations suggest that the module is a key driver of the biological differences observed between your groups."
                                                  )
                                                )
                                              ),
                                              br(),br(),
                                              div(
                                                style = "display: flex; gap: 20px;",
                                                div(
                                                  style = "flex: 1;",
                                                  plotlyOutput(ns("plot_MMvsGS"), height = "700px")
                                                ),
                                                div(
                                                  style = "flex: 1;",
                                                  plotlyOutput(ns("plot_MMvsConnectivity"), height = "700px")
                                                )
                                              )
                                              
                                            )
                                )
                       ),
                       tabPanel(title="WGCNA QC",
                                tabsetPanel(id=ns("WGCNA_QC"),
                                            tabPanel(title="Soft-Thresholding Power", 
                                                     br(),
                                                     # tags$style("#soft_threshold_table_note { font-size: 18px; color: steelblue; }"),
                                                     # htmlOutput(ns("soft_threshold_table_note")),
                                                     tags$details(
                                                       tags$summary(
                                                         tags$strong("Click to view Soft-Thresholding Power Guide", 
                                                                     style = "color: #28a745; cursor: pointer;")
                                                       ),
                                                       tags$div(
                                                         style = "padding: 15px; background-color: #f4faf6; border-left: 5px solid #28a745; margin-top: 10px;",
                                                         tags$p("The soft-thresholding power (β) is a critical parameter that emphasizes strong correlations and suppresses noise to achieve a 'scale-free' network topology. The table below showed the result for all β candidates tested. ", tags$strong("The row is highlighted for the picked power β.")),
                                                         
                                                         tags$p(tags$strong("How to Choose the Optimal Power:")),
                                                         tags$ul(
                                                           tags$p("In practice, soft‑thresholding power (β) is chosen with a balanced good scale‑free fit and reasonable connectivity. WGCNA recommends choosing the lowest power where R² ≥ 0.8 (or 0.9 for stricter networks). And typical good range of mean connectivity is 20–200, depending on dataset size."),
                                                           tags$li(tags$strong("The R² Criterion:"), " Aim for the lowest power where the ", tags$code("SFT.R.sq"), " reaches a plateau, ideally above 0.80 or 0.90."),
                                                           tags$li(tags$strong("Mean Connectivity:"), " Ensure the connectivity is not too low. A typical range is between 20 and 200 for biological relevance."),
                                                           tags$li(tags$strong("Highlighted Row:"), " The blue highlighted row in the table indicates the power currently used by the system (e.g., Power 6 in this dataset).")
                                                         ),
                                                         
                                                         tags$p(tags$strong("Practical Decision Rules:")),
                                                         tags$ul(
                                                           tags$li("If R² rises and then plateaus, pick the lowest power at the start of the plateau."),
                                                           tags$li("If R² rises very slowly, a power between 6 and 10 is generally acceptable for signed networks."),
                                                           tags$li("Avoid very high powers (>20) as they can cause the network to become too sparse, leading to a collapse in connectivity."),
                                                           tags$li("If estimated power is NA, use any number in the typical β range 6 to 12 for signed networks.")
                                                         ),
                                                         
                                                         tags$p(tags$strong("Key Plots Interpretation:")),
                                                         tags$ul(
                                                           tags$li(
                                                             tags$strong("Scale-Free Topology Fit:"),
                                                             " Biological networks are typically 'scale-free.' This plot helps you choose the ", tags$em("Soft Threshold"), 
                                                             " (Power). You are looking for the lowest power where the curve levels off (usually reaching a R^2 fit of 0.8 or 0.9)."
                                                           ),
                                                           tags$li(
                                                             tags$strong("Mean Connectivity:"), 
                                                             " As the power increases, the average connectivity decreases. This plot ensures the network remains connected enough to be biologically meaningful without being overly noisy."
                                                           )
                                                         ),
                                                         
                                                         tags$p(tags$strong("Important Note:")),
                                                         tags$p(
                                                           "If the Scale-Free Topology plot does not reach a plateau above 0.8, it may indicate that the data has a strong batch effect or that the samples do not follow a typical co-expression structure."
                                                         )
                                                       )
                                                     ),
                                                     br(),
                                                     DT::dataTableOutput(ns("soft_threshold_table")),
                                                     br(),
                                                     # tags$style("#soft_threshold_diagnose_plot_note { font-size: 18px; color: steelblue; }"),
                                                     # htmlOutput(ns("soft_threshold_diagnose_plot_note")),
                                                     # br(),
                                                     plotOutput(ns("soft_threshold_diagnose_plot"), height=800)
                                            ),
                                            tabPanel(title="WGCNA Network",
                                                     br(),
                                                     tags$details(
                                                       tags$summary(
                                                         tags$strong("Click to view Data Variance & Sample Clustering Guide", 
                                                                     style = "color: #28a745; cursor: pointer;")
                                                       ),
                                                       tags$div(
                                                         style = "padding: 15px; background-color: #f4faf6; border-left: 5px solid #28a745; margin-top: 10px;",
                                                         tags$p("Before building a network, it is essential to assess the global structure of the dataset. These plots help confirm that the most variable genes are being captured and that the samples group logically."),
                                                         
                                                         tags$p(tags$strong("1. Distribution of Gene Variance:")),
                                                         tags$ul(
                                                           tags$li("This histogram shows the variability of gene expression across all samples."),
                                                           tags$li("WGCNA typically focuses on the top most variable genes (e.g., top 5,000 or 10,000) to reduce noise and computational load. The ", tags$strong("Select Top N Genes"), " input in the sidebar allows you to adjust this threshold.")
                                                         ),
                                                         
                                                         tags$p(tags$strong("2. Sample Clustering to Detect Outliers:")),
                                                         tags$ul(
                                                           tags$li("This dendrogram groups samples based on their global expression profiles."),
                                                           tags$li(tags$strong("Biological Sanity Check:"), " Ideally, biological replicates (e.g., all WT samples) should cluster together."),
                                                           tags$li(tags$strong("Outlier Identification:"), " If a single sample is separated by a very long vertical branch from the rest of the data, it may be an outlier that should be excluded to improve module detection.")
                                                         ),
                                                         
                                                         tags$p(tags$strong("Pro Tip:")),
                                                         tags$p(
                                                           "If you see a sample that is clearly an outlier, you may need to filter it out in the ", tags$strong("Dataset"), " tab and re-run the WGCNA analysis to achieve a cleaner network structure."
                                                         )
                                                       )
                                                     ),
                                                     br(),
                                                     plotOutput(ns("gene_variance_distribution_plot"), height=400),
                                                     plotOutput(ns("sample_cluster_plot"), height=800)
                                            )
                                )
                       ),
                       tabPanel(title="Help", htmlOutput('help_WGCNA'))
           )
    )
  )
}

wgcna_server <- function(id, parent_session) {
  shiny::moduleServer(id,
                      function(input, output, session) {
                        ns <- shiny::NS(id)
                        working_project <- reactiveVal()
                        MEs_updated <- reactiveVal()  # updated MEs data.frame
                        MEs_name_updated <- reactiveVal()  # updated MEs data.frame colnames
                        df_gene_clusters <- reactiveVal()
                        powers <- reactiveVal()
                        trait_data <- reactiveVal() 
                        moduleTraitCor <- reactiveVal()
                        df_hub_gene <- reactiveVal() 
                        module_trait_plot <- reactiveVal(NULL)

                        observe({
                          req(ProjectInfo)
                          working_project(ProjectInfo$ProjectID)
                        })

                        observe({
                          DataIn <- DataReactive()
                          req(DataIn$data_wide)
                          data_wide <- na.omit(DataIn$data_wide)
                          default_n_gene <- min(10000, nrow(data_wide))
                          updateNumericInput(session, "WGCNAtopNum", 
                                             label= "Select Top N Genes, where N is :",  value=default_n_gene, min=250L, step=25L, max = default_n_gene)
                          working_project(ProjectInfo$ProjectID)
                        })
                        
                        # reset all results for new project
                        observeEvent(working_project(), {
                          output$Dendrogram <- NULL
                          output$gene_cluster <- NULL
                          output$MEs <- NULL
                          output$Eigenene_Network <- NULL
                          output$soft_threshold_table_note <- NULL
                          output$soft_threshold_table <- NULL
                          output$soft_threshold_diagnose_plot_note <- NULL
                          output$soft_threshold_diagnose_plot <- NULL
                          output$gene_variance_distribution_plot <- NULL
                          output$sample_cluster_plot <- NULL
                          module_trait_plot(NULL)
                          output$hub_gene_table <- DT::renderDT({NULL})
                          output$plot_MMvsGS <- renderPlotly({NULL})
                          output$plot_MMvsConnectivity <- renderPlotly({NULL})
                        })
                        
                        observeEvent(input$plotwgcna, {
                          # output$MEs <- NULL
                          output$Eigenene_Network <- NULL
                          output$soft_threshold_table_note <- NULL
                          output$soft_threshold_table <- NULL
                          output$soft_threshold_diagnose_plot_note <- NULL
                          output$soft_threshold_diagnose_plot <- NULL
                          output$gene_variance_distribution_plot <- NULL
                          output$sample_cluster_plot <- NULL
                          module_trait_plot(NULL)
                          output$hub_gene_table <- DT::renderDT({NULL})
                          output$plot_MMvsGS <- renderPlotly({NULL})
                          output$plot_MMvsConnectivity <- renderPlotly({NULL})
                        })
                        
                        observeEvent(list(input$WGCNA_trait, input$WGCNA_module), {
                          output$hub_gene_table <- DT::renderDT({NULL})
                          output$plot_MMvsGS <- renderPlotly({NULL})
                          output$plot_MMvsConnectivity <- renderPlotly({NULL})
                        })
                        
                        # Display precomputed results if exists.
                        observeEvent(list(working_project(),input$WGCNAgenelable,parent_session$input$menu == "wgcna"), {
                          req(ProjectInfo, DataReactive(), parent_session$input$menu == "wgcna")
                          ProjectID <- ProjectInfo$ProjectID
                          wgcnafile <- ProjectInfo$file3
                          load_wgcna_file <- file.path(dirname(wgcnafile), sub("^wgcna", "load", basename(wgcnafile)))
                          
                          gene_label <- input$WGCNAgenelable
                          DataIn = DataReactive()
                          req(DataIn$ProteinGeneName)
                          ProteinGeneName  <- DataIn$ProteinGeneName
                          
                          wgcna_out <- tryCatch(WGCNAReactive(), error = function(e) NULL)
                          if (! is.null(wgcna_out)) {
                            wgcna <- wgcna_out$netwk
                            if ("sft" %in% names(wgcna_out)) {
                              sft <- wgcna_out$sft
                            }
                          } else if(file.exists(wgcnafile) & file.exists(load_wgcna_file)) {
                            load(wgcnafile)
                            load(load_wgcna_file)
                            wgcna <- netwk
                          }
                          
                          if(exists('wgcna')){
                            mergedColors = labels2colors(wgcna$colors)
                            
                            MEs <- wgcna$MEs
                            moduleColors <- wgcna$colors
                            MEs_updated(rename_MEs_blockwise(MEs, moduleColors))
                            MEs_name_updated(colnames(MEs_updated()))
                            
                            # showTab(session = session, inputId = "WGCNA_tabset", target = "Module-Trait Relationships")
                            output$Dendrogram <- renderPlot({
                              withProgress(message = "Creating plot using pre-calculated data", value = 0, {
                                plotDendroAndColors(
                                  wgcna$dendrograms[[1]],
                                  mergedColors[wgcna$blockGenes[[1]]],
                                  "Module colors",
                                  dendroLabels = FALSE,
                                  hang = 0.03,
                                  addGuide = TRUE,
                                  guideHang = 0.05 )
                              })
                            })

                            # t0: merge WGCNA output with ProteinGeneName so that genes can be shown as UniqueID or Gene.Name
                            t2 <- get_wgcna_table(wgcna, MEs_name_updated(), ProteinGeneName, gene_label)
                            df_gene_clusters(t2)
                            
                            output$gene_cluster <- DT::renderDT({
                              DT::datatable(t2,
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
                            
                            
                            data_wide <- DataIn$data_wide
                            data_wide <- na.omit(DataIn$data_wide)
                            diff <- apply(data_wide, 1, sd, na.rm = TRUE)/(rowMeans(data_wide) + median(rowMeans(data_wide)))
                            
                            output$gene_variance_distribution_plot <- renderPlot({
                              plot_gene_variance_distribution(diff)
                            })
                            
                            output$sample_cluster_plot <- renderPlot({
                              plot_sample_clustering(data_wide)
                            })
                            output$MEs <- DT::renderDT({
                              DT::datatable(
                                MEs_updated(),  extensions = 'Buttons', escape = FALSE, selection = 'none', class = 'cell-border strip hover',
                                options = list(    dom = 'lBfrtip', pageLength = 15,
                                                   buttons = list(
                                                     list(extend = "csv", text = "Download Page", filename = "Page_results",
                                                          exportOptions = list(modifier = list(page = "current"))),
                                                     list(extend = "csv", text = "Download All", filename = "All_Results",
                                                          exportOptions = list(modifier = list(page = "all")))
                                                   )
                                )) %>% 
                                formatSignif(columns=names(MEs_updated()), digits=3)
                            })
                            
                            output$Eigenene_Network <- renderPlot({
                              MEs <- MEs_updated()
                              validate(need(ncol(MEs) > 2,"Eigengene Network requires at least 2 module eigengenes."))
                              
                              plotEigengeneNetworks(MEs, "Eigengene Network", 
                                                    marDendro = c(2,3,2,1),
                                                    marHeatmap =c(6,8,2,1),
                                                    plotDendrograms = TRUE, 
                                                    plotHeatmaps = TRUE)
                            })
##############                            
                            wgcna_gct <- reactive({
                              DataIn <- DataReactive()
                              req(DataIn$data_wide)
                              data_wide <- na.omit(DataIn$data_wide)
                              ProteinGeneName <- DataIn$ProteinGeneName
                              MetaData <- DataIn$MetaData
                              
                              module_map <- tibble::tibble(
                                ME_name = MEs_name_updated(),
                                color = sub("^ME\\d+_", "", MEs_name_updated()),
                              )
                              
                              t2 <- tibble::tibble(
                                UniqueID = names(wgcna$colors),
                                color = labels2colors(wgcna$colors)
                              ) %>% 
                                dplyr::left_join(module_map, by = "color") %>%
                                dplyr::inner_join(ProteinGeneName, by = "UniqueID") %>%
                                dplyr::select(-color)
                              
                              data.in <- DataIn$data_wide[t2$UniqueID,]
                              
                              gct_data <- create_gct_object(data.in, t2, MetaData)
                              gct_data
                            })
                            
                            observeEvent(input$wgcna_gct, {
                              saved_gcts$wgcna_gct <- wgcna_gct()
                            })                               
##############                            
                            if (exists("sft")) {
                              powers(sft$fitIndices$Power)
                              
                              output$soft_threshold_table_note <- renderText({
                                HTML("In WGCNA, soft‑thresholding power (β) controls how strongly you emphasize large correlations and suppress small ones when building the adjacency matrix. The table below showed the result for all β candidates tested. <b>The row is highlighted for the picked power β.<b>")
                              })
                              
                              output$soft_threshold_table <- DT::renderDT({
                                datatable(
                                  sft$fitIndices,
                                  rownames = FALSE, extensions = 'Buttons', escape = FALSE, selection = 'none', class = 'cell-border strip hover',
                                  options = list(dom = 'lBfrtip', pageLength = 15,
                                                 buttons = list(
                                                   list(extend = "csv", text = "Download Page", filename = "Page_results", 
                                                        exportOptions = list(modifier = list(page = "current"))),
                                                   list(extend = "csv", text = "Download All", filename = "All_Results",
                                                        exportOptions = list(modifier = list(page = "all")))
                                                   )
                                                 )
                                  ) %>% 
                                  formatSignif(
                                    columns = setdiff(names(Filter(is.numeric, sft$fitIndices)), "Power"),
                                    digits = 3
                                  ) %>% formatStyle(columns="Power", target="row", backgroundColor = styleEqual(picked_power, "skyblue", default='white'))
                              })
                              
                              output$soft_threshold_diagnose_plot_note <- renderUI({
                                HTML("
                                In practice, soft‑thresholding power (β) is chosen with a balanced good scale‑free fit and reasonable connectivity.
                                WGCNA recommends choosing the lowest power where R² ≥ 0.8 (or 0.9 for stricter networks). 
                                And typical good range of mean connectivity is 20–200, depending on dataset size.
                                <br><br>
                                
                                <b>Practical decision rules:</b><br>
                                -- If R² rises and then plateaus, pick the lowest β at the plateau.<br>
                                -- If R² rises slowly but connectivity is still good, pick β = 6–10 (signed).<br>
                                -- If R² is flat and connectivity collapses at high β, pick a moderate β (6–8) to avoid oversparsifying.<br>
                                -- If R² keeps rising up to the max tested power, increase the tested range (e.g., 1–30), but still avoid β > 20 unless absolutely necessary.<br>
                                -- Ensure mean connectivity is not too low. It should > ~ 20 (for typical datasets).<br>
                                -- If estimated power is NA, use any number in the typical β range (6–12) for signed networks.
                                     ")
                              })
                              
                              output$soft_threshold_diagnose_plot <- renderPlot({
                                plot_soft_threshold_diagnose(sft, powers())
                              })
                            } else if (input$WGCNA_tabset == 'WGCNA_QC' && input$WGCNA_QC == "Soft-Thresholding Power") {
                              showNotification("Pre-calculated wgcna file does not contain the soft_thresholding power testing result. No results loaded.", duration = 5, type = "warning")
                            }
                          } else if (parent_session$input$menu == "wgcna") {
                            # hideTab(session = session, inputId = "WGCNA_tabset", target = "Module-Trait Relationships")
                            print("no pre-computed wgcna file available and cannot load wgcna results")
                            showNotification("Cannot find pre-calculated wgcna file, no WGCNA results loaded.", duration = 5, type = "warning")
                          }
                        })
                        
                        # Run WGCNA, get dataExpr, network, picked_power, sft(if exists or computed) 
                        # use eventReactive to control reactivity of WGCNAReactive;
                        # otherwise, whenever an input change, WGCNAReactive will be re-calculated
                        # and its re-calculation could take a long time.
                        WGCNAReactive <- eventReactive(input$plotwgcna, {
                          withProgress(message = "Running WGCNA", detail = 'This may take a while...', value = 0.2, {
                            # what if the user-imported data doesn't have $data_wide, $ProjectID..etc?
                            req(ProjectInfo, DataReactive()$data_wide, ProjectInfo$ProjectID)
                            # showTab(session = session, inputId = "WGCNA_tabset", target = "Module-Trait Relationships")
                            
                            DataIn = DataReactive()
                            data_wide <- na.omit(DataIn$data_wide)
                            ProjectID <- ProjectInfo$ProjectID
                            wgcnafile <- ProjectInfo$file3
                            
                            diff <- apply(data_wide, 1, sd, na.rm = TRUE)/(rowMeans(data_wide) + median(rowMeans(data_wide)))
                            data_wide=data_wide[order(diff, decreasing=TRUE), ]                            
                            
                            if (nrow(data_wide)>10000 ) {
                              data_wide=data_wide[1:10000, ] 
                              cat("reduce gene size to 10K for project ", ProjectID, "\n")
                            } 

                            print(paste0("**** dim of dataExpr after-preprocssing is ****", dim(data_wide)))
                            
                            default_n_gene <- nrow(data_wide)
                            
                            load_wgcna_file <- file.path(dirname(wgcnafile), sub("^wgcna", "load", basename(wgcnafile)))
                            
                            if (file.exists(load_wgcna_file) & default_n_gene==input$WGCNAtopNum){
                              # Scenario 1: If file exist and the number of genes selected rename the same, load 
                              # pre-computed result and TOM file (blockwiseModules(loadTom = T)) to 
                              # reduce running time
                              # The load_*.RData contains two objects, dataExpr and picked_power, so that
                              # the app doesn't need to recalculate either from scratch
                              load(load_wgcna_file)
                              netwk <- get_wgcna_netwk(dataExpr, picked_power, 1, input$mergeCutHeight, input$WGCNAtopNum, ProjectID)
                            } else if (file.exists(load_wgcna_file) & (default_n_gene - input$WGCNAtopNum)/default_n_gene < 0.1) {
                              # Scenario 2: If file exist and the number of genes selected is within 10% of 
                              # the default number of genes, load pre-computed result
                              # but do not load TOM file (blockwiseModules(loadTom = F))
                              load(load_wgcna_file)
                              dataExpr= dataExpr[,1L:input$WGCNAtopNum]
                              netwk <- get_wgcna_netwk(dataExpr, picked_power, 2, input$mergeCutHeight, input$WGCNAtopNum, ProjectID)
                            } else {
                              # Scenario 3: Not scenario 1 or 2, and recalculate everything
                              print(paste0("**** compute everything from scratch ****"))
                              ProteinGeneName  <- DataIn$ProteinGeneName
                              topNum <- as.numeric(input$WGCNAtopNum)
                              gene_label <- input$WGCNAgenelable

                              data_wide <- clean_expression_data(data_wide)
                              dataExpr = as.data.frame(t(data_wide))
                              dataExpr= dataExpr[,1L:topNum]
                              
                              # Choose a set of soft-thresholding powers
                              powers(c(1L:10L, seq(from = 12L, to = 20L, by = 2L)))
                              t2 <- Sys.time()
                              cor <- stats::cor
                              sft <- WGCNA::pickSoftThreshold(dataExpr, dataIsExpr = TRUE, powerVector = powers(), corFnc = cor, corOptions = list(use = 'p'), networkType = "signed")
                              
                              # Generating adjacency and TOM similarity matrices based on the selected softpower
                              if (!is.na(sft$powerEstimate)){
                                print("**** Pick power from sft$powerEstmate **** ")
                                picked_power <- softPower <- sft$powerEstimate
                              } else {
                                print("**** Use 6 as default if automatic selection fails **** ")
                                picked_power <- 6L
                              }
                              
                              t3 <- Sys.time()
                              cat(paste0("scenario 3 computing softpower: ", round(difftime(t3, t2, units='mins'),2), " min\n"))
                              netwk <- get_wgcna_netwk(dataExpr, picked_power, 3, input$mergeCutHeight, input$WGCNAtopNum, ProjectID)
                            }
                            out <- list(netwk = netwk, picked_power = picked_power, dataExpr = dataExpr)
                            if (exists("sft")) {
                              out$sft <- sft
                            }
                            return(out)
                          })
                        })
                        
                        observeEvent(input$Eigengene, {
                          saved_table$Eigengene <- MEs_updated()
                        })
                        
                        # Show WGCNA QC and result tables and plots #####
                        # use input$WGCNAReactive() as event handler to ensure observeEvent() depends on it only
                        # and does not directly depends on input$, which ensure WGCNAReactive() will be calculated first.
                        observeEvent(input$plotwgcna,{
                          req(DataReactive(),WGCNAReactive())
                          wgcna_out <- WGCNAReactive()
                          wgcna <- wgcna_out$netwk
                          picked_power <- wgcna_out$picked_power
                          DataIn = DataReactive()
                          mergedColors = labels2colors(wgcna$colors)
                          
                          MEs <- wgcna$MEs
                          moduleColors <- wgcna$colors # here are numbers
                          MEs_updated(rename_MEs_blockwise(MEs, moduleColors))
                          MEs_name_updated(colnames(MEs_updated()))

                          output$Dendrogram <- renderPlot({
                            plotDendroAndColors(
                              wgcna$dendrograms[[1]],
                              mergedColors[wgcna$blockGenes[[1]]],
                              "Module colors",
                              dendroLabels = FALSE,
                              hang = 0.03,
                              addGuide = TRUE,
                              guideHang = 0.05 )
                          })
                          
                          # generate table showing clustered genes #
                          ProteinGeneName  <- DataIn$ProteinGeneName
                          gene_label <- input$WGCNAgenelable
                          
                          t2 <- get_wgcna_table(wgcna, MEs_name_updated(), ProteinGeneName, gene_label)
                          df_gene_clusters(t2)

                          output$gene_cluster <- DT::renderDT({
                            DT::datatable(t2,
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
                          
                          output$MEs <- DT::renderDT({
                            DT::datatable(MEs_updated(),  
                                          extensions = 'Buttons', 
                                          escape = FALSE, 
                                          selection = 'none', 
                                          class = 'cell-border strip hover',
                                          options = list(dom = "Blfrtip", 
                                                         buttons = c("csv", "excel", "print")
                                                         )
                                          ) %>% 
                              formatSignif(columns=names(MEs_updated()), digits=3)
                          })
                          
                          output$Eigenene_Network <- renderPlot({
                            MEs <- MEs_updated()
                            validate(need(ncol(MEs) > 2,"Eigengene Network requires at least 2 module eigengenes."))
                            plotEigengeneNetworks(MEs,
                                                  "Eigengene Network",
                                                  marDendro = c(2,3,2,1),
                                                  marHeatmap = c(6,8,2,1),
                                                  plotDendrograms = TRUE,
                                                  plotHeatmaps = TRUE
                                                  )
                          })
                          
                          ##############                            
                          wgcna_gct <- reactive({
                            DataIn <- DataReactive()
                            req(DataIn$data_wide)
                            data_wide <- na.omit(DataIn$data_wide)
                            ProteinGeneName <- DataIn$ProteinGeneName
                            MetaData <- DataIn$MetaData
                            
                            module_map <- tibble::tibble(
                              ME_name = MEs_name_updated(),
                              color = sub("^ME\\d+_", "", MEs_name_updated()),
                            )
                            
                            t2 <- tibble::tibble(
                              UniqueID = names(wgcna$colors),
                              color = labels2colors(wgcna$colors)
                            ) %>% 
                              dplyr::left_join(module_map, by = "color") %>%
                              dplyr::inner_join(ProteinGeneName, by = "UniqueID") %>%
                              dplyr::select(-color)
                            
                            data.in <- DataIn$data_wide[t2$UniqueID,]
                            
                            gct_data <- create_gct_object(data.in, t2, MetaData)
                            gct_data
                          })
                          
                          observeEvent(input$wgcna_gct, {
                            saved_gcts$wgcna_gct <- wgcna_gct()
                          })                               
                          ##############                            
                          
                          if ("sft" %in% names(wgcna_out)) {
                            sft <- wgcna_out$sft
                            output$soft_threshold_table_note <- renderText({
                              HTML("In WGCNA, soft‑thresholding power (β) controls how strongly you emphasize large correlations and suppress small ones when building the adjacency matrix. The table below showed the result for all β candidates tested. <b>The row is highlighted for the picked power β.<b>")
                            })
                            
                            output$soft_threshold_table <- DT::renderDT({
                              datatable(
                                sft$fitIndices,
                                rownames = FALSE, extensions = 'Buttons', escape = FALSE, selection = 'none', class = 'cell-border strip hover',
                                options = list(dom = 'lBfrtip', pageLength = 15,
                                               buttons = list(
                                                 list(extend = "csv", text = "Download Page", filename = "Page_results", 
                                                      exportOptions = list(modifier = list(page = "current"))),
                                                 list(extend = "csv", text = "Download All", filename = "All_Results",
                                                      exportOptions = list(modifier = list(page = "all")))
                                               )
                                )
                              ) %>% 
                                formatSignif(
                                  columns = setdiff(names(Filter(is.numeric, sft$fitIndices)), "Power"),
                                  digits = 3
                                ) %>% formatStyle(columns="Power", target="row", backgroundColor = styleEqual(picked_power, "skyblue", default='white'))
                            })
                            
                            output$soft_threshold_diagnose_plot_note <- renderUI({
                              HTML("
                                In practice, soft‑thresholding power (β) is chosen with a balanced good scale‑free fit and reasonable connectivity.
                                WGCNA recommends choosing the lowest power where R² ≥ 0.8 (or 0.9 for stricter networks). 
                                And typical good range of mean connectivity is 20–200, depending on dataset size.
                                <br><br>
                                
                                <b>Practical decision rules:</b><br>
                                -- If R² rises and then plateaus, pick the lowest β at the plateau.<br>
                                -- If R² rises slowly but connectivity is still good, pick β = 6–10 (signed).<br>
                                -- If R² is flat and connectivity collapses at high β, pick a moderate β (6–8) to avoid oversparsifying.<br>
                                -- If R² keeps rising up to the max tested power, increase the tested range (e.g., 1–30), but still avoid β > 20 unless absolutely necessary.<br>
                                -- Ensure mean connectivity is not too low. It should > ~ 20 (for typical datasets).<br>
                                -- If estimated power is NA, use any number in the typical β range (6–12) for signed networks.
                                     ")
                            })
                            
                            output$soft_threshold_diagnose_plot <- renderPlot({
                              plot_soft_threshold_diagnose(sft, powers())
                            })
                          } else if (input$WGCNA_QC == "Soft-Thresholding Power") {
                            showNotification("Pre-calculated wgcna file does not contain the soft_thresholding power testing result. No results loaded.", duration = 5, type = "warning")
                            # output$soft_threshold_table <- NULL
                            # output$soft_threshold_diagnose_plot <- NULL
                          }
                        })
                        
                        observeEvent(input$runORA_trigger, {
                          gene_list <- input$runORA_trigger
                          updateNavbarPage(parent_session, inputId = "menu", selected = "gsea")
                          updateTabsetPanel(parent_session, inputId = "GS-geneset_tabset", selected = "Over-Representation Analysis (ORA)")
                          parent_session$onFlushed(function() {
                            updateSelectInput(parent_session, inputId = "GS-ORA_input_type", selected = "Gene List")
                            updateTextAreaInput(parent_session, inputId = "GS-ORA_list", value = gene_list)
                          }, once = TRUE)
                        })
                        
                        # Update Module-Trait Relationships menu
                        observe({
                          DataIn = DataReactive()
                          MetaData = DataIn$MetaData
                          req(MetaData)
                          attributes=sort(setdiff(colnames(MetaData), c("sampleid", "Order", "ComparePairs") ))
                          updateSelectizeInput(session, "WGCNA_trait_var", choices=attributes, selected="group")
                        })
                        
                        output$attribute_settings_ui <- renderUI({
                          req(input$WGCNA_trait_var, DataReactive())
                          
                          DataIn = DataReactive()
                          MetaData = DataIn$MetaData
                          # Only keep categorical attributes
                          categorical_attrs <- input$WGCNA_trait_var[
                            sapply(input$WGCNA_trait_var, function(a) is_categorical(MetaData[[a]]))
                          ]
                          
                          # Build UI for each categorical attribute
                          lapply(categorical_attrs, function(attr) {
                            base_choices <- unique(MetaData[[attr]])
                            wellPanel(
                              h4(attr),
                              selectInput(
                                inputId = ns(paste0("base_", attr)),
                                label   = paste("Select base level for", attr),
                                choices = base_choices
                              )
                            )
                          })
                        })
                        
                        # Create trait_data() and calculate moduleTraitCor() whenever wgcna result (MEs_updated()) changes or trait information (input$plot_module_trait) changes
                        observeEvent(list(MEs_updated(),input$plot_module_trait), {
                          req(DataReactive(), MEs_updated())
                          DataIn <- DataReactive()
                          MetaData <- DataIn$MetaData
                          attrs <- input$WGCNA_trait_var
                          # Only categorical attributes have base_* inputs
                          categorical_attrs <- attrs[
                            sapply(attrs, function(a) is_categorical(MetaData[[a]]))
                          ]
                          base_levels <- sapply(categorical_attrs, function(a) input[[paste0("base_", a)]])
                          trait_data(build_traits_matrix(MetaData, attrs, base_levels))

                          nSamples <- nrow(trait_data())
                          moduleTraitCor(cor(MEs_updated(), trait_data(), use = "p"))
                        })
                        
                        # Create Module-Trait correlation heatmap object and reset hub gene results
                        observeEvent(input$plot_module_trait, {
                          req(MEs_updated(),trait_data())
                          p <- compute_module_trait_heatmap(
                            MEs = MEs_updated(),
                            trait_mat = trait_data()
                          ) 
                          module_trait_plot(p)
                          # reset the dependent data
                          output$hub_gene_table <- DT::renderDT({NULL})
                          output$plot_MMvsGS <- renderPlotly({NULL})
                          output$plot_MMvsConnectivity <- renderPlotly({NULL})
                        })
                        
                        output$module_trait_hmap <- renderPlotly({
                          module_trait_plot()
                        })
                        
                        observeEvent(event_data("plotly_click", source = "moduleTraitHeatmap"), {
                          req(input$plot_module_trait,module_trait_plot())
                          click <- event_data("plotly_click", source = "moduleTraitHeatmap")
                          req(click)
                          updateSelectInput(session, "WGCNA_trait",  selected = click$x)
                          updateSelectInput(session, "WGCNA_module", selected = click$y)
                          updateTabsetPanel(session, "Module_Trait", selected = "Hub Gene Identification Table")
                          }, ignoreNULL = TRUE)

                        # Update trait and module menu whenever the trait_data() changes.
                        observeEvent(list(trait_data(),input$plot_module_trait), {
                          updateSelectInput(session, "WGCNA_trait", choices=names(trait_data()), selected=character(0))
                        })
                        observeEvent(list(MEs_updated(),input$plot_module_trait), {
                          updateSelectInput(session, "WGCNA_module", choices=names(MEs_updated()), selected=character(0))
                        })
                        
                        observeEvent(input$plot_module_hub, {
                          withProgress(message = 'Processing...', value = 0, {
                            req(moduleTraitCor())
                            if (is.null(input$WGCNA_trait) || input$WGCNA_trait == "") { 
                              showNotification("Selecting a trait is required to run the analysis.", type = "error") 
                              return() 
                            }
                            
                            # Retrive pre-computed wgcna result (load_wgcna_file) or on-the-fly result(WGCNAReactive()) 
                            # compute hub gene results
                            wgcna_out <- tryCatch(WGCNAReactive(), error = function(e) NULL)
                            if (is.null(wgcna_out)) {
                              req(ProjectInfo,DataReactive())
                              ProjectID <- ProjectInfo$ProjectID
                              wgcnafile <- ProjectInfo$file3
                              load(wgcnafile)
                              wgcna <- netwk
                              load_wgcna_file <- file.path(dirname(wgcnafile), sub("^wgcna", "load", basename(wgcnafile)))
                              load(load_wgcna_file)
                            } else {
                              picked_power <- wgcna_out$picked_power
                              dataExpr <- wgcna_out$dataExpr
                              wgcna <- wgcna_out$netwk
                            }
                            gene_label <- input$WGCNAgenelable
                            DataIn = DataReactive()
                            req(DataIn$ProteinGeneName)
                            ProteinGeneName  <- DataIn$ProteinGeneName
                            
                            selected_trait <- input$WGCNA_trait     # trait is required 
                            
                            # module is optional, if not selected, then identify the one with strongest correlation
                            if (!is.null(input$WGCNA_module) && input$WGCNA_module != "") {
                              selected_trait_module <- input$WGCNA_module
                            } else {
                              # Get the module most strongly associated with the trait
                              module_trait_cor <- moduleTraitCor()
                              selected_trait_correlations <- abs(module_trait_cor[, selected_trait])
                              selected_trait_module <- MEs_name_updated()[which.max(selected_trait_correlations)]
                            }
                            
                            # Convert numeric module label to color name
                            parts <- strsplit(selected_trait_module, "_")[[1]]
                            numeric_label   <- sub("ME", "", parts[1])   # number
                            module_name <- parts[2]                  # color
                            
                            # Get genes in this module
                            module_genes <- tibble::tibble(
                              UniqueID = names(wgcna$colors),
                              color = labels2colors(wgcna$colors)
                            ) %>%
                              dplyr::filter(color == module_name) %>%
                              pull(UniqueID) %>% 
                              strsplit(",") %>% 
                              unlist() %>% 
                              trimws()  
                            
                            module_genes <- module_genes[module_genes != ""]
                            module_genes <- intersect(module_genes, colnames(dataExpr))
                            
                            # Calculate connectivity 
                            # First, calculate adjacency matrix for genes in this module
                            adjacency_matrix <- adjacency(
                              dataExpr[, module_genes],
                              power = picked_power,
                              type = "signed"
                            )
                            # Calculate connectivity: sum of connection weights for each gene (subtract 1 to exclude self-connection)
                            connectivity <- rowSums(adjacency_matrix) - 1
                            
                            # Calculate module membership for these genes
                            # MM measures how correlated each gene is with the module eigengene
                            gene_module_membership <- cor(dataExpr[, module_genes],
                                                          MEs_updated()[, selected_trait_module],
                                                          use = "pairwise.complete.obs")
                            
                            # Calculate gene significance (correlation with trait)
                            nSamples = nrow(trait_data())
                            gene_trait_cor <- cor(dataExpr, trait_data()[ , selected_trait], use = "pairwise.complete.obs")
                            gene_trait_pvalue <- corPvalueStudent(as.numeric(gene_trait_cor), nSamples)
                            
                            # Add gene names to the p-value vector
                            names(gene_trait_cor) <- colnames(dataExpr)
                            names(gene_trait_pvalue) <- colnames(dataExpr)
                            
                            if (! all(names(dataExpr) %in% ProteinGeneName[, gene_label])) {
                              current_label <- setdiff(c("UniqueID","Gene.Name"), gene_label)
                              id_to_gene <- setNames(ProteinGeneName[ , gene_label], ProteinGeneName[ , current_label])
                              module_genes_display = id_to_gene[module_genes]
                            } else {
                              module_genes_display = module_genes
                            }
                            
                            # Combine all metrics for hub gene identification
                            hub_gene_info <- data.frame(
                              Gene = module_genes_display,
                              Connectivity = connectivity,
                              ModuleMembership = as.numeric(gene_module_membership),
                              MM_pvalue = corPvalueStudent(as.numeric(gene_module_membership), nSamples),
                              GeneSignificance = gene_trait_cor[module_genes],
                              GS_pvalue = gene_trait_pvalue[module_genes],
                              stringsAsFactors = FALSE
                            ) %>% 
                              dplyr::filter(Gene != "")
                            
                            # Sort by connectivity to identify hubs
                            hub_gene_info <- hub_gene_info[order(-hub_gene_info$Connectivity), ]
                            
                            df_hub_gene(hub_gene_info)
                          })
                          
                          output$hub_gene_table <- DT::renderDT({
                              DT::datatable(df_hub_gene(),
                                            rownames = FALSE, 
                                            extensions = 'Buttons', 
                                            escape = FALSE, 
                                            selection = 'none', 
                                            class = 'cell-border strip hover',
                                            options = list(dom = "Blfrtip", 
                                                           buttons = c("csv", "excel", "print")
                                            )
                              ) %>% 
                                formatSignif(columns=names(Filter(is.numeric, df_hub_gene())), digits=3)
                            })
                            
                          observeEvent(input$hub_gene, {
                            saved_table$hub_gene <- df_hub_gene()
                          })
                          
                          output$plot_MMvsGS <- renderPlotly({
                              module_gene_info <- df_hub_gene()
                              
                              # Sort
                              module_gene_info <- module_gene_info[order(-abs(module_gene_info$ModuleMembership)), ]
                              
                              # x and y
                              x <- abs(module_gene_info$ModuleMembership)
                              y <- abs(module_gene_info$GeneSignificance)
                              
                              # Regression
                              fit <- lm(y ~ x)
                              x_seq <- seq(min(x), max(x), length.out = 100)
                              y_pred <- predict(fit, newdata = data.frame(x = x_seq))
                              
                              # Correlation
                              mm_gs_cor <- cor(x, y, use = "pairwise.complete.obs")
                              
                              char_df <- module_gene_info
                              
                              # Format numeric columns to 3 decimals (or whatever you want)
                              num_cols <- sapply(char_df, is.numeric)
                              char_df[num_cols] <- lapply(char_df[num_cols], function(x) sprintf("%.3f", x))
                              
                              # Convert everything to character
                              char_df <- data.frame(lapply(char_df, as.character), stringsAsFactors = FALSE)
                              
                              cols <- colnames(char_df)
                              
                              hover_text <- apply(char_df, 1, function(row) {
                                paste0("<b>", cols, ":</b> ", row, collapse = "<br>")
                              })
                              
                              # Build hover template dynamically
                              cols <- colnames(module_gene_info)
                              hover_lines <- paste0(
                                "<b>", cols, ":</b> %{customdata[", seq_along(cols) - 1, "]}", 
                                collapse = "<br>"
                              )
                              hover_template <- paste0(hover_lines, "<extra></extra>")
                              
                              plot_ly() %>%
                                add_markers(
                                  x = x,
                                  y = y,
                                  type = "scatter",
                                  mode = "markers", # Added this to ensure clean rendering
                                  marker = list(color = module_name, size = 8),
                                  text = hover_text,                       # <--- Pass your pre-built R strings here
                                  hovertemplate = "%{text}<extra></extra>", # <--- Just tell Plotly to show the text
                                  name = 'Gene'
                                ) %>%
                                add_lines(
                                  x = x_seq,
                                  y = y_pred,
                                  line = list(color = "red", width = 2),
                                  name = "Regression"
                                ) %>%
                                layout(
                                  title = paste("MM vs GS in", selected_trait_module, "module"),
                                  xaxis = list(title = paste("Module Membership in", selected_trait_module, "module")),
                                  yaxis = list(title = paste("Gene Significance for", selected_trait)),
                                  annotations = list(
                                    list(
                                      x = min(x),
                                      y = max(y),
                                      text = paste("cor =", round(mm_gs_cor, 3)),
                                      xanchor = "left",
                                      yanchor = "top",
                                      showarrow = FALSE,
                                      font = list(size = 14)
                                    )
                                  )
                                )
                            })
                          
                          output$plot_MMvsConnectivity <- renderPlotly({
                            hub_gene_info <- df_hub_gene()
                            
                            # Identify top 10% connected
                            top_cutoff <- quantile(hub_gene_info$Connectivity, 0.9)
                            group <- ifelse(hub_gene_info$Connectivity > top_cutoff,
                                            "Top 10% connected", "Other genes")
                            
                            # Build full-row hover text
                            char_df <- hub_gene_info
                            num_cols <- sapply(char_df, is.numeric)
                            char_df[num_cols] <- lapply(char_df[num_cols], function(x) sprintf("%.3f", x))
                            char_df <- data.frame(lapply(char_df, as.character), stringsAsFactors = FALSE)
                            cols <- colnames(char_df)
                            
                            hover_text <- apply(char_df, 1, function(row) {
                              paste0("<b>", cols, ":</b> ", row, collapse = "<br>")
                            })
                            
                            plot_ly(
                              data = hub_gene_info,
                              x = ~ModuleMembership,
                              y = ~Connectivity,
                              type = "scatter",
                              mode = "markers",
                              color = ~group,
                              colors = c("Other genes" = "black", "Top 10% connected" = "red"),
                              text = hover_text,
                              hovertemplate = "%{text}<extra></extra>",
                              marker = list(size = 8)
                            ) %>%
                              layout(
                                title = paste("Hub Gene Identification in", selected_trait_module, "Module"),
                                xaxis = list(title = "Module Membership"),
                                yaxis = list(title = "Connectivity (Intramodular)"),
                                legend = list(
                                  x = 0,
                                  y = 1,
                                  xanchor = "left",
                                  yanchor = "top",
                                  orientation = "v",
                                  font = list(size = 16)
                                )
                              )
                          })
                        })
                      })
}