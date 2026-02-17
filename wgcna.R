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
    dplyr::left_join(module_map, by = "color") %>%   # <-- add module label
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
  numeric_traits <- MetaData[, sapply(MetaData, is.numeric), drop = FALSE]
  if (ncol(factor_traits) > 0) {
    dummy_traits <- model.matrix(~ . - 1, data = factor_traits)
    dummy_traits <- as.data.frame(dummy_traits)
  } else {
    dummy_traits <- NULL
  }
  cbind(numeric_traits, dummy_traits)
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
    hoverinfo = "text"
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
                              strong("Running the data could take 3-10 minutes",style="color:red"), span("depending on data size and complexity; once Run starts, please refrain from clicking the button repeatedly.",style="color:red", inline = TRUE),
             ),
             conditionalPanel(ns = ns, "input.WGCNA_tabset=='Module-Trait Relationships'",
                              conditionalPanel(ns = ns, "input.Module_Trait=='Module Trait heatmap'",
                                               selectizeInput(ns("WGCNA_trait_var"),label="Select attributes", choices = NULL,multiple = TRUE,options = list(placeholder = "Choose one or more attributes")),
                                               uiOutput(ns("attribute_settings_ui")),
                                               actionButton(ns("plot_module_trait"),"Plot Module Trait Correlation")
                              ),
                              conditionalPanel(ns = ns, "input.Module_Trait=='Hub Gene Identification Table' || input.Module_Trait == 'Hub Gene Identification Plot'",
                                               selectInput(ns("WGCNA_trait"),label="Select a trait", choices = NULL ,multiple = FALSE),
                                               actionButton(ns("plot_module_hub"),"Run")
                              )
             )
           )
    ),
    column(9,
           tabsetPanel(id=ns("WGCNA_tabset"),
                       tabPanel(title="WGCNA Result",
                                tabsetPanel(id=ns("WGCNA_Result"),
                                            tabPanel(title="Dendrogram", plotOutput(ns("Dendrogram"), height=800)),
                                            tabPanel(title="Gene Clusters",
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
                                            tabPanel(title="Eigengene table", DT::dataTableOutput(ns("MEs"))),
                                            tabPanel(title = "Eigengene Network",plotOutput(ns("Eigenene_Network"), height = "1000px"))
                                )
                       ),
                       tabPanel(title="Module-Trait Relationships",
                                tabsetPanel(id=ns("Module_Trait"),
                                            tabPanel(title="Module Trait heatmap", 
                                                     plotlyOutput(ns("module_trait_hmap"), height=800)
                                            ),
                                            tabPanel(title="Hub Gene Identification Table",
                                                     DT::dataTableOutput(ns("hub_gene_table"))
                                            ),
                                            tabPanel(
                                              title = "Hub Gene Identification Plot",
                                              div(
                                                style = "padding-top: 30px; padding-bottom: 40px;",
                                                plotlyOutput(ns("plot_MMvsGS"), height = "1000px")
                                              ),
                                              div(
                                                style = "padding-top: 40px;",
                                                plotlyOutput(ns("plot_MMvsConnectivity"), height = "1000px")
                                              )
                                            )
                                )
                       ),
                       tabPanel(title="WGCNA QC",
                                tabsetPanel(id=ns("WGCNA_QC"),
                                            tabPanel(title="Soft-Thresholding Power", 
                                                     DT::dataTableOutput(ns("soft_threshold_table")),
                                                     plotOutput(ns("soft_threshold_diagnose_plot"), height=800)
                                            ),
                                            tabPanel(title="WGCNA Network",
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
                        wgcna_run_control<-reactiveVal(0)
                        
                    
                        observe({
                          req(ProjectInfo)
                          working_project(ProjectInfo$ProjectID)
                        })
                        
                        observe({
                          DataIn <- DataReactive()
                          data_wide <- na.omit(DataIn$data_wide)
                          default_n_gene <- min(10000, nrow(data_wide))
                          
                          updateNumericInput(session, "WGCNAtopNum", 
                                             label= "Select Top N Genes, where N is :",  value=default_n_gene, min=250L, step=25L, max = default_n_gene)
                        })
                        
                        observeEvent(list(working_project(),input$WGCNAgenelable,parent_session$input$menu == "wgcna"), {
                          req(ProjectInfo,DataReactive(),parent_session$input$menu == "wgcna")
                          ProjectID <- ProjectInfo$ProjectID
                          wgcnafile <- ProjectInfo$file3
                          gene_label <- input$WGCNAgenelable
                          
                          load_wgcna_file <- paste("data/wgcna_data/load_", ProjectID, ".RData", sep = "")
                          
                          if(file.exists(wgcnafile) & file.exists(load_wgcna_file)){
                            DataIn = DataReactive()
                            req(DataIn$ProteinGeneName)
                            ProteinGeneName  <- DataIn$ProteinGeneName
                            load(wgcnafile)
                            load(load_wgcna_file)
                            wgcna <- netwk
                            mergedColors = labels2colors(wgcna$colors)
                            
                            MEs <- wgcna$MEs
                            moduleColors <- wgcna$colors
                            MEs_updated(rename_MEs_blockwise(MEs, moduleColors))
                            MEs_name_updated(colnames(MEs_updated()))

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
                              DT::datatable(
                                t2,
                                escape = FALSE,
                                selection = "none",
                                colnames=c("Cluster", "Number of genes", "Action","Genes in cluster")
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

                            if (exists("sft")) {
                              output$soft_threshold_table <- DT::renderDT({
                                DT::datatable(
                                  sft$fitIndices, 
                                  rownames = FALSE,  extensions = 'Buttons', escape = FALSE, selection = 'none', class = 'cell-border strip hover',
                                  options = list(    dom = 'lBfrtip', pageLength = 15,
                                                     buttons = list(
                                                       list(extend = "csv", text = "Download Page", filename = "Page_results",
                                                            exportOptions = list(modifier = list(page = "current"))),
                                                       list(extend = "csv", text = "Download All", filename = "All_Results",
                                                            exportOptions = list(modifier = list(page = "all")))
                                                     )
                                  )) %>% 
                                  formatSignif(columns=names(Filter(is.numeric, sft$fitIndices)), digits=3)
                              })
                              
                              output$soft_threshold_diagnose_plot <- renderPlot({
                                plot_soft_threshold_diagnose(sft, powers())
                              })
                            } else if (input$WGCNA_tabset == 'WGCNA_QC' && input$WGCNA_QC == "Soft-Thresholding Power") {
                              showNotification("Pre-calculated wgcna file does not contain the soft_thresholding power testing result. No results loaded.", duration = 5, type = "warning")
                              output$soft_threshold_table <- NULL
                              output$soft_threshold_diagnose_plot <- NULL
                            }
                          } else if (parent_session$input$menu == "wgcna") {
                            print("no pre-computed wgcna file available and cannot load wgcna results")
                            showNotification("Cannot find pre-calculated wgcna file, no WGCNA results loaded. ", duration = 5, type = "warning")
                            output$Dendrogram <- NULL
                            output$gene_cluster <- NULL
                          }
                        })
                        
                        # use eventReactive to control reactivity of WGCNAReactive;
                        # otherwise, whenever an input change, WGCNAReactive will be re-calculated
                        # and its re-calculation could take a long time.
                        WGCNAReactive <- eventReactive(input$plotwgcna, {
                          withProgress(message = "Running WGCNA", detail = 'This may take a while...', value = 0.2, {
                            # what if the user-imported data doesn't have $data_wide, $ProjectID..etc?
                            req(ProjectInfo, DataReactive(), ProjectInfo$ProjectID)

                            DataIn = DataReactive()
                            data_wide <- DataIn$data_wide
                            ProjectID <- ProjectInfo$ProjectID
                            
                            data_wide <- na.omit(DataIn$data_wide)
                            diff <- apply(data_wide, 1, sd, na.rm = TRUE)/(rowMeans(data_wide) + median(rowMeans(data_wide)))
                            data_wide=data_wide[order(diff, decreasing=TRUE), ]                            

                            if (nrow(data_wide)>10000 ) {
                              data_wide=data_wide[1:10000, ] 
                              cat("reduce gene size to 10K for project ", ProjectID, "\n")
                            } 
                            # dataExpr_ori <- data_wide
                            
                            print(paste0("**** dim of dataExpr after-preprocssing is ****", dim(data_wide)))
                            
                            default_n_gene <- nrow(data_wide)
                            
                            # Note: if launching app from the server, the path for `load_`  files should be
                            # paste0("/mnt/depts/dept04/compbio/projects/xOmicsShiny/data/wgcna_data/TOM
                            load_wgcna_file <- paste("data/wgcna_data/load_", ProjectID, ".RData", sep = "")
                            
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
                              # Scenario 3: Not scenario 1 or 2, and recalcuate everything
                              print(paste0("**** compute everything from scratch ****"))
                              ProteinGeneName  <- DataIn$ProteinGeneName

                              ## Top number of genes
                              topNum <- as.numeric(input$WGCNAtopNum)
                              # Gene Label
                              gene_label <- input$WGCNAgenelable
                              dataExpr <- data_wide %>%
                                na.omit() %>%
                                tibble::rownames_to_column("UniqueID") %>%
                                dplyr::left_join(ProteinGeneName, by = "UniqueID") %>%
                                dplyr::select(-id, -Gene.Name, -Protein.ID) %>%
                                tibble::column_to_rownames("UniqueID")
                              
                              # Replace Inf or NA with NA
                              dataExpr[dataExpr == "Inf" | is.na(dataExpr)] <- NA
                              
                              gene.names <- rownames(dataExpr)                              
                              SubGeneNames=gene.names[1L:topNum]
                              
                              # Ensure all columns are numeric before transposing; otherwise cell values may
                              # become character, causing problems in WGCNA::blockwiseModules, as happened to
                              # the Mouse_microglia_RNA-Seq data
                              dataExpr <-  dataExpr %>%
                                dplyr::select(tidyselect::where(is.numeric))
                              
                              dataExpr = as.data.frame(t(dataExpr))
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
                        
                        observeEvent(input$plotwgcna, {
                          wgcna_run_control(wgcna_run_control()+1)
                        })
                        
                        
                        #### generate dendrogram and gene cluster table #####
                        # use input$WGCNAReactive() as event handler to ensure observeEvent() depends on it only
                        # and does not directly depends on input$, which ensure WGCNAReactive() will be calculated first.
                        observeEvent(wgcna_run_control(),{
                          req(DataReactive())
                          wgcna_out <- WGCNAReactive()
                          wgcna <- wgcna_out$netwk
                          picked_power <- wgcna_out$picked_power
                          DataIn = DataReactive()
                          mergedColors = labels2colors(wgcna$colors)
                          
                          MEs <- wgcna$MEs
                          moduleColors <- wgcna$colors
                          MEs_updated <- rename_MEs_blockwise(MEs, moduleColors)
                          MEs_name_updated(colnames(MEs_updated))

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
                            DT::datatable(
                              t2,
                              escape = FALSE,
                              selection = "none",
                              colnames=c("Cluster", "Number of genes", "Action","Genes in cluster")
                            )
                          })
                          
                          MEs <- wgcna$MEs
                          moduleColors <- wgcna$colors
                          MEs_updated(rename_MEs_blockwise(MEs, moduleColors))
                          MEs_name_updated <- colnames(MEs_updated())
                          
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
                            
                            plotEigengeneNetworks(MEs,
                                                  "Eigengene Network",
                                                  marDendro = c(2,3,2,1),
                                                  marHeatmap = c(6,8,2,1),
                                                  plotDendrograms = TRUE,
                                                  plotHeatmaps = TRUE
                                                  )
                          })
                          
                          
                          if ("sft" %in% names(wgcna_out)) {
                            sft <- wgcna_out$sft
                            output$soft_threshold_table <- DT::renderDT({
                              DT::datatable(
                                sft$fitIndices, rownames = FALSE, extensions = 'Buttons', escape = FALSE, selection = 'none', class = 'cell-border strip hover',
                                options = list(    dom = 'lBfrtip', pageLength = 15,
                                                   buttons = list(
                                                     list(extend = "csv", text = "Download Page", filename = "Page_results",
                                                          exportOptions = list(modifier = list(page = "current"))),
                                                     list(extend = "csv", text = "Download All", filename = "All_Results",
                                                          exportOptions = list(modifier = list(page = "all")))
                                                   )
                                )) %>% 
                                formatSignif(columns=names(Filter(is.numeric, sft$fitIndices)), digits=3)
                            })
                            
                            output$soft_threshold_diagnose_plot <- renderPlot({
                              plot_soft_threshold_diagnose(sft, powers())
                            })
                          } else if (input$WGCNA_QC == "Soft-Thresholding Power") {
                            showNotification("Pre-calculated wgcna file does not contain the soft_thresholding power testing result. No results loaded.", duration = 5, type = "warning")
                            output$soft_threshold_table <- NULL
                            output$soft_threshold_diagnose_plot <- NULL
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
                        
                        observe({
                          DataIn = DataReactive()
                          MetaData = DataIn$MetaData
                          attributes=sort(setdiff(colnames(MetaData), c("sampleid", "Order", "ComparePairs") ))
                          updateSelectizeInput(session, "WGCNA_trait_var", choices=attributes, selected="group")
                          
                        })
                        
                        output$attribute_settings_ui <- renderUI({
                          req(input$WGCNA_trait_var)
                          
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
                        
                        WGCNA_mm_trait_Reactive <- reactive({
                          req(DataReactive(), MEs_updated(),input$plot_module_trait, wgcna_run_control())

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
                          moduleTraitPvalue <- corPvalueStudent(moduleTraitCor(), nSamples)

                          p <- make_plotly_heatmap(
                            cor_mat = moduleTraitCor(),
                            p_mat = moduleTraitPvalue,
                            row_labels = colnames(MEs_updated()),
                            col_labels = colnames(trait_data())
                          )
                          p
                        })

                        output$module_trait_hmap <- renderPlotly({
                          WGCNA_mm_trait_Reactive()
                        })
                        
                        observeEvent(trait_data(), {
                          updateSelectInput(session, "WGCNA_trait", choices=names(trait_data()), selected=character(0))
                        })
                        
                        observeEvent(input$plot_module_hub, {
                          browser()
                          req(WGCNA_mm_trait_Reactive())
                          
                          
                          wgcna_out <- tryCatch(WGCNAReactive(), error = function(e) NULL)
                          if (is.null(wgcna_out)) {
                            req(ProjectInfo,DataReactive())
                            ProjectID <- ProjectInfo$ProjectID
                            load_wgcna_file <- paste("data/wgcna_data/load_", ProjectID, ".RData", sep = "")
                            load(load_wgcna_file)
                          } else {
                            picked_power <- wgcna_out$picked_power
                            dataExpr <- wgcna_out$dataExpr
                          }
                          gene_label <- input$WGCNAgenelable
                          DataIn = DataReactive()
                          req(DataIn$ProteinGeneName)
                          ProteinGeneName  <- DataIn$ProteinGeneName
                          if (! all(names(dataExpr) %in% ProteinGeneName[, gene_label])) {
                            current_label <- setdiff(c("UniqueID","Gene.Name"), gene_label)
                            id_to_gene <- setNames(ProteinGeneName[ , gene_label], ProteinGeneName[ , current_label])
                            new_names <- id_to_gene[colnames(dataExpr)]
                            colnames(dataExpr) <- new_names
                          }
                          
                          selected_trait <- input$WGCNA_trait
                          module_trait_cor <- moduleTraitCor()
                          # Get the module most strongly associated with T_AKO
                          selected_trait_correlations <- abs(module_trait_cor[, selected_trait])
                          selected_trait_module <- MEs_name_updated()[which.max(selected_trait_correlations)]

                          # Convert numeric module label to color name
                          parts <- strsplit(selected_trait_module, "_")[[1]]
                          numeric_label   <- sub("ME", "", parts[1])   # number
                          module_name <- parts[2]                  # color
                          
                          # Get genes in this module
                          module_genes <- df_gene_clusters() %>%
                            dplyr::filter(color == selected_trait_module) %>%
                            dplyr::pull(gene_group) %>% 
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
                          
                          # Combine all metrics for hub gene identification
                          hub_gene_info <- data.frame(
                            Gene = module_genes,
                            Connectivity = connectivity,
                            ModuleMembership = as.numeric(gene_module_membership),
                            MM_pvalue = corPvalueStudent(as.numeric(gene_module_membership), nSamples),
                            GeneSignificance = gene_trait_cor[module_genes],
                            GS_pvalue = gene_trait_pvalue[module_genes],
                            stringsAsFactors = FALSE
                          )
                          
                          # Sort by connectivity to identify hubs
                          hub_gene_info <- hub_gene_info[order(-hub_gene_info$Connectivity), ]
                          
                          df_hub_gene(hub_gene_info)
                          
                          output$hub_gene_table <- DT::renderDT({
                            DT::datatable(
                              df_hub_gene(),rownames = FALSE, extensions = 'Buttons', escape = FALSE, selection = 'none', class = 'cell-border strip hover',
                              options = list(    dom = 'lBfrtip', pageLength = 15,
                                                 buttons = list(
                                                   list(extend = "csv", text = "Download Page", filename = "Page_results",
                                                        exportOptions = list(modifier = list(page = "current"))),
                                                   list(extend = "csv", text = "Download All", filename = "All_Results",
                                                        exportOptions = list(modifier = list(page = "all")))
                                                 )
                              )) %>% 
                              formatSignif(columns=names(Filter(is.numeric, df_hub_gene())), digits=3)
                          })
                          
                          # output$plot_MMvsGS <- renderPlot({
                          #   module_gene_info <- df_hub_gene()
                          #   module_gene_info <- module_gene_info[order(-abs(module_gene_info$ModuleMembership)), ]
                          #   
                          #   plot(abs(module_gene_info$ModuleMembership), 
                          #        abs(module_gene_info$GeneSignificance),
                          #        xlab = paste("Module Membership in", selected_trait_module, "module"),
                          #        ylab = paste("Gene Significance for", selected_trait), 
                          #        main = paste("MM vs GS in", selected_trait_module, "module"),
                          #        pch = 20, 
                          #        col = module_name,
                          #        cex = 1.2)
                          #   
                          #   # Add regression line
                          #   abline(lm(abs(module_gene_info$GeneSignificance) ~ abs(module_gene_info$ModuleMembership)), 
                          #          col = "red", lwd = 2)
                          #   
                          #   # Calculate and display correlation
                          #   mm_gs_cor <- cor(abs(module_gene_info$ModuleMembership), 
                          #                    abs(module_gene_info$GeneSignificance),
                          #                    use = "pairwise.complete.obs")
                          #   text(x = 0.2, y = max(abs(module_gene_info$GeneSignificance)) * 0.95,
                          #        labels = paste("cor =", round(mm_gs_cor, 3)),
                          #        pos = 4, cex = 1.2)
                          # })
                          
                          output$plot_MMvsGS <- renderPlotly({
                            module_gene_info <- df_hub_gene()
                            module_gene_info <- module_gene_info[order(-abs(module_gene_info$ModuleMembership)), ]
                            
                            x <- abs(module_gene_info$ModuleMembership)
                            y <- abs(module_gene_info$GeneSignificance)
                            
                            # Regression line
                            fit <- lm(y ~ x)
                            x_seq <- seq(min(x), max(x), length.out = 100)
                            y_pred <- predict(fit, newdata = data.frame(x = x_seq))
                            
                            # Correlation
                            mm_gs_cor <- cor(x, y, use = "pairwise.complete.obs")
                            
                            plot_ly() %>%
                              add_markers(
                                x = x,
                                y = y,
                                marker = list(color = module_name, size = 8),
                                text = module_gene_info$Gene,
                                hoverinfo = "text"
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
                          
                          # output$plot_MMvsConnectivity <- renderPlot({
                          #   hub_gene_info <- df_hub_gene()
                          #   plot(hub_gene_info$ModuleMembership, 
                          #        hub_gene_info$Connectivity,
                          #        xlab = "Module Membership",
                          #        ylab = "Connectivity (Intramodular)",
                          #        main = paste("Hub Gene Identification in", selected_trait_module, "Module"),
                          #        pch = 20,
                          #        col = ifelse(hub_gene_info$Connectivity > quantile(hub_gene_info$Connectivity, 0.9),
                          #                     "red", "black"))
                          #   legend("topleft", 
                          #          legend = c("Top 10% connected", "Other genes"),
                          #          col = c("red", "black"),
                          #          pch = 20)
                          # })
                          output$plot_MMvsConnectivity <- renderPlotly({
                            hub_gene_info <- df_hub_gene()
                            
                            # Identify top 10% connected
                            top_cutoff <- quantile(hub_gene_info$Connectivity, 0.9)
                            colors <- ifelse(hub_gene_info$Connectivity > top_cutoff, "red", "black")
                            
                            plot_ly(
                              x = hub_gene_info$ModuleMembership,
                              y = hub_gene_info$Connectivity,
                              type = "scatter",
                              mode = "markers",
                              marker = list(color = colors, size = 8),
                              text = hub_gene_info$Gene,
                              hoverinfo = "text"
                            ) %>%
                              layout(
                                title = paste("Hub Gene Identification in", selected_trait_module, "Module"),
                                xaxis = list(title = "Module Membership"),
                                yaxis = list(title = "Connectivity (Intramodular)"),
                                legend = list(orientation = "h"),
                                shapes = list()  # placeholder if you want to add lines later
                              ) %>%
                              add_trace(
                                x = NA, y = NA, mode = "markers",
                                marker = list(color = "red"),
                                name = "Top 10% connected"
                              ) %>%
                              add_trace(
                                x = NA, y = NA, mode = "markers",
                                marker = list(color = "black"),
                                name = "Other genes"
                              )
                          })
                          
                        })
                      }
  )
}