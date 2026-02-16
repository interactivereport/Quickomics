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

get_wgcna_table <-function(wgcna, ProteinGeneName, gene_label) {
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
  dplyr::select(color, n_gene, copy, gene_group)
  return(t2)
}

wgcna_ui <- function(id) {
  ns <- shiny::NS(id)
  fluidRow(
    rclipboard::rclipboardSetup(),
    column(3,
           wellPanel(
             uiOutput(ns('loadedprojects')),
             radioButtons(ns("WGCNAgenelable"),label="Select Gene Label",inline = TRUE, choices=c("Gene.Name","UniqueID"), selected="Gene.Name"),
             numericInput(ns("WGCNAtopNum"), label= "Select Top N Genes, where N is :",  value=250L, min=250L, step=25L, max = 10000L),
             numericInput(ns("mergeCutHeight"), label= "Dendrogram Cut Height for Merging:",  value=0.25, min= 0, max = 1.0, step = 0.01),
             actionButton(ns("plotwgcna"),"Run"),
             br(),
             strong("Running the data could take 3-10 minutes",style="color:red"), span("depending on data size and complexity; once Run starts, please refrain from clicking the button repeatedly.",style="color:red", inline = TRUE),
           )
    ),
    column(9,
           tabsetPanel(id=ns("WGCNA_tabset"),
                       tabPanel(title="Dendrogram", value=ns("Dendrogram"),
                                plotOutput(ns("Dendrogram"), height=800)
                       ),
                       tabPanel(title="Gene Clusters",
                                DT::dataTableOutput(ns("gene_cluster")),
                                tags$script(HTML(sprintf("
                                  $(document).on('click', '[id^=runORA_]', function() {
                                    var genes = $(this).attr('data-genes');
                                    Shiny.setInputValue('%s', genes, {priority: 'event'});
                                  });
                                ", ns("runORA_trigger"))))
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
                        working_project=reactiveVal()
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
                          wgcnafile <- ProjectInfo$file3
                          gene_label <- input$WGCNAgenelable

                          if(file.exists(wgcnafile)){
                            DataIn = DataReactive()
                            req(DataIn$ProteinGeneName)
                            ProteinGeneName  <- DataIn$ProteinGeneName
                            
                            load(wgcnafile)
                            wgcna <- netwk
                            mergedColors = labels2colors(wgcna$colors)

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
                            t2 <- get_wgcna_table(wgcna, ProteinGeneName, gene_label)
                            
                            output$gene_cluster <- DT::renderDT({
                              DT::datatable(
                                t2,
                                escape = FALSE,
                                selection = "none",
                                colnames=c("Color of cluster", "Number of genes", "Action","Genes in cluster")
                              )
                            })
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
                              powers <- c(1L:10L, seq(from = 12L, to = 20L, by = 2L))
                              t2 <- Sys.time()
                              cor <- stats::cor
                              sft <- WGCNA::pickSoftThreshold(dataExpr, dataIsExpr = TRUE, powerVector = powers, corFnc = cor, corOptions = list(use = 'p'), networkType = "signed")
                              
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
                            netwk
                          })
                        })
                        
                        #### generate dendrogram and gene cluster table #####
                        # use input$WGCNAReactive() as event handler to ensure observeEvent() depends on it only
                        # and does not directly depends on input$, which ensure WGCNAReactive() will be calculated first.
                        observeEvent(WGCNAReactive(),{
                          req(DataIn = DataReactive())
                          wgcna <- WGCNAReactive()
                          DataIn = DataReactive()
                          mergedColors = labels2colors(wgcna$colors)
                          
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
                          
                          t2 <- get_wgcna_table(wgcna, ProteinGeneName, gene_label)
                          
                          output$gene_cluster <- DT::renderDT({
                            DT::datatable(
                              t2,
                              escape = FALSE,
                              selection = "none",
                              colnames=c("Color of cluster", "Number of genes", "Action","Genes in cluster")
                            )
                          })
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
                      }
  )
}