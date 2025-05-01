###########################################################################################################
## Omics Visualization R Shiny App
##
##This software belongs to Biogen Inc. All right reserved.
##
##@file: genesetmodule.R
##@Developer : Benbo Gao (benbo.gao@Biogen.com); Xinmin Zhang (xinmin@bioinforx.com)
##@Date : 7/3/2023
##@version 3.0

##########################################################################################################
## Gene Set Enrichment
##########################################################################################################
# pkgs: "pathview", "biomaRt",  "ComplexHeatmap",  "DT", "stringr",  "dplyr", "tidyr"

# this module require results_long, tests, ProteinGeneName, Species, data_long(for heatmap)
library(tibble)
library(pathview)
library(biomaRt)
#library(enrichplot)
#library(enrichR)
library(ComplexHeatmap)
library(fgsea)
#library(wikiprofiler)
library(svgPanZoom)
source("wikiprofiler_Fixed.R")
library(org.Hs.eg.db); library(org.Mm.eg.db); library(org.Rn.eg.db)

#Use metabaser only when the gmt file exists
if (file.exists("db/human/metabase_maps_genesymbols.gmt")) {
  library(metabaser) 
  #Now launch metabaser and connect to database
  if (!metabase.alive()){
    library(metabaser) 
    load("db/metabase_config.RData") 
    metabase.connect(dbname, uid, pwd, host=host, port=port, type=type, driver=driver) 
  }
  MetabaseOutput <- function(outputId, width = "100%", height = "1200px") {
    htmlwidgets::shinyWidgetOutput(outputId, "cyBasicMap", width, height, package = "metabaser")
  }
  renderMetabase <- function(expr, env = parent.frame(), quoted = FALSE) {
    if (!quoted) { expr <- substitute(expr) } # force quoted
    htmlwidgets::shinyRenderWidget(expr, "cyBasicMapOutput", env, quoted = TRUE)
  }
}



ORAEnrichment <- function(deGenes,universe, gsets, logFC, Dir="Both"){
  deGenes = deGenes[which(deGenes %in% universe)]
  tmp = rep(NA, length(gsets))
  #ora.stats = data.frame(p.value=tmp, p.adj = tmp, DeGeneNum=tmp,UpGene= tmp, DownGene=tmp, SetNum = tmp, N_q=tmp, SetNumAll=tmp, DeGene_in_Set=tmp)
  ora.stats = data.frame(p.value=tmp, p.adj = tmp, DeGeneNum=tmp,DE_UpGene= tmp, DE_DownGene=tmp, SetNum = tmp, N_q=tmp, Fold_Enrich=tmp,
                         SetNumAll=tmp, Total_DEG=tmp, Total_Gene=tmp,  DeGene_in_Set=tmp)
  
  totalDE = length(deGenes)
  if (Dir=="Up") {totalDE = sum(logFC > 0)
  } else if (Dir=="Down") {totalDE = sum(logFC < 0) }
  n = length(universe) - totalDE
  
  for (j in 1:length(gsets)){
    gset = gsets[[j]]
    DEinS = intersect(gset, deGenes)
    logFCinS = logFC[DEinS]
    totalDEinS = length(intersect(gset, deGenes))
    totalSinUniverse = length(intersect(gset, universe))
    
    N_q=totalDEinS- 0.5
    if (Dir=="Up") {N_q=length(logFCinS[logFCinS > 0])-0.5; DEinS<-names(logFCinS[logFCinS > 0])
    } else if (Dir=="Down") {N_q=length(logFCinS[logFCinS < 0])-0.5; DEinS<-names(logFCinS[logFCinS < 0]) }
    
    ora.stats[j, "p.value"] = phyper(q = N_q, m=totalDE, n = n, k = totalSinUniverse, lower.tail = FALSE)
    ora.stats[j, "DeGeneNum"] = totalDEinS
    ora.stats[j, "SetNum"] = totalSinUniverse  #previous versions used length(gset)
    ora.stats[j, "DE_UpGene"] = length(logFCinS[logFCinS > 0])
    ora.stats[j, "DE_DownGene"] = length(logFCinS[logFCinS < 0])
    ora.stats[j, "N_q"]=N_q
    ora.stats[j, "SetNumAll"]=length(gset)
    ora.stats[j, "Fold_Enrich"]= round( (totalDEinS/totalDE) / (totalSinUniverse/length(universe) )*100)/100
    ora.stats[j, "Total_DEG"]=totalDE
    ora.stats[j, "Total_Gene"]=length(universe)
    ora.stats[j, "DeGene_in_Set"]=paste(DEinS, collapse = ",")
    
  }
  ora.stats[, "p.adj"] = p.adjust(ora.stats[, "p.value"], method = "BH")
  ora.stats[, "DE_Direction"]=Dir
  ora.stats<-ora.stats%>%mutate(GeneSet= names(gsets))%>%
    arrange(p.value, dplyr::desc(N_q))%>% rownames_to_column('rank')%>%dplyr::select(-N_q)%>%relocate(GeneSet)
  
  return(ora.stats)
}


geneset_ui <- function(id) {
  ns <- shiny::NS(id)
  fluidRow(
    rclipboard::rclipboardSetup(),
    column(3,
           wellPanel(
             uiOutput(ns('loadedprojects')),
             conditionalPanel(ns = ns, "input.geneset_tabset=='Over-Representation Analysis (ORA)'",
                              radioButtons(ns("ORA_input_type"),label="Genes Used for ORA", choices=c("Subset from Comparison","Gene List"),inline = TRUE, selected="Subset from Comparison"),
                              conditionalPanel(ns = ns, "input.ORA_input_type=='Gene List'",
                                               textAreaInput(ns("ORA_list"), "Enter Gene List", "", cols = 5, rows=6),
                                               radioButtons(ns("ORA_universe"),label="Background Genes", choices=c("Genes in current project","Genes from Gene Sets"),inline = TRUE, selected="Genes in current project"))
             ),
             conditionalPanel(ns = ns, "input.ORA_input_type!='Gene List' || input.geneset_tabset!='Over-Representation Analysis (ORA)' ",
                              selectInput(ns("geneset_test"), label="Select Comparison for Gene Set Analysis", choices=NULL),
                              conditionalPanel(ns = ns, "input.geneset_tabset!='MetaBase Pathway View' && input.geneset_tabset!='Wikipathways View'",
                                               uiOutput(ns('Second_Comparison'))),
                              conditionalPanel(ns = ns, "input.geneset_tabset=='Over-Representation Analysis (ORA)' || input.geneset_tabset=='Gene Expression'",
                                               column(width=6,numericInput(ns("geneset_FCcut"), label= "Fold Change Cutoff", value = 1.2, min=1, step=0.1)),
                                               column(width=6,numericInput(ns("geneset_pvalcut"), label= "P Value Cutoff", value=0.01, min=0, step=0.001)),
                                               radioButtons(ns("geneset_psel"), label= "P value or P.adj Value?", choices= c("Pval"="Pval","Padj"="Padj"),inline = TRUE),
                                               radioButtons(ns("geneset_direction"), label= "Up- or Down-Regulated Genes?", choices= c("Both", "Up", "Down"),inline = TRUE),
                                               span(textOutput(ns("filteredgene1")), style = "color:red; font-size:15px; font-family:arial; font-style:italic"),
                                               span(textOutput(ns("filteredgene2")), style = "color:red; font-size:15px; font-family:arial; font-style:italic")
                              )),
             radioButtons(ns("MSigDB_species"), label= "Gene Set Species",choices=c("human","mouse", "rat"), inline = TRUE, selected = "human"),
             span(textOutput(ns("SpeciesGSEAInfo")), style = "color:red; font-size:13px; font-family:arial; font-style:italic"),
             selectInput(ns("map_genes"), "Gene Symbol Mapping", choices = c("Homologous Genes", "Change to UPPER case (human)", 
                                                                              "Change to Title Case (mouse/rat)", "No Change (as it is)"), selected="Homologous Genes"),
             conditionalPanel(ns = ns, "input.geneset_tabset=='Gene Set Enrichment Analysis (GSEA)' || input.geneset_tabset=='Over-Representation Analysis (ORA)' || input.geneset_tabset=='Gene Expression'",
                              conditionalPanel(ns = ns, "input.geneset_tabset=='Gene Set Enrichment Analysis (GSEA)'",
                                               radioButtons(ns("gene_rank_method"), label= "Rank Genes By",choices=c("logFC","-log(P-Value)"), inline = TRUE, selected = "logFC"),
                                               column(width=6,numericInput(ns("gsetMin"), label= "GeneSet Min Size",  value = 15, min=5, step=1)),
                                               column(width=6,numericInput(ns("gsetMax"), label= "GeneSet Max Size",  value = 1000, min=100, step=1)),
                                               sliderInput(ns("gsea_FDR"), "GSEA Adjusted P-Value Cutoff", min = 0, max = 1, step = 0.01, value = 0.25),
                                               checkboxInput(ns("gsea_collapase"), "Collapse Gene Sets",  FALSE, width="90%"),
                                               ),
                              conditionalPanel(ns = ns, "input.geneset_tabset=='Over-Representation Analysis (ORA)'",
                                               sliderInput(ns("ora_pvalue"), "ORA p.adj Cutoff:", min = 0, max = 1, step = 0.01, value = 0.05),
                                               checkboxInput(ns("ora_collapase"), "Collapse Gene Sets",  FALSE, width="90%"),
                                               ),
                              conditionalPanel(ns = ns, "input.MSigDB_species=='human'",
                                               checkboxGroupInput(ns("MSigDB_species_human_GSEA"), label= "Human Collections",
                                                                  choices= NULL, selected = NULL) ),
                              conditionalPanel(ns = ns, "input.MSigDB_species=='mouse'",
                                               checkboxGroupInput(ns("MSigDB_species_mouse_GSEA"), label= "Mouse Collections",
                                                                  choices= NULL, selected = NULL) ),
                              conditionalPanel(ns = ns, "input.MSigDB_species=='rat'",
                                               checkboxGroupInput(ns("MSigDB_species_rat_GSEA"), label= "Rat Collections",
                                                                  choices= NULL, selected = NULL) ),
                              checkboxInput(ns("use_customset"), "Add Custom GeneSet",  FALSE, width="90%"),
                              conditionalPanel(ns = ns, "input.use_customset==1",
                                               radioButtons(ns("custom_set_option"), label= "Custom GeneSet From",choices=c("Text Box","File"), inline = TRUE, selected = "Text Box"),
                                               conditionalPanel(ns = ns, "input.custom_set_option=='Text Box'",
                                                 textAreaInput(ns("CustomGeneset"), "Enter gene symbols", "", cols = 5, rows=6) ),
                                               conditionalPanel(ns = ns, "input.custom_set_option=='File'",
                                                                p("Prepare your own geneset in GMT format. One geneset per row;  each row has name, description, and the gene symbols in the gene set separated by tabs. Description is ignored, you can just put NA for it.", 
                                                                  style = "color:red; font-size:11px; font-family:arial; font-style:italic"),
                                                                tags$a(href="CustomGeneSets.gmt", "Download example custom GMT file with three human gene sets"),
                                                                fileInput(ns("custom_gmt_file"), "Choose geneset GMT file") )

                              )
             ),
             conditionalPanel(ns = ns, "input.geneset_tabset=='Gene Set Heatmap'",
                              radioButtons(ns("gs_heatmap_DE"),label="From the Gene Set, Show",inline = FALSE, choices=c("All Genes", "DE/LeadingEdge Genes Only"), 
                                           selected="All Genes"),
                              column(width=6,sliderInput(ns("hxfontsize_gsh"), "Column Font Size:", min = 2, max = 24, step = 1, value = 12)),
                              column(width=6,sliderInput(ns("hyfontsize_gsh"), "Row Font Size:", min = 2, max = 24, step = 1, value = 10)),
                              column(width=6,sliderInput(ns("htfontsize_gsh"), "Title Font Size:", min = 10, max = 32, step = 1, value = 14)),
                              column(width=6,sliderInput(ns("hlfontsize_gsh"), "Legend Font Size:", min = 2, max = 20, step = 1, value = 10)),
                              radioButtons(ns("gs_heatmap_label"),label="Gene Label",inline = TRUE, choices=c("UniqueID", "Gene.Name"), selected="Gene.Name"),
                              sliderInput(ns("geneset_heatmap_width"), "Heatmap Width:", min = 200, max = 3000, step = 50, value = 900),
                              sliderInput(ns("geneset_heatmap_height"), "Heatmap Height:", min = 200, max = 3000, step = 50, value = 800)
             ),
             
             conditionalPanel(ns = ns, "input.geneset_tabset=='KEGG Pathway View' || input.geneset_tabset=='MetaBase Pathway View' ",
                              radioButtons(ns("kegg_more_tests"), label= "Add more comparisons?", choices= c("Yes", "No"),selected="No", inline = TRUE),
                              conditionalPanel(ns = ns, "input.kegg_more_tests=='Yes'",
                                               selectInput(ns("geneset_test2"), label="2nd Comparison", choices=NULL),
                                               selectInput(ns("geneset_test3"), label="3rd Comparison", choices=NULL),
                                               selectInput(ns("geneset_test4"), label="4th Comparison", choices=NULL),
                                               selectInput(ns("geneset_test5"), label="5th Comparison", choices=NULL)),
                              conditionalPanel(ns = ns, "input.geneset_tabset=='KEGG Pathway View'",
                                 selectInput(ns("kegg_logFC"), label= "Gene log2FC Range:", choices= c(0.5, 1, 2, 3), selected=1),
                                 selectInput(ns("kegg_logFC_cpd"), label= "Compound log2FC Range:", choices= c(0.5, 1, 2, 3), selected=1),
                                radioButtons(ns("kegg_mapsample"), label= "Map Symbols to KEGG Nodes?", choices= c("Yes"=TRUE, "No"=FALSE),inline = TRUE)),
                              conditionalPanel(ns = ns, "input.geneset_tabset=='MetaBase Pathway View'",
                                               selectInput(ns("obj_style"), label="Network Object Style", choices=c("icon", "polygon", "none"), selected="icon")),
                              )
           )
    ),
    column(9,
           tabsetPanel(id=ns("geneset_tabset"),
                       tabPanel(title="Gene Set Enrichment Analysis (GSEA)",
                                actionButton(ns("compute_gsea"), "Compute/Refresh", style="color: #0961E3; background-color: #F6E98C ; border-color: #2e6da4"),
                                tags$br(),tags$hr(),
                                DT::dataTableOutput(ns("MSigDB_GSEA"))),
                       tabPanel(title="Over-Representation Analysis (ORA)",
                                actionButton(ns("compute_ora"), "Compute/Refresh", style="color: #0961E3; background-color: #F6E98C ; border-color: #2e6da4"),
                                tags$br(),tags$hr(),
                                conditionalPanel(ns = ns, "input.ORA_input_type!='Gene List'",
                                                 DT::dataTableOutput(ns("MSigDB_ORA"))),
                                conditionalPanel(ns = ns, "input.ORA_input_type=='Gene List'",
                                                 DT::dataTableOutput(ns("MSigDB_ORA_list"))) ),					
                       # need to use label
                       tabPanel(title="Gene Expression", fluidRow( column(4,textInput(ns('analysis_type_1'), 'Analysis Type')%>%disabled()),
                                                                   column(4,textInput(ns('x1'), 'Gene Set')) ),
                                conditionalPanel(ns = ns, "input.analysis_type_1=='GSEA'",  uiOutput(ns('ui.GSEA_GEX'))),
                                conditionalPanel(ns = ns, "input.analysis_type_1=='ORA'",  uiOutput(ns('ui.ORA_GEX')) ),
                                DT::dataTableOutput(ns("Expression"))),
                       tabPanel(title="Gene Set Heatmap", fluidRow( column(4,textInput(ns('analysis_type_2'), 'Analysis Type')%>%disabled()),
                                                                    column(4,textInput(ns('x2'), 'Gene Set')) ),
                                actionButton(ns("genesetheatmap"), "Save to output"),
                                uiOutput(ns("plot.geneset.heatmap"))),
                       tabPanel(title="KEGG Pathway View",
                              p("Select a KEGG Pathway by either clicking its name from the results table in the GSEA/ORA tab, or choose/search from the dropdown list below."),
                                selectizeInput(ns("sel_kegg_set"), label="KEGG Pathway for Visualization", choices = NULL, multiple = FALSE, width="600px", 
                                               options = list(placeholder =	'Type to search')),
                                actionButton(ns("keggSave"), "Save to output"),plotOutput(ns('keggView'))),
                       tabPanel(title="MetaBase Pathway View",
                                p("Select a MetaBase Pathway by either clicking its name from the results table in the GSEA/ORA tab, or choose/search from the dropdown list below."),
                                selectizeInput(ns("sel_metabase_set"), label="MetaBase Pathway for Visualization", choices = NULL, multiple = FALSE, width="600px", 
                                               options = list(placeholder =	'Type to search')),
                                #actionButton(ns("metabaseSave"), "Save to output"),
                                uiOutput(ns("plot.metabase"))),
                       tabPanel(title="Wikipathways View",
                                p("Select a wikipathway by either clicking its name from the results table in the GSEA/ORA tab, or choose/search from the dropdown list below."),
                                selectizeInput(ns("sel_wikipathways_set"), label="Wikipathways for Visualization", choices = NULL, multiple = FALSE, width="600px", 
                                               options = list(placeholder =	'Type to search')),
                                #actionButton(ns("metabaseSave"), "Save to output"),
                                svgPanZoomOutput(ns("wikipathways_plot"),width = "100%", height = "800px") ),
                       tabPanel(title="Help", htmlOutput('help_geneset'))
           )
    )
  )
}

geneset_server <- function(id) {
  shiny::moduleServer(id,
                      function(input, output, session) {
                        ns <- shiny::NS(id)
                        
                        #Use metabaser only when the gmt file exists
                        if (!file.exists("db/human/metabase_maps_genesymbols.gmt")) {
                          removeTab(session=session, inputId = "geneset_tabset", target = "MetaBase Pathway View")
                        }
                        if (!file.exists("db/human/kegg.gmt")) {
                          removeTab(session=session, inputId = "geneset_tabset", target = "KEGG Pathway View")
                        }
                        system="xOmicsShiny"
                        if (exists("system_info")){
                          if (!is.null(system_info)) {
                            if (system_info=="quickomics" ) { system="QuickOmics"}
                            if (system_info=="rnaview" ) { system="RNAView"}
                          }
                        }
                        #cat(system_info, date(), "\n")

                        ###Modify UI and input data based on QuickOmics or xOmicsShiny system
                        if (system %in% c("QuickOmics", "RNAView")) {
                          #shinyjs::hide(id = "comparison_2")
                          working_project=reactiveVal()
                          observe({
                            req(ProjectInfo)
                            working_project(ProjectInfo$ProjectID)
                          })
                          #browser()
                        } else if (system=="xOmicsShiny") {
                          output$loadedprojects <- renderUI({
                            req(length(working_project()) > 0)
                            radioButtons(ns("current_dataset"), label = "Change Working Dataset", choices=DS_names(), inline = F, selected=working_project())
                          })
                          
                          output$Second_Comparison<-renderUI({
                            tagList(
                              tags$hr(style="border-color: RoyalBlue;"),
                              checkboxInput(ns("comparison_2"), "Add 2nd Comparison for Gene + Compound Analysis",  FALSE, width="95%"),
                              conditionalPanel(ns = ns, "input.comparison_2==1",
                                               p("Select this option only when you have genes in one comparison, and compounds from the same expriment in another comparison. Compounds should have Gene.Name in KEGG compound ID format(e.g. C12345, C05519).", 
                                                 style = "color:red; font-size:11px; font-family:arial; font-style:italic"),
                                               selectInput(ns("dataset_2"), label="Dataset", choices=c("None", DS_names()), selected="None"),
                                               selectInput(ns("geneset_test_2nd"), label="Select 2nd Comparison", choices=NULL)
                              ),
                              tags$hr(style="border-color: RoyalBlue;"),
                            )
                          })
                          

                          observeEvent(input$dataset_2, {
                            req(input$dataset_2)
                            if (input$dataset_2!="None"){
                              DataIn_2 = DataInSets[[input$dataset_2]]
                              tests_2 = DataIn_2$tests
                              updateSelectizeInput(session,'geneset_test_2nd',choices=tests_2, selected=tests_2[1])
                            }
                          })
                          
                          observeEvent(input$current_dataset, {
                            working_project(input$current_dataset)
                          })
                          
                          gsea_control<-reactiveVal(0)
                          ora_control<-reactiveVal(0)
                          group_order <- reactiveVal()
                          ProjectInfo<-reactiveValues(ProjectID=NULL, Name=NULL, Species=NULL, ShortName=NULL, file1=NULL, file2=NULL, Path=NULL)
                          observeEvent(working_project(), {
                            DataIn=DataInSets[[working_project()]]
                            ProjectInfo$ProjectID=DataIn$ProjectID
                            ProjectInfo$Species=DataIn$Species
                            group_order(DataIn$group_order)
                          })
                          
                          DataReactive <-reactive({
                            req(length(working_project()) > 0)
                            req(DataInSets[[working_project()]])
                            Data1<-DataInSets[[working_project()]]
                            # append 2nd comparison if selected
                            if (!is.null(input$comparison_2)){
                              if (input$comparison_2==1) {
                                if (!is.null(input$dataset_2)) {
                                  if (input$dataset_2!="None" && input$geneset_test_2nd!=""){
                                    Data2<-DataInSets[[input$dataset_2]]
                                    if (!is.null(Data2)){
                                      #browser()
                                      ProteinGeneName1<-Data1$ProteinGeneName%>%dplyr::select(UniqueID, Gene.Name, Protein.ID)
                                      ProteinGeneName2<-Data2$ProteinGeneName%>%dplyr::select(UniqueID, Gene.Name, Protein.ID)
                                      ProteinGeneName_combined<-rbind(ProteinGeneName1, ProteinGeneName2)%>%dplyr::filter(!duplicated(UniqueID))
                                      Data1$ProteinGeneName=ProteinGeneName_combined
                                      #for 2nd comparison, only use new uniqueIDs
                                      results_long1_selTest<-Data1$results_long%>%dplyr::filter(test==input$geneset_test)
                                      results_long2_selTest<-Data2$results_long%>%dplyr::filter(test==input$geneset_test_2nd, !(UniqueID %in% results_long1_selTest$UniqueID) )%>%
                                        mutate(test=input$geneset_test)%>%dplyr::select(UniqueID, Gene.Name, test, logFC, P.Value, Adj.P.Value)
                                      results_long2_otherTests<-Data2$results_long%>%dplyr::select(UniqueID, Gene.Name, test, logFC, P.Value, Adj.P.Value)%>%
                                        dplyr::filter(test!=input$geneset_test_2nd,  !(UniqueID %in% ProteinGeneName1$UniqueID), 
                                                      test %in% unique(Data1$results_long$test) ) #add other tests if test name match
                                      #cat("Add", nrow(results_long2_selTest), "for", input$geneset_test_2nd, " |all other tests", nrow(results_long2_otherTests), "\n") #debug
                                      if (nrow(results_long2_selTest)==0) {results_long2_selTest=NULL}
                                      if (nrow(results_long2_otherTests)==0) {results_long2_otherTests=NULL}
                                      results_long_combined<-rbind(Data1$results_long%>%dplyr::select(UniqueID, Gene.Name, test, logFC, P.Value, Adj.P.Value), results_long2_selTest, results_long2_otherTests)
                                      Data1$results_long=results_long_combined
                                      #combined data_long
                                      data_long1<-Data1$data_long
                                      data_long2<-Data2$data_long%>%dplyr::filter(group %in%data_long1$group,  !(UniqueID %in% results_long1_selTest$UniqueID) )
                                      data_long_combined<-rbind(data_long1, data_long2)
                                      Data1$data_long=data_long_combined
                                    }
                                    
                                  }

                                }
                              }
                            } 
                            return(Data1)
                          })
                        }
                        ###
                        
                        #active_tests<-reactiveVal(NULL)
                        observeEvent(working_project(), {
                          req(working_project())
                          req(DataReactive())
                          DataIn = DataReactive()
                          tests=DataIn$tests
                          updateSelectizeInput(session,'geneset_test',choices=tests, selected=tests[1])
                         # active_tests(tests)
                          #cat("Update tests.", working_project(),  tests, "\n")
                          tests_more=c("None", tests)
                          updateSelectizeInput(session,'geneset_test2',choices=tests_more, selected="None")
                          updateSelectizeInput(session,'geneset_test3',choices=tests_more, selected="None")
                          updateSelectizeInput(session,'geneset_test4',choices=tests_more, selected="None")
                          updateSelectizeInput(session,'geneset_test5',choices=tests_more, selected="None")
                        })

                        
                        observe({
                          if (!is.null(gmt_file_info)) {  #update gmt file choices
                            gmt_info <- read.csv(gmt_file_info, fileEncoding="UTF-8-BOM")
                            gmt_use<-gmt_info%>%filter(Show=="YES")%>%mutate(label=str_c(Short_name, " ", label_name, " (", N_sets, ")"))
                            gmt_choice=gmt_use$gmt_file_name
                            names(gmt_choice)=gmt_use$label
                            selH=(gmt_use$Species=="human")
                            updateCheckboxGroupInput(session, "MSigDB_species_human_GSEA", choices = gmt_choice[selH], selected = gmt_choice[selH][1] )
                            selM=(gmt_use$Species=="mouse")
                            updateCheckboxGroupInput(session, "MSigDB_species_mouse_GSEA", choices = gmt_choice[selM], selected = gmt_choice[selM][1])
                            selR=(gmt_use$Species=="rat")
                            updateCheckboxGroupInput(session, "MSigDB_species_rat_GSEA", choices = gmt_choice[selR], selected = gmt_choice[selR][1])
                          }
                        })
                        #browser() 
                        observeEvent(working_project(),{
                          req(ProjectInfo)
                          if (!is.null(ProjectInfo$Species)) {
                            prj_species=tolower(str_trim(ProjectInfo$Species))
                            if (prj_species %in% c("human","mouse", "rat") ) {
                              updateRadioButtons(session, "MSigDB_species", selected = prj_species)
                             # cat("Speceis choice updatef, species is", ProjectInfo$Species, "\n")
                            } else if (prj_species %in% c("cho","chinese hamster", "mus musculus") ) {
                              updateRadioButtons(session, "MSigDB_species", selected = "mouse")
                              # cat("Speceis choice updatef, species is", ProjectInfo$Species, "\n")
                            } 
                          }
                        })
                        
                        observeEvent(c(input$MSigDB_species, ProjectInfo$Species), {
                          req(input$MSigDB_species, ProjectInfo$Species)
                          if (!is.null(ProjectInfo$Species)) {
                            species_info=str_c("Current project is for ",ProjectInfo$Species)
                            if (ProjectInfo$Species!=input$MSigDB_species) {
                              updateSelectizeInput(session,'map_genes', selected="Homologous Genes")
                              species_info=str_c(species_info, ", but gene set is from ", input$MSigDB_species, ". Choose how to map to ", 
                                    input$MSigDB_species, " genes.")
                            } else {
                              updateSelectizeInput(session,'map_genes', selected="No Change (as it is)")
                            }
                            output$SpeciesGSEAInfo <- renderText({species_info})
                          }
                          #Update MetaBase Pathway View tab
                          if (file.exists("db/human/metabase_maps_genesymbols.gmt")) {
                            if (input$MSigDB_species!="human") {
                              hideTab(session=session, inputId = "geneset_tabset", target = "MetaBase Pathway View")
                            } else {
                              showTab(session=session, inputId = "geneset_tabset", target = "MetaBase Pathway View")
                            }
                          }
                          if (system=="RNAView") {
                            hideElement(id="kegg_logFC_cpd")
                          }
                        })
                        
                        observeEvent(input$MSigDB_species, {
                          kegg_file=str_c("db/", input$MSigDB_species, "/kegg.gmt")
                          if (file.exists(kegg_file)) {
                            path1<-gmtPathways(kegg_file)
                            kegg_names=names(path1)
                            kegg_choices= c('Type to Search' = '', kegg_names)
                            #browser()
                            updateSelectizeInput(session, "sel_kegg_set", choices = kegg_choices, selected="Type to Search", server = TRUE)
                          }
                        })
                        
                        observe({
                          metabase_file="db/human/metabase_maps_genesymbols.gmt"
                          if (file.exists(metabase_file)) {
                            path1<-gmtPathways(metabase_file)
                            path_names=names(path1)
                            path_choices= c('Type to Search' = '', path_names)
                            updateSelectizeInput(session, "sel_metabase_set", choices = path_choices, selected="Type to Search", server = TRUE)
                          }
                        })
                        
                        observeEvent(input$MSigDB_species, {
                          wiki_file=str_c("db/", input$MSigDB_species, "/Wikipathways_Symbol.gmt")
                          if (file.exists(wiki_file)) {
                            path1<-gmtPathways(wiki_file)
                            wiki_names=names(path1)
                            wiki_choices= c('Type to Search' = '', wiki_names)
                            #browser()
                            updateSelectizeInput(session, "sel_wikipathways_set", choices = wiki_choices, selected="Type to Search", server = TRUE)
                          }
                        })
                        
                        observeEvent(c(input$ORA_input_type, input$geneset_tabset),{
                          req(input$ORA_input_type, input$geneset_tabset)
                          if (input$ORA_input_type=='Gene List' && input$geneset_tabset=='Over-Representation Analysis (ORA)')  {
                            hideTab(inputId="geneset_tabset", target="Gene Expression")
                            hideTab(inputId="geneset_tabset", target="Gene Set Heatmap")
                            if (file.exists("db/human/kegg.gmt")) { 
                              hideTab(inputId="geneset_tabset", target="KEGG Pathway View")
                            }
                          } else if (input$ORA_input_type!='Gene List' || input$geneset_tabset=='Gene Set Enrichment Analysis (GSEA)') {
                            showTab(inputId="geneset_tabset", target="Gene Expression")
                            showTab(inputId="geneset_tabset", target="Gene Set Heatmap")
                            if (file.exists("db/human/kegg.gmt")) { 
                              showTab(inputId="geneset_tabset", target="KEGG Pathway View")
                            }
                          }
                        } ) 
                        
                        ##################################### GSEA/ORA
                        
                        DataGenesetReactive_GSEA <- reactive({
                          DataIn = DataReactive()
                          results_long = DataIn$results_long
                          ProteinGeneName = DataIn$ProteinGeneName
                          
                          comp_sel = input$geneset_test
                          
                          
                          GSEAgene <- results_long %>% dplyr::filter(test==comp_sel, !is.na(Gene.Name), Gene.Name!="", !is.na(logFC) ) %>% 
                            dplyr::select(one_of(c("Gene.Name","logFC", "P.Value", "Adj.P.Value", "UniqueID", "test"))) %>%
                            dplyr::distinct(., Gene.Name,.keep_all = TRUE)
                          #browser()#bebug
                          
                          GSEA.terminals.df <- GSEAgene %>% arrange(dplyr::desc(logFC))
                          if (input$map_genes!="No Change (as it is)" && !(ProjectInfo$Species==input$MSigDB_species && input$map_genes=="Homologous Genes") ) {
                            if ( input$map_genes=="Change to UPPER case (human)" )  mapped_symbols<-toupper(GSEA.terminals.df$Gene.Name)
                            if ( input$map_genes=="Change to Title Case (mouse/rat)" )  mapped_symbols<-str_to_title(GSEA.terminals.df$Gene.Name)
                            if (ProjectInfo$Species!=input$MSigDB_species && input$map_genes=="Homologous Genes" ) {
                              mapped_symbols<-homolog_mapping(GSEA.terminals.df$Gene.Name, ProjectInfo$Species, input$MSigDB_species, homologs) }	  
                            GSEA.terminals.df<-GSEA.terminals.df%>%mutate(Gene.Name.Ori=Gene.Name, Gene.Name=mapped_symbols)%>% 
                              dplyr::distinct(., Gene.Name,.keep_all = TRUE)%>%dplyr::filter(!is.na(Gene.Name), Gene.Name!="") 
                          }
                          
                          if (input$gene_rank_method=="logFC") {
                            rank_list <- GSEA.terminals.df$logFC
                          } else { #by -log(P-value)
                            GSEA.terminals.df <-GSEA.terminals.df %>% mutate(minus_logPvalue=ifelse(logFC<0, log(P.Value), 0-log(P.Value)), minus_logPvalue=round(minus_logPvalue*100)/100 )%>%arrange( dplyr::desc(minus_logPvalue))
                            rank_list <- GSEA.terminals.df$minus_logPvalue
                          }
                          names(rank_list) <- GSEA.terminals.df$Gene.Name
                          
                          #browser()#bebug
                          
                          return(list("gene_list" = rank_list, "GSEA.terminals.df" = GSEA.terminals.df ))
                        }
                        )
                        
                        gsets_Reactive <- reactive ({
                          if (input$MSigDB_species == "human")  {
                            selectedGMT=input$MSigDB_species_human_GSEA
                          } else if (input$MSigDB_species == "mouse") {
                            selectedGMT=input$MSigDB_species_mouse_GSEA
                          } else if (input$MSigDB_species == "rat") {
                            selectedGMT=input$MSigDB_species_rat_GSEA
                          }
                          gsets_GSEA=NULL
                          for (f in selectedGMT) {
                            path1<-gmtPathways(str_c("db/", f) ) #read individual GTM file
                            gsets_GSEA=c(gsets_GSEA, path1)
                          }
                          if (input$use_customset) {
                            if (input$custom_set_option=="Text Box"){
                              CustomGeneset <- input$CustomGeneset
                              if(grepl("\n", CustomGeneset)) {
                                CustomGeneset <-  stringr::str_split(CustomGeneset, "\n")[[1]]
                              } else if (grepl(",", CustomGeneset)) {
                                CustomGeneset <-  stringr::str_split(CustomGeneset, ",")[[1]]
                              }
                              CustomGeneset <- gsub(" ", "", CustomGeneset, fixed = TRUE)
                              CustomGeneset <- unique(CustomGeneset[!is.na(CustomGeneset)])
                              if (input$MSigDB_species == "human")  {CustomGeneset <- toupper(CustomGeneset)}
                              validate(need(length(CustomGeneset)>2, message = "Please enter at least three custom genes"))
                              gsets_GSEA=c(gsets_GSEA, list(CustomGeneset=CustomGeneset))
                            }
                            if (input$custom_set_option=="File"){
                              req(input$custom_gmt_file)
                              CustomGeneset=gmtPathways(input$custom_gmt_file$datapath)
                              #cat("loaded", length(CustomGeneset), "custom gene sets.\n")
                              #browser() 
                              gsets_GSEA=c(gsets_GSEA, CustomGeneset)
                            }

                            #browser() #debug
                          }
                          return(gsets_GSEA)
                        })
                        
                        observeEvent(input$compute_gsea, {
                          gsea_control(gsea_control()+1)
                        })
                        
                        gsea_results<- eventReactive (gsea_control(), { 
                          withProgress(message = 'Running GSEA, pleaes be patient...', value = 0, {
                            getresults <- DataGenesetReactive_GSEA()
                            gsets_GSEA<-gsets_Reactive()
                            validate(need(length(gsets_GSEA)>0,"Please select at least one gene set."))
                            #browser()#bebug
                            logFC_list <- 	getresults$gene_list
                            ############################
                            gsMin = input$gsetMin
                            gsMax = input$gsetMax
                            ############################
                            output <- fgseaMultilevel(gsets_GSEA,logFC_list,minSize=gsMin,maxSize=gsMax)
                            if (input$gsea_collapase) {
                              collapsed_set <- collapsePathways(output[order(pval)][padj <= input$gsea_FDR], gsets_GSEA,logFC_list )
                              output= output[pathway %in% collapsed_set$mainPathways]
                            }
                            
                            res = output %>%                   # Using dplyr functions
                              mutate_if(is.numeric,
                                        signif,
                                        digits = 3)  #some may fail for very small numbers (<1e-20), added formatSignif in datatable as well.
                            
                            res <- res %>% dplyr::rename(GeneSet=pathway)%>%
                              dplyr::select(-ES, -log2err) %>%  #keep leadingEdge
                              arrange(padj, dplyr::desc(abs(NES))) %>%dplyr::filter(padj<=input$gsea_FDR)%>%
                              rownames_to_column('rank') %>%
                              relocate(GeneSet)
                          } )
                          return(res)
                        } )
                        
                        output$MSigDB_GSEA <-  DT::renderDT(server=FALSE,{ withProgress(message = 'Processing...', value = 0, {
                          res<-gsea_results()
                          validate(need(nrow(res)>0,"No results. Try to increase p.adj cutoff, or try a different comparison."))
                          res$Action<-vapply(1:nrow(res), function(i){
                            as.character(
                              rclipButton(
                                paste0("clipbtn_", i), 
                                label = "Copy Leading Edge Genes", 
                                clipText = paste(unlist(res[i,  "leadingEdge"]), collapse=","), 
                                #icon = icon("clipboard"),
                                icon = icon("copy", lib = "glyphicon"),
                                class = "btn-primary btn-sm"
                              )
                            )
                          }, character(1L))
                          res<-res%>%dplyr::relocate(Action, .before=leadingEdge)
                          
                          DT::datatable(
                            res,  extensions = 'Buttons', escape = FALSE, selection = 'none', class = 'cell-border strip hover',
                            options = list(    dom = 'lBfrtip', pageLength = 15,
                                               buttons = list(
                                                 list(extend = "csv", text = "Download Page", filename = "Page_results",
                                                      exportOptions = list(modifier = list(page = "current"))),
                                                 list(extend = "csv", text = "Download All", filename = "All_Results",
                                                      exportOptions = list(modifier = list(page = "all")))
                                               )
                            )	)%>%formatSignif(columns=c('pval', "padj"), digits=3)   %>% formatStyle(1, cursor = 'pointer',color='blue')
                        })
                        })
                        
                        observeEvent(input$MSigDB_GSEA_cell_clicked, {
                          info = input$MSigDB_GSEA_cell_clicked
                          if (is.null(info$value) || info$col != 1) return()
                          analysis_type = 'GSEA'
                          updateTabsetPanel(session, 'geneset_tabset', selected = 'Gene Expression')
                          updateTextInput(session, 'x1', value = info$value)
                          updateTextInput(session, 'x2', value = info$value)
                          #updateTextInput(session, 'x3', value = info$value)
                          updateSelectizeInput(session, "sel_kegg_set",selected=info$value)
                          updateSelectizeInput(session, "sel_metabase_set",selected=info$value)
                          updateSelectizeInput(session, "sel_wikipathways_set",selected=info$value)
                          updateTextInput(session, 'analysis_type_1', value = analysis_type)
                          updateTextInput(session, 'analysis_type_2', value = analysis_type)
                         # updateTextInput(session, 'analysis_type_3', value = analysis_type)
                        })
                        
                        
                        
                        
                        ####################### ORA
                        DataGenesetReactive_ORA <- reactive({
                          req(DataReactive(), input$ORA_input_type, input$map_genes)
                          req( input$geneset_test %in% DataReactive()$tests)
                          if (input$ORA_input_type!='Gene List') {
                            DataIn = DataReactive()
                            results_long = DataIn$results_long
                            #ProteinGeneName = DataIn$ProteinGeneName
                            
                            comp_sel = input$geneset_test
                            
                            all_genes <-results_long %>% dplyr::filter(test == comp_sel) %>%dplyr::filter(!is.na(Gene.Name),  Gene.Name!="")%>%
                              dplyr::select(one_of(c("Gene.Name"))) %>% 
                              collect %>% .[["Gene.Name"]] %>% unique()
                            
                            absFCcut = log2(as.numeric(input$geneset_FCcut))
                            pvalcut = as.numeric(input$geneset_pvalcut)
                            terminals.df <-  results_long %>% dplyr::filter(test == comp_sel) %>%dplyr::filter(!is.na(Gene.Name),  Gene.Name!="")%>%
                              dplyr::select(one_of(c("Gene.Name","logFC","P.Value", "Adj.P.Value", "UniqueID", "test")))%>%arrange(P.Value)  %>%
                              dplyr::distinct(., Gene.Name,.keep_all = TRUE)
                            
                            #browser()#bebug
                            #terminals.df <-  filteredgene%>%arrange(P.Value)  #old way
                            if (input$map_genes!="No Change (as it is)" && !(ProjectInfo$Species==input$MSigDB_species && input$map_genes=="Homologous Genes") ) {
                              if (input$map_genes=="Change to UPPER case (human)")  {
                                mapped_symbols<-toupper(terminals.df$Gene.Name); all_genes=toupper(all_genes)
                              } else if ( input$map_genes=="Change to Title Case (mouse/rat)" ) {
                                mapped_symbols<-str_to_title(terminals.df$Gene.Name); all_genes=str_to_title(all_genes)
                              } else if (ProjectInfo$Species!=input$MSigDB_species && input$map_genes=="Homologous Genes" ) {
                                mapped_symbols<-homolog_mapping(terminals.df$Gene.Name, ProjectInfo$Species, input$MSigDB_species, homologs)
                                all_genes <-homolog_mapping(all_genes , ProjectInfo$Species, input$MSigDB_species, homologs)
                              } 
                              terminals.df<-terminals.df%>%mutate(Gene.Name.Ori=Gene.Name, Gene.Name=mapped_symbols)%>%
                                dplyr::distinct(., Gene.Name,.keep_all = TRUE)%>%dplyr::filter(!is.na(Gene.Name),  Gene.Name!="")
                              }
                            
                            filteredgene <-  terminals.df %>%  mutate(psel=input$geneset_psel, P.stat=ifelse(psel == "Padj",  Adj.P.Value,  P.Value)) %>%
                              dplyr::filter(abs(logFC) > absFCcut,  P.stat < pvalcut) %>%
                              dplyr::select(one_of(c("Gene.Name","logFC","P.Value", "Adj.P.Value", "UniqueID", "test")))
                            
                            sig_genes <- filteredgene$logFC
                            names(sig_genes) <- filteredgene$Gene.Name
                            sig_genes_Dir=sig_genes
                            if (input$geneset_direction=="Up") {
                              sig_genes_Dir=sig_genes[sig_genes>0]
                            } else if (input$geneset_direction=="Down") {
                              sig_genes_Dir=sig_genes[sig_genes<0]
                            }
                            #output$geneset_filteredgene <- renderText({ paste("Selected Genes:",length(sig_genes_Dir), " (",input$geneset_direction, " Direction)",  sep="")})
                            return(list("sig_genes" = sig_genes, "all_genes" = all_genes, "terminals.df"=terminals.df, "sig_genes_Dir" = sig_genes_Dir))
                          }
                        })
                        
                        
                        
                        DataGenesetReactive_ORA_list <- reactive({
                          if (input$ORA_input_type=='Gene List') {
                            ORA_list <- input$ORA_list
                            ORA_list=str_replace_all(ORA_list, " ", "")
                            if(grepl("\n",ORA_list)) {
                              ORA_list <-  stringr::str_split(ORA_list, "\n")[[1]]
                            } else if(grepl(",",ORA_list)) {
                              ORA_list <-  stringr::str_split(ORA_list, ",")[[1]]
                            }
                            ORA_list<-ORA_list [ORA_list !=""]
                            validate(need(length(ORA_list)>1, message = "Please input at least 2 valid genes."))
                              DataIn = DataReactive()
                              ProteinGeneName = DataIn$ProteinGeneName
                              all_genes <- dplyr::filter(ProteinGeneName, !is.na(Gene.Name), Gene.Name!="") %>%
                                dplyr::select(one_of(c("Gene.Name"))) %>% collect %>% .[["Gene.Name"]] %>% unique()
                            if (input$map_genes!="No Change (as it is)" && !(ProjectInfo$Species==input$MSigDB_species && input$map_genes=="Homologous Genes") ) {
                              if (input$map_genes=="Change to UPPER case (human)")  {
                                ORA_list<-toupper(ORA_list); all_genes=toupper(all_genes)
                              } else if ( input$map_genes=="Change to Title Case (mouse/rat)" ) {
                                ORA_list<-str_to_title(ORA_list); all_genes=str_to_title(all_genes)
                              } else if (ProjectInfo$Species!=input$MSigDB_species && input$map_genes=="Homologous Genes" ) {
                                ORA_list<-homolog_mapping(ORA_list, ProjectInfo$Species, input$MSigDB_species, homologs)
                                all_genes <-homolog_mapping(all_genes , ProjectInfo$Species, input$MSigDB_species, homologs)
                              } 
                            }
                            if (input$ORA_universe!="Genes in current project") {
                              gsets_GSEA<-gsets_Reactive()
                              all_genes<-unlist(gsets_GSEA)%>%unique()
                            }
                            sig_genes=rep(1, length(ORA_list))
                            names(sig_genes)=ORA_list
                            return(list("sig_genes" = sig_genes, "all_genes" = all_genes, "terminals.df"=NULL, "sig_genes_Dir" = NULL))
                          } 
                        })
                        observe({
                          req(working_project())
                          req(input$geneset_test)
                          req(DataReactive())
                          req(DataGenesetReactive_ORA())
                          req(DataGenesetReactive_ORA()$sig_genes)
                          tmpdat=DataGenesetReactive_ORA()$sig_genes
                          output$filteredgene1 <- renderText({paste("Genes Pass Cutoff (DEGs):",length(tmpdat),sep="")})
                          output$filteredgene2 <- renderText({paste("Genes Up: ", sum(tmpdat>0), "; Genes Down: ", sum(tmpdat<0),sep="")})
                          DEGs=names(tmpdat)
                          updateTextAreaInput(session, "ORA_list", value=paste(DEGs, collapse="\n"))
                        })

                        
                        observeEvent(input$compute_ora, {
                          ora_control(ora_control()+1)
                        })
                        
                        ora_results<- eventReactive (ora_control(), { 
                          withProgress(message = 'Running ORA...', value = 0, {
                            gsets_GSEA<-gsets_Reactive()
                            validate(need(length(gsets_GSEA)>0,"Please select at least one gene set."))
                            getresults <- DataGenesetReactive_ORA()
                            logFC_list <- 	getresults$sig_genes
                            all_genes <- 	getresults$all_genes
                            gsa <- ORAEnrichment (deGenes=names(logFC_list),universe=all_genes, gsets_GSEA, logFC =logFC_list, Dir=input$geneset_direction)
                            if (input$ora_collapase) {
                              Dir=input$geneset_direction
                              logFC_list1=logFC_list
                              if (Dir=="Up") {logFC_list1=logFC_list[logFC_list>0]
                              } else if (Dir=="Down") {logFC_list1=logFC_list[logFC_list<0]}
                              collapsed_set <- collapsePathwaysORA(foraRes=gsa%>%mutate(pathway=GeneSet)%>%filter(p.adj<= input$ora_pvalue),
                                                    pathways=gsets_GSEA,genes=names(logFC_list1), universe=all_genes, pval.threshold=0.05) 
                              #browser()
                              gsa<-gsa%>%dplyr::filter(GeneSet %in% collapsed_set$mainPathways)
                            }
                            ############################
                            #browser()
                            res <- 	gsa %>%dplyr::filter( p.adj <= input$ora_pvalue)                    # Using dplyr functions
                           #   mutate_if(is.numeric,signif,digits = 3)
                          })
                        })
                        
                        ora_results_list<- eventReactive (ora_control(), { 
                          withProgress(message = 'Running ORA...', value = 0, {
                            gsets_GSEA<-gsets_Reactive()
                            validate(need(length(gsets_GSEA)>0,"Please select at least one gene set."))
                            getresults <- DataGenesetReactive_ORA_list()
                            logFC_list <- 	getresults$sig_genes
                            all_genes <- 	getresults$all_genes
                            
                            gsa <- ORAEnrichment (deGenes=names(logFC_list),universe=all_genes, gsets_GSEA, logFC =logFC_list, Dir="Both")%>%
                              dplyr::select(-DE_UpGene, -DE_DownGene, -DE_Direction)
                            if (input$ora_collapase) {
                              collapsed_set <- collapsePathwaysORA(foraRes=gsa%>%mutate(pathway=GeneSet)%>%filter(p.adj<= input$ora_pvalue),
                                                                   pathways=gsets_GSEA,genes=names(logFC_list), universe=all_genes, pval.threshold=0.05) 
                              gsa<-gsa%>%dplyr::filter(GeneSet %in% collapsed_set$mainPathways)
                            }
                            ############################
                            res <- 	gsa %>%dplyr::filter( p.adj <= input$ora_pvalue)                   # Using dplyr functions
                             #mutate_if(is.numeric,signif,digits = 3)
                          })
                        })
                        
                        
                        output$MSigDB_ORA <-DT::renderDT(server=FALSE,{ withProgress(message = 'Processing...', value = 0, {
                          res<-ora_results()
                          validate(need(nrow(res)>0,"No results. Try to increase p.adj cutoff, or ajust DEG cutoff."))
                          res$Action<-vapply(1:nrow(res), function(i){
                            as.character(
                              rclipButton(
                                paste0("clipbtn_", i), 
                                label = "Copy DE Gene Names", 
                                clipText = res[i,  "DeGene_in_Set"], 
                                #icon = icon("clipboard"),
                                icon = icon("copy", lib = "glyphicon"),
                                class = "btn-primary btn-sm"
                              )
                            )
                          }, character(1L))
                          res<-res%>%dplyr::relocate(Action, .before=DeGene_in_Set)
                  
                          DT::datatable(
                            res,  extensions = 'Buttons', escape = FALSE, selection = 'none', class = 'cell-border strip hover',
                            options = list(    dom = 'lBfrtip', pageLength = 15,
                                               buttons = list(
                                                 list(extend = "csv", text = "Download Page", filename = "Page_results",
                                                      exportOptions = list(modifier = list(page = "current"))),
                                                 list(extend = "csv", text = "Download All", filename = "All_Results",
                                                      exportOptions = list(modifier = list(page = "all")))
                                               ))	)%>%formatSignif(columns=c('p.value', "p.adj"), digits=3)  %>% formatStyle(1, cursor = 'pointer',color='blue')
                        }) })
                        
                        output$MSigDB_ORA_list <-DT::renderDT(server=FALSE,{ withProgress(message = 'Processing...', value = 0, {
                          res<-ora_results_list()
                          validate(need(nrow(res)>0,"No results. Try to increase p.adj cutoff, or use a different list."))
                          res$Action<-vapply(1:nrow(res), function(i){
                            as.character(
                              rclipButton(
                                paste0("clipbtn_", i), 
                                label = "Copy DE Gene Names", 
                                clipText = res[i,  "DeGene_in_Set"], 
                                #icon = icon("clipboard"),
                                icon = icon("copy", lib = "glyphicon"),
                                class = "btn-primary btn-sm"
                              )
                            )
                          }, character(1L))
                          res<-res%>%dplyr::relocate(Action, .before=DeGene_in_Set)
                          DT::datatable(
                            res,  extensions = 'Buttons',  escape = FALSE, 
                            options = list(    dom = 'lBfrtip', pageLength = 15,
                                               buttons = list(
                                                 list(extend = "csv", text = "Download Page", filename = "Page_results",
                                                      exportOptions = list(modifier = list(page = "current"))),
                                                 list(extend = "csv", text = "Download All", filename = "All_Results",
                                                      exportOptions = list(modifier = list(page = "all")))
                                               ))	)%>%formatSignif(columns=c('p.value', "p.adj"), digits=3)
                        }) })
                        
                        observeEvent(input$MSigDB_ORA_cell_clicked, {
                          if (input$ORA_input_type!='Gene List') {
                            info = input$MSigDB_ORA_cell_clicked
                            if (is.null(info$value) || info$col != 1) return()
                            analysis_type = 'ORA'
                            updateTabsetPanel(session, 'geneset_tabset', selected = 'Gene Expression')
                            updateTextInput(session, 'x1', value = info$value)
                            updateTextInput(session, 'x2', value = info$value)
                            #updateTextInput(session, 'x3', value = info$value)
                            updateSelectizeInput(session, "sel_kegg_set",selected=info$value)
                            updateSelectizeInput(session, "sel_metabase_set",selected=info$value)
                            updateSelectizeInput(session, "sel_wikipathways_set",selected=info$value)
                            updateTextInput(session, 'analysis_type_1', value = analysis_type)
                            updateTextInput(session, 'analysis_type_2', value = analysis_type)
                            #updateTextInput(session, 'analysis_type_3', value = analysis_type)
                          }
                        })
                        
                        
                        output$ui.ORA_GEX <- renderUI({
                          tagList(
                            DT::dataTableOutput(ns('selected_ORA_set')),
                            h5("The table below contains all genes from the selected gene set. Selected DE genes for ORA are shown in green."))
                        })
                        
                        output$ui.GSEA_GEX <- renderUI({
                          tagList(
                            DT::dataTableOutput(ns('selected_GSEA_set')),
                            h5("GSEA Plot"),
                            plotOutput(ns('GSEA_plot')),
                            textOutput(ns("leadingEdge_genes")), 
                            h5("The table below contains all genes from the selected gene set. Leading edge genes are shown in green."))
                        })
                        
                        sel_GSEA_set<- reactive({
                          analysis_type = input$analysis_type_1
                          if (analysis_type=="GSEA"){
                            res<-gsea_results()
                            ID = input$x1
                            if (ID!=""){
                              sel_res<-res%>%filter(GeneSet==ID)
                              sel_resT<-sel_res%>%dplyr::select(-leadingEdge)
                              LE_Genes=sel_res$leadingEdge[[1]]
                              return(list(table=sel_resT, genes=LE_Genes))
                            } else {return(list(table=NULL, genes=NULL)) }
                          }  else {return(list(table=NULL, genes=NULL)) }
                        })
                        
                        output$selected_GSEA_set <- DT::renderDT(server=FALSE, {
                          validate(need(sel_GSEA_set()$table, "select a gene set from GSEA table")) 
                          DT::datatable(sel_GSEA_set()$table, rownames=F, options = list(dom="t", ordering=F))%>%formatSignif(columns=c('pval', "padj"), digits=3)
                        } )
                        
                        output$leadingEdge_genes <- renderText({
                          genes=sel_GSEA_set()$genes
                          paste("Leading Edge Genes (the core group of genes that accounts for the gene set's enrichment signal):", paste(genes, collapse = ", "))
                        })
                        
                        output$GSEA_plot = renderPlot({
                          getresults <- DataGenesetReactive_GSEA()
                          gsets_GSEA<-gsets_Reactive()
                          ID = input$x1
                          logFC_list <- 	getresults$gene_list
                          plotEnrichment(gsets_GSEA[[ID]],logFC_list) + labs(title=ID)
                        }, height=400) 
                        
                        sel_ORA_set<- reactive({
                          analysis_type = input$analysis_type_1
                          if (analysis_type=="ORA"){
                            res<-ora_results()
                            ID = input$x1
                            if (ID!=""){
                              sel_res<-res%>%filter(GeneSet==ID)
                              #cat("ORA selected", ID, "\n")
                              return(list(table=sel_res))
                            } else {return(list(table=NULL)) }
                          } else {return(list(table=NULL)) }
                        })
                        
                        output$selected_ORA_set <- DT::renderDT(server=FALSE, {
                          validate(need(sel_ORA_set()$table, "select a gene set from ORA table")) 
                          DT::datatable(sel_ORA_set()$table, rownames=F, options = list(dom="t", ordering=F))%>%formatSignif(columns=c('p.value', "p.adj"), digits=3)
                        } )
                        
                        keggView_out <- reactive({withProgress(message = 'Making KEGG Pathway View...', value = 0, {
                          #analysis_type = input$analysis_type_2
                          #ID = input$x3
                          ID=input$sel_kegg_set
                          validate(need(ID!="", message = "Please select a KEGG pathway to map logFC data to it."))
                          validate(need(str_detect(ID, "^(hsa|mmu|rno)\\d{5}"), message = "Only works on human/mouse/rat KEGG pathways."))
                          species=input$MSigDB_species
                          if (species=="rat") {species="rno"}
                          if (species=="mouse") {species="mmu"}
                          if (species=="human") {species="hsa"}
                          pid <- strsplit(ID,"_")[[1]][1]
                          img.file <- paste(pid,"pathview","png",sep=".")
                          dataIn=DataReactive()
                          results_long=dataIn$results_long
                          ProteinGeneName = dataIn$ProteinGeneName
                          comp_sel = input$geneset_test
                          
                          data_results=dataIn$data_results
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
                          ##Make FCdata  from results_long
                          res1<-results_long%>%dplyr::filter(test==comp_sel)%>%dplyr::select(UniqueID, logFC); names(res1)[2]=comp_sel
                          FC_df<-data.frame(UniqueID=unique(results_long$UniqueID))%>%left_join(ProteinGeneName%>%dplyr::select(UniqueID, Gene.Name))%>%
                            left_join(res1)
                          if (length(tests)>1) {
                            for (i in 2:length(tests)){
                              comp_add=tests[i]
                              res1<-results_long%>%dplyr::filter(test==comp_add)%>%dplyr::select(UniqueID, logFC); names(res1)[2]=comp_add
                              FC_df<-FC_df%>%left_join(res1)
                            }
                          }
                          
                          if (input$map_genes!="No Change (as it is)" && !(ProjectInfo$Species==input$MSigDB_species && input$map_genes=="Homologous Genes") ) {
                            if ( input$map_genes=="Change to UPPER case (human)" )  mapped_symbols<-toupper(FC_df$Gene.Name)
                            if ( input$map_genes=="Change to Title Case (mouse/rat)" )  mapped_symbols<-str_to_title(FC_df$Gene.Name)
                            if (ProjectInfo$Species!=input$MSigDB_species && input$map_genes=="Homologous Genes" ) {
                              mapped_symbols<-homolog_mapping(FC_df$Gene.Name, ProjectInfo$Species, input$MSigDB_species, homologs) }	  
                              FC_df<-FC_df%>%mutate(Gene.Name.Ori=Gene.Name, Gene.Name=mapped_symbols)
                          }
                          
                          selCol=which(names(FC_df) %in% tests)
                          sel_gene=which(!is.na(FC_df$Gene.Name))
                          FCdata=data.matrix(FC_df[sel_gene, selCol])
                          all.names=FC_df$Gene.Name[sel_gene]
                          rownames(FCdata)=all.names
                          if (ncol(FCdata)>1) {	img.file <- paste(pid,"pathview.multi.png",sep=".")}
                          #Separte genes and compounds
                          is.compound=str_detect(all.names, "C\\d{5}$")
                          FCdata_gene=FCdata[!is.compound, , drop=F]
                          FCdata_compound=FCdata[is.compound, , drop=F]
                          #browser()
                          if (sum(!is.compound)>0) {
                            tmp <- try (pathview(gene.data=FCdata_gene, cpd.data=FCdata_compound, pathway.id=pid, kegg.dir="./kegg", kegg.native = T, gene.idtype="SYMBOL", species=species,
                                                 low = list(gene = "green", cpd = "blue"), mid = list(gene = "gray90", cpd = "gray90"), high = list(gene = "red", cpd ="yellow"),
                                                 same.layer = F, map.symbol=as.logical(input$kegg_mapsample), limit=list(gene=as.numeric(input$kegg_logFC), cpd=as.numeric(input$kegg_logFC_cpd)) ) )
                            if (class(tmp) == "try-error") { #if no gene mapped, the above command will fail, try the one below with compound only
                              tmp <- pathview(cpd.data=FCdata_compound, pathway.id=pid, kegg.dir="./kegg", kegg.native = T,  species=species,
                                              low = list(gene = "green", cpd = "blue"), mid = list(gene = "gray90", cpd = "gray90"), high = list(gene = "red", cpd ="yellow"),
                                              same.layer = F, limit=list(gene=as.numeric(input$kegg_logFC), cpd=as.numeric(input$kegg_logFC_cpd)) )
                            
                            }
                          } else {
                            tmp <- pathview(cpd.data=FCdata_compound, pathway.id=pid, kegg.dir="./kegg", kegg.native = T,  species=species,
                                            low = list(gene = "green", cpd = "blue"), mid = list(gene = "gray90", cpd = "gray90"), high = list(gene = "red", cpd ="yellow"),
                                            same.layer = F, limit=list(gene=as.numeric(input$kegg_logFC), cpd=as.numeric(input$kegg_logFC_cpd)) )
                          }
                            
                          #browser() #debug
                          return(img.file)
                        })
                        })
                        output$wikipathways_plot <- renderSvgPanZoom({withProgress(message = 'Making WikiPathways View...', value = 0, {
                          ID=input$sel_wikipathways_set
                          validate(need(ID!="", message = "Please select a Wikipathway to map logFC data to it."))
                          validate(need(str_detect(ID, "WP\\d+$"), message = "Only works on human/mouse/rat KEGG pathways."))
                          wiki_ID=str_extract(ID, "WP\\d+$")
                          species=input$MSigDB_species
                          dataIn=DataReactive()
                          results_long=dataIn$results_long
                          ProteinGeneName = dataIn$ProteinGeneName
                          comp_sel = input$geneset_test
                          
                          ##Make FCdata  from results_long
                          FC_df<-results_long%>%dplyr::filter(test==comp_sel)%>%dplyr::select(UniqueID, logFC, Gene.Name)  
                          if (input$map_genes!="No Change (as it is)" && !(ProjectInfo$Species==input$MSigDB_species && input$map_genes=="Homologous Genes") ) {
                            if ( input$map_genes=="Change to UPPER case (human)" )  mapped_symbols<-toupper(FC_df$Gene.Name)
                            if ( input$map_genes=="Change to Title Case (mouse/rat)" )  mapped_symbols<-str_to_title(FC_df$Gene.Name)
                            if (ProjectInfo$Species!=input$MSigDB_species && input$map_genes=="Homologous Genes" ) {
                              mapped_symbols<-homolog_mapping(FC_df$Gene.Name, ProjectInfo$Species, input$MSigDB_species, homologs) }	  
                            FC_df<-FC_df%>%mutate(Gene.Name.Ori=Gene.Name, Gene.Name=mapped_symbols)%>%dplyr::filter(!is.na(Gene.Name), Gene.Name!="")
                          }
                          FCdata=FC_df$logFC; names(FCdata)=FC_df$Gene.Name
                          
                          cat(comp_sel, wiki_ID, length(FCdata), "\n")
                          p1 <- wpplot(wiki_ID)
                          #p2 <- p1 |> wp_bgfill(FCdata, low='darkgreen', high='firebrick', legend_x = .9, legend_y = .95) #|> is the base R "pipe" operator. It was new in version 4.1.0.
                          p2 <- wp_bgfill(p1, FCdata, low='darkgreen', high='firebrick', legend_x = .9, legend_y = .95)
                          svgPanZoom(paste(p2$svg, collapse="\n"),  controlIconsEnabled = T)  
                        } )
                        })
                        
                        ##MetabaseR map, only when human gene set is selected
                        Data_metabase <- reactive({
                          req(DataReactive())
                          if (input$sel_metabase_set!="") {
                            cat(input$sel_metabase_set, "\n") #for debug
                            #get data
                            dataIn=DataReactive()
                            results_long=dataIn$results_long
                            ProteinGeneName = dataIn$ProteinGeneName
                            if (input$map_genes!="No Change (as it is)" && !(ProjectInfo$Species==input$MSigDB_species && input$map_genes=="Homologous Genes") ) {
                              if ( input$map_genes=="Change to UPPER case (human)" )  mapped_symbols<-(ProteinGeneName$Gene.Name)
                              if ( input$map_genes=="Change to Title Case (mouse/rat)" )  mapped_symbols<-str_to_title(ProteinGeneName$Gene.Name)
                              if (ProjectInfo$Species!=input$MSigDB_species && input$map_genes=="Homologous Genes" ) {
                                mapped_symbols<-homolog_mapping(ProteinGeneName$Gene.Name, ProjectInfo$Species, input$MSigDB_species, homologs) }	  
                              ProteinGeneName<-ProteinGeneName%>%mutate(Gene.Name.Ori=Gene.Name, Gene.Name=mapped_symbols)
                            }
                            #gene entrezID
                            Name_list<-ProteinGeneName%>%dplyr::filter(!is.na(Gene.Name), Gene.Name!="", Gene.Name!="NA")%>%
                                dplyr::select(Gene.Name)%>%unlist%>%unname
                            Gene2ID=mapIds(org.Hs.eg.db, keys=Name_list, column="ENTREZID", keytype="SYMBOL")
                            #browser()
                            df_ID<-data.frame(Gene.Name=names(Gene2ID), EntrezID=Gene2ID)%>%filter(!is.na(Gene.Name), Gene.Name!="", 
                                                                                                   !duplicated(Gene.Name) )
                            ProteinGeneName1<-ProteinGeneName%>%dplyr::select(UniqueID, Gene.Name)%>%left_join(df_ID)
                            results_long<-results_long%>%left_join(ProteinGeneName1%>%dplyr::select(UniqueID, EntrezID))
                            #get logFC data from data_results
                            tests=input$geneset_test
                            if (input$kegg_more_tests=="Yes") {
                              if (input$geneset_test2!="None") {tests=c(tests, input$geneset_test2)}
                              if (input$geneset_test3!="None") {tests=c(tests, input$geneset_test3)}
                              if (input$geneset_test4!="None") {tests=c(tests, input$geneset_test4)}
                              if (input$geneset_test5!="None") {tests=c(tests, input$geneset_test5)}
                            }
                            data_plot<-list()
                            for (t in tests){
                              data1<-results_long%>%filter(test==t, !is.na(EntrezID))%>%dplyr::select(EntrezID, logFC, Adj.P.Value, UniqueID)
                              list1=list(data1); names(list1)[1]=t
                              data_plot=c(data_plot, list1)
                            }
                            return(data_plot)
                          }
                    
                        })
                        
                        output$plot.metabase=renderUI({
                          req(working_project())
                          req(input$geneset_test)
                          req(DataReactive())
                          if (input$sel_metabase_set!="") {
                            tmp<-try(view.map(input$sel_metabase_set)) #test if the metabaseR map can be loaded, some like "Oxidative phosphorylation" has issues
                            if (class(tmp)[1]=="try-error") {
                              cat(input$sel_metabase_set, "can't be plotted\n")
                                tagList(tags$div(
                                  tags$p("This Metabase pathway cannot be drawn. Please select another pathway.")
                                )) 
                            } else {
                              withProgress(message = 'Plotting metabase pathway...', value = 0, {
                                cat(input$sel_metabase_set, "\n") #for debug
                                suppressWarnings(metabase_map<-  view.map(input$sel_metabase_set,
                                                                          datasets = Data_metabase(),  bg_image=FALSE,  input.types="gene", nwobj_style=input$obj_style ) )
                                output$my_widget <- renderMetabase(metabase_map)
                                plot_height=metabase_map$height+50
                                plot_width=metabase_map$width+50
                                MetabaseOutput(ns("my_widget"), width = plot_width, height = plot_height)
                              })
                            }
                            }
                        })
                        
                        observeEvent(input$keggSave, {
                          #ID = input$x3
                          ID=input$sel_kegg_set
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
                          analysis_type = input$analysis_type_2
                          ID = input$x2
                          #validate(need(ID!="", message = "Select one geneset by clicking a GeneSet name from 'Gene Set Enrichment Analysis (GSEA)' or 'Over-Representation Analysis (ORA)' tab."))
                          
                          DataIn = DataReactive()
                          data_long = DataIn$data_long
                          ProteinGeneName = DataIn$ProteinGeneName
                          
                          if (input$gs_heatmap_DE=="All Genes") {
                            gsets_GSEA<-gsets_Reactive()
                            GenesetSig=gsets_GSEA[[ID]]
                          } else if (analysis_type == "GSEA"){
                            res<-gsea_results()
                            GenesetSig<-res%>%dplyr::filter(GeneSet==ID)%>%dplyr::select(leadingEdge)%>%unlist
                          } else {
                            res<-ora_results()
                            GenesetSig<-res%>%dplyr::filter(GeneSet==ID)%>%dplyr::select(DeGene_in_Set)
                            GenesetSig<-str_split(GenesetSig, ",")[[1]]
                          }
                          #browser()
                          if (analysis_type == "GSEA") {
                            getresults <- DataGenesetReactive_GSEA()
                            terminals.df <- getresults$GSEA.terminals.df
                          } else if (analysis_type == "ORA") {
                            getresults <- DataGenesetReactive_ORA()
                            terminals.df <- getresults$terminals.df
                          }
                          
                          
                          terminalsdf.set <- dplyr::filter(terminals.df, Gene.Name %in% GenesetSig)
                          terminals_id<-terminalsdf.set%>%dplyr::select(UniqueID)%>%unlist%>%unname
                          
                          #browser()
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
                          ID = input$x2
                          validate(need(ID!="", message = "Select one geneset by clicking a GeneSet name from 'Gene Set Enrichment Analysis (GSEA)' or 'Over-Representation Analysis (ORA)' tab."))
                          #grid.draw(genesetheatmap_out()$gtable)
                          draw(genesetheatmap_out(), merge_legend=T,  auto_adjust = FALSE)
                        }) 
                        
                        output$plot.geneset.heatmap=renderUI({
                          plotOutput(ns('SetHeatMap'), height=input$geneset_heatmap_height, width = input$geneset_heatmap_width)
                        })
                        
                        
                        observeEvent(input$genesetheatmap, {
                          #ID = input$x3
                          # saved_plots$genesetheatmap[[ID]] <- genesetheatmap_out() #this only works on R4.0
                          saved_plots$genesetheatmap<- genesetheatmap_out() #this works on R3.5 - 3.6
                          
                          #saved_plots$genesetheatmap <- genesetheatmap_out()$gtable
                        }
                        )
                        
                        geneset_Expression_data<-reactive({ 
                          analysis_type = input$analysis_type_1
                          ID = input$x1
                          gsets_GSEA<-gsets_Reactive()
                          GenesetSig=gsets_GSEA[[ID]]
                          if (analysis_type == "GSEA") {
                            getresults <- DataGenesetReactive_GSEA()
                            terminals.df <- getresults$GSEA.terminals.df
                            genes=sel_GSEA_set()$genes
                            NES=sel_GSEA_set()$table$NES
                            if (input$gene_rank_method=="logFC") {
                              if (NES>0) {  #sort so leading edge genes are shown first
                                terminals.df <-terminals.df %>%arrange(dplyr::desc(logFC))
                              } else { terminals.df <-terminals.df %>%arrange(logFC)}
                            } else {
                              if (NES>0) {  #sort so leading edge genes are shown first
                                terminals.df <-terminals.df %>%arrange(dplyr::desc(minus_logPvalue))
                              } else { terminals.df <-terminals.df %>%arrange(minus_logPvalue)}
                            }
                          } else if (analysis_type == "ORA") {
                            getresults <- DataGenesetReactive_ORA()
                            terminals.df <- getresults$terminals.df
                            genes=intersect(names(getresults$sig_genes_Dir), GenesetSig)
                          }
                          #cat(date(), analysis_type, ID, "\n")\
                          terminalsdf.set <- dplyr::filter(terminals.df, Gene.Name %in% GenesetSig)
                          terminalsdf.set[,sapply(terminalsdf.set,is.numeric)] <- signif(terminalsdf.set[,sapply(terminalsdf.set,is.numeric)],3)
                          return(list("terminalsdf.set"=terminalsdf.set, "genes"=genes))
                        })
                        
                        output$Expression <-  DT::renderDT(server=FALSE,{
                          ID = input$x1
                          validate(need(ID!="", message = "Select one geneset by clicking a GeneSet name from 'Gene Set Enrichment Analysis (GSEA)' or 'Over-Representation Analysis (ORA)' tab."))
                          analysis_type = input$analysis_type_1
                          terminalsdf.set<-geneset_Expression_data()$terminalsdf.set
                          genes<-geneset_Expression_data()$genes
                          #cat(analysis_type, ID, nrow(terminalsdf.set), "table", length(genes), "DEGs\n")
                          DT::datatable(terminalsdf.set,  extensions = 'Buttons', options = list( 
                            dom = 'lBfrtip', pageLength = 15,
                            buttons = list(
                              list(extend = "csv", text = "Download Page", filename = "Page_results",
                                   exportOptions = list(modifier = list(page = "current"))),
                              list(extend = "csv", text = "Download All", filename = "All_Results",
                                   exportOptions = list(modifier = list(page = "all")))
                            ))) %>% formatStyle(columns="Gene.Name", target="row",   backgroundColor = styleEqual(genes, "lightgreen", default='white'))
                        })
                      } )
}

