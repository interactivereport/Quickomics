###########################################################################################################
## Proteomics Visualization R Shiny App
##
##This software belongs to Biogen Inc. All right reserved.
##
##@file: ui.R
##@Developer : Benbo Gao (benbo.gao@Biogen.com)
##@Date : 9/6/2019
##@version 3.0
###########################################################################################################
source("global.R",local = TRUE)$value

# ui must be a function(request) (not a bare fluidPage() object) for Shiny's
# bookmarking (enableBookmarking = "server" on shinyApp() below) to be able
# to restore a saved session's URL-encoded state ID.
ui <- function(request) {
fluidPage(
  theme = shinytheme("cerulean"),
  windowTitle = "Quickomics",
  tagList(tags$head(tags$style(type = 'text/css','.navbar-brand{display:none;}')),
          useShinyjs(),
          rclipboard::rclipboardSetup(),
          tags$head(
            tags$link(rel = "stylesheet", type = "text/css",
                      href = "datatables/jquery.dataTables.min.css"),
            tags$script(src = "datatables/jquery.dataTables.min.js"),
            tags$script(src = "multidrag.js")
          ),
          tags$head(
            uiOutput("dynamic_sidebar_css"),
            tags$style(HTML("
    /* Top-right, fixed -- clear of the logo (top-left) so there's no overlap.
       Always visible now: no toggle button, no hidden/collapsed state. */
    #sidebar_width_control {
      position: fixed;
      top: 6px;
      right: 15px;
      z-index: 1000;
      background: rgba(255,255,255,0.95);
      padding: 8px 12px;
      border-radius: 4px;
      border: 1px solid #ccc;
      box-shadow: 0 1px 3px rgba(0,0,0,0.15);
    }
    #sidebar_width_control label { font-size: 12px; margin-bottom: 2px; }
    /* Venn Diagram Intersection Output tables: forces the Intersection
       column to wrap (at the <br> tags inserted server-side) inside its
       narrow columnDefs width, instead of overflowing. */
    .wrap-cell { white-space: normal !important; word-break: break-word; }
  "))
          ),
          div(id = "sidebar_width_control",
              sliderInput("sidebar_width_pct", "Set the Side Menu Width:", min = 15, max = 50, value = 25, step = 1, ticks = FALSE, width = "180px")
          ),

          titlePanel(
            fluidRow(
              column(4, img(height =75 , src = "Quickomics.png")),
              # padding-right reserves space for #sidebar_width_control (a
              # position:fixed box anchored to the viewport's top-right corner,
              # so it occupies roughly the rightmost 220px regardless of window
              # width) -- without this, a long project name can render right
              # underneath it as the window narrows and this column shrinks.
              column(8,  h2(strong(textOutput('project')), align = 'left'), style = "padding-right: 230px;")
            ),
            windowTitle = "Quickomics" ),
          
          
          navbarPage(title ="", id="menu",                     
                     ##########################################################################################################
                     ## Select Dataset
                     ##########################################################################################################
                     tabPanel("Dataset",
                              fluidRow(
                                column(12,
                                       tabsetPanel(id="Tables",
                                                   #tabPanel(title="Introduction",htmlOutput('intro')),
                                                   #tabPanel(title="Project Table", DT::dataTableOutput('projecttable')),
                                                   tabPanel(title="Select Dataset",
                                                            radioButtons("select_dataset",label="Select data set", choices=c("Saved Projects","Upload RData File", "Upload Data Files (csv)"),inline = F, selected="Saved Projects"),
                                                            conditionalPanel("input.select_dataset=='Saved Projects'",
                                                                             selectInput("sel_project", label="Available Dataset",
                                                                                         choices=c("", projects), selected=NULL)),
                                                            conditionalPanel("input.select_dataset=='Upload RData File'",
                                                                             fileInput("file1", "Choose data file"),
                                                                             fileInput("file2", "(Optional) Choose network file"),
                                                                             uiOutput('ui.action') ),
                                                            conditionalPanel("input.select_dataset=='Upload Data Files (csv)'",
                                                                             h5("Use the Upload Files tab to the right to create your own data set.") )
                                                   ),
                                                   tabPanel(title="Project Overview", htmlOutput("summary"),
                                                            tableOutput('group_table'), tags$br(),
                                                            textInput("exp_unit", "Expression Data Units", value="Expression Level", width="300px"),tags$br(),
                                                            uiOutput('comp_info') ),
                                                   tabPanel(title="Sample Table", value = "sample_table", actionButton("sample", "Save to output"), dataTableOutput('sample')),
                                                   # uiOutput("conditinal_sample_table"),
                                                   tabPanel(title="Result Table", actionButton("results", "Save to output"), dataTableOutput('results')),
                                                   tabPanel(title="Data Table", value = "data_table", actionButton("data_wide", "Save to output"), dataTableOutput('data_wide')),
                                                   tabPanel(title="Protein Gene Names", actionButton("ProteinGeneName", "Save to output"), dataTableOutput('ProteinGeneName')),
                                                   tabPanel(title="Upload Files",
                                                            uiOutput('upload.files.ui'), #see process_uploaded_files.R for details.
                                                            tags$br(),
                                                            textOutput('upload.message')
                                                   ),
                                                   tabPanel(title="Help", htmlOutput('help_input'))
                                       )
                                )
                              )
                     ),
                     
                     ##########################################################################################################
                     ## Groups and Samples
                     ##########################################################################################################
                     tabPanel("Groups and Samples",
                              #  fluidRow(
                              #    column(4,
                              #           wellPanel(
                              #             tags$style(mycss),
                              #             actionButton("reset_group", "Reset Groups and Samples"),
                              #             tags$p("Use the boxes below to remove or add groups/samples. Use the tools at right side for advanced selection and ordering."),
                              #             selectInput("QC_groups_var", label="Select a Attribute (Covariate)", choices=NULL),
                              #             selectizeInput("QC_groups", label="Select Groups", choices=NULL, multiple=TRUE),
                              #             checkboxInput("QC_comp2sample", "Use Samples in Subset/Comparison?",  FALSE, width="90%"),
                              #             conditionalPanel("input.QC_comp2sample==1",
                              #                              uiOutput("QC_samples_from_comp")	),
                              #             selectizeInput("QC_samples", label="Select Samples", choices=NULL,multiple=TRUE),
                              #             column(width=12,uiOutput("selectGroupSample")),
                              #             tags$hr()
                              #           )),
                              #    column(8,
                              #      tags$br(),
                              #      tags$hr(style="border-color: RoyalBlue;"),
                              #      uiOutput('reorder_group'),
                              #      uiOutput('sample_choose_order')),
                              # 
                              # )
                              fluidRow(
                                column(1,
                                ),
                                column(10,
                                       wellPanel(
                                         actionButton("reset_all_types", "Reset All"),
                                         tagList(tags$div(
                                           tags$p("Tip: Use Ctrl (or Cmd on Mac) + click to pick multiple values.")
                                         )),  
                                         # dynamic UI region
                                         div(
                                           uiOutput('ui_all_types'),
                                           tags$hr(style="border-color: black;"),
                                           uiOutput("ui_source_s"),
                                           uiOutput("ui_dest_s"),
                                           tags$hr(style="border-color: black;"),
                                           uiOutput("ui_source_test"),
                                           uiOutput("ui_dest_test"),
                                           uiOutput("reset_test"),
                                           tags$hr(style="border-color: black;")
                                         ),
                                         
                                         # static summary region (must NOT be inside dynamic UI)
                                         div(
                                           span(textOutput("selectGroupSample"),
                                                style = "color:red; font-size:20px; font-family:arial; font-style:italic"),
                                           tags$hr(style="border-color: grey;"),
                                           span(verbatimTextOutput("summaryDetail"),
                                                style = "color:red; font-size:20px; font-family:arial; font-style:italic")
                                         )
                                       )
                                ),
                                column(1,
                                )
                              )
                     ),
                     
                     ##########################################################################################################
                     ## QC Plots
                     ##########################################################################################################
                     tabPanel("QC Plots", value = "QC_Plots",
                              fluidRow(
                                column(3,
                                       wellPanel(
                                         tags$style(mycss),
                                         column(width=12,uiOutput("selectGroupSampleQC")),
                                         conditionalPanel("input.groupplot_tabset=='PCA Plot' || input.groupplot_tabset=='PCA 3D Interactive' || input.groupplot_tabset=='PCA 3D Plot' || input.sd_tabset=='Distance Heatmap' ",
                                                          selectInput("PCAcolorby", label="Color By", choices=NULL)),
                                         conditionalPanel("input.groupplot_tabset=='PCA Plot' || input.groupplot_tabset=='PCA 3D Interactive'",
                                                          selectInput("PCAshapeby", label="Shape By", choices=NULL)),
                                         conditionalPanel("input.groupplot_tabset=='PCA Plot'",
                                                          selectInput("PCAsizeby", label="Size By", choices=NULL),
                                                          selectizeInput("pcnum",label="Select Principal Components", choices=1:10, multiple=TRUE, selected=1:2, options = list(maxItems = 2)),
                                                          radioButtons("ellipsoid", label="Plot Ellipsoid (>3 per Group)", inline = TRUE, choices = c("No" = FALSE,"Yes" = TRUE)),
                                                          radioButtons("mean_point", label="Show Mean Point", inline = TRUE, choices = c("No" = FALSE,"Yes" = TRUE)),
                                                          radioButtons("rug", label="Show Marginal Rugs", inline = TRUE, choices = c("No" = FALSE,"Yes" = TRUE)),
                                                          selectInput("PCAcolpalette", label= "Select palette", choices=c("Accent(8)"="Accent","Dark2(8)"="Dark2","Paired(12)"="Paired","Pastel1(9)"="Pastel1",
                                                                                                                          "Pastel2(8)"="Pastel2","Set1(9)"="Set1","Set2(8)"="Set2","Set3(12)"="Set3", 
                                                                                                                          "NPG(10)"="npg", "AAAS(10)"="aaas", "NEJM(8)"="nejm", "Lancet(9)"="lancet",
                                                                                                                          "JAM(7)"="jama", "D3(10)"="d3", "Uchiago(9)"="uchicago"), selected="Dark2"),
                                                          sliderInput("PCAdotsize", "Dot Size (when size by not used):", min = 1, max = 20, step = 1, value = 4),
                                                          radioButtons("PCA_subsample", label="Label Samples:", inline = TRUE, choices = c("All","None", "Subset"), selected = "All"),
                                                          conditionalPanel("input.PCA_subsample!='None'",
                                                                           sliderInput("PCAfontsize", "Label Font Size:", min = 1, max = 20, step = 1, value = 10),
                                                                           radioButtons("PCA_label",label="Select Sample Label",inline = TRUE, choices="")),
                                                          conditionalPanel("input.PCA_subsample=='Subset'",
                                                                           actionButton("PCA_refresh_sample", "Reload Sample IDs"),
                                                                           textAreaInput("PCA_list", "List of Samples to Label", "", cols = 5, rows=6))
                                         ),
                                         conditionalPanel("input.groupplot_tabset=='Covariates'",
                                                          selectizeInput("covar_variates", label="Select Covariates:", choices=NULL, multiple = TRUE),
                                                          numericInput("covar_PC_cutoff", label= "Principle Component (PC) Cutoff (% explained variance)",  value = 5, min=0, max=100, step=1),
                                                          numericInput("covar_FDR_cutoff", label= "Choose FDR Value Cutoff", value=0.1, min=0, max=1, step=0.001),
                                                          sliderInput("covar_ncol", label= "Number of Columns for Plots", min = 1, max = 6, step = 1, value = 3)
                                         ),
                                         conditionalPanel("input.groupplot_tabset=='PCA 3D Plot'",
                                                          radioButtons("ellipsoid3d", label="Plot Ellipsoid (>3 per Group)", inline = TRUE, choices = c("No" = "No","Yes" = "Yes")),
                                                          radioButtons("dotlabel", label="Dot Label", inline = TRUE, choices =  c("No" = "No","Yes" = "Yes"))
                                         ),
                                         conditionalPanel("input.groupplot_tabset=='Dendrograms'",
                                                          sliderInput("DendroCut", label="tree cut number:", min = 2, max = 10, step = 1, value = 4),
                                                          sliderInput("DendroFont", label= "Label Font Size:", min = 0.5, max = 4, step = 0.5, value = 1),
                                                          radioButtons("dendroformat", label="Select Plot Format", inline = TRUE, choices = c("tree" = "tree","horizontal" = "horiz", "circular" = "circular"), selected="circular")
                                         ),
                                         conditionalPanel("input.groupplot_tabset=='Sample-sample Distance' && input.sd_tabset=='Distance Heatmap'",
                                                          radioButtons("sd_dist_cluster", label="Apply Clustering", inline=TRUE,
                                                                       choices=c("Yes"="TRUE","No"="FALSE"), selected="TRUE"),
                                                          radioButtons("sd_dist_more", label="Show More Options", inline=TRUE,
                                                                       choices=c("Yes","No"), selected="No"),
                                                          conditionalPanel("input.sd_dist_more=='Yes'",
                                                                           tags$hr(),
                                                                           # h6("Font Sizes"),
                                                                           fluidRow(
                                                                             column(width=6, sliderInput("sd_dist_row_font",  "Row Labels:",    min=4, max=20, step=1, value=10)),
                                                                             column(width=6, sliderInput("sd_dist_col_font",  "Col Labels:",    min=4, max=20, step=1, value=10))
                                                                           ),
                                                                           fluidRow(
                                                                             column(width=6, sliderInput("sd_dist_lgd_title", "Legend Title:",  min=4, max=20, step=1, value=12)),
                                                                             column(width=6, sliderInput("sd_dist_lgd_labels","Legend Labels:", min=4, max=20, step=1, value=10))
                                                                           ),
                                                                           # fluidRow(
                                                                           #   column(width=6, sliderInput("sd_dist_ann_title", "Annot. Titles:", min=4, max=20, step=1, value=10)),
                                                                           #   column(width=6, sliderInput("sd_dist_ann_labels","Annot. Labels:", min=4, max=20, step=1, value=10))
                                                                           # ),
                                                                           tags$hr(),
                                                                           sliderInput("sd_dist_height", "Plot Height:", min=200, max=1500, step=50, value=800)
                                                          )
                                         ),
                                         conditionalPanel("input.groupplot_tabset=='Sample-sample Distance' && input.sd_tabset=='Correlation Heatmap'",
                                                          radioButtons("sd_cor_cluster", label="Apply clustering", inline=TRUE,
                                                                       choices=c("Yes"="TRUE","No"="FALSE"), selected="FALSE"),
                                                          selectizeInput("sd_cor_annotate_by", label="Annotate By",
                                                                         choices=NULL, multiple=TRUE,
                                                                         options=list(placeholder="None (no annotation)")),
                                                          selectInput("sd_cor_label_by", label="Axis Label By", choices=NULL),
                                                          radioButtons("sd_cor_more", label="Show More Options", inline=TRUE,
                                                                       choices=c("Yes","No"), selected="No"),
                                                          conditionalPanel("input.sd_cor_more=='Yes'",
                                                                           tags$hr(),
                                                                           h6("Heatmap Color"),
                                                                           fluidRow(
                                                                             column(width=4, colourInput("sd_cor_lowColor",  "Low",  "#ADD8E6")),
                                                                             column(width=4, colourInput("sd_cor_midColor",  "Mid",  "#FFFF00")),
                                                                             column(width=4, colourInput("sd_cor_highColor", "High", "#FF0000"))
                                                                           ),
                                                                           tags$hr(),
                                                                           h6("Annotation Color Setting"),
                                                                           radioButtons("sd_cor_annot_color", label=NULL, inline=FALSE,
                                                                                        choices=c("Auto-Set by Rand. Seed","Select Palette","Upload Colors"),
                                                                                        selected="Auto-Set by Rand. Seed"),
                                                                           conditionalPanel("input.sd_cor_annot_color=='Auto-Set by Rand. Seed'",
                                                                                            numericInput("sd_cor_hm_seed", label="Random Seed", min=1, max=5000, value=123, step=1)
                                                                           ),
                                                                           conditionalPanel("input.sd_cor_annot_color=='Select Palette'",
                                                                                            selectizeInput("sd_cor_cat_pal", label="Category Palettes (one per attribute)",
                                                                                                           choices=c("Dark2(8)"="Dark2","Accent(8)"="Accent","Set1(9)"="Set1","Set2(8)"="Set2","Set3(12)"="Set3","NPG(10)"="npg","NEJM(8)"="nejm","Lancet(9)"="lancet","JAMA(7)"="jama","D3(10)"="d3","Uchicago(9)"="uchicago"), multiple=TRUE),
                                                                                            selectInput("sd_cor_num_pal", label="Numeric Annotation Palette",
                                                                                                        choices=c("Dark2(8)"="Dark2","Accent(8)"="Accent","Set1(9)"="Set1","Set2(8)"="Set2","Set3(12)"="Set3","NPG(10)"="npg","NEJM(8)"="nejm","Lancet(9)"="lancet","JAMA(7)"="jama","D3(10)"="d3","Uchicago(9)"="uchicago"), selected="Set1")
                                                                           ),
                                                                           conditionalPanel("input.sd_cor_annot_color=='Upload Colors'",
                                                                                            fileInput("sd_cor_annot_color_file", "Upload Color CSV", accept=".csv")
                                                                           ),
                                                                           tags$hr(),
                                                                           h6("Font Sizes"),
                                                                           fluidRow(
                                                                             column(width=6, sliderInput("sd_cor_row_font",  "Row Labels:",    min=4, max=20, step=1, value=10)),
                                                                             column(width=6, sliderInput("sd_cor_col_font",  "Col Labels:",    min=4, max=20, step=1, value=10))
                                                                           ),
                                                                           fluidRow(
                                                                             column(width=6, sliderInput("sd_cor_lgd_title", "Legend Title:",  min=4, max=20, step=1, value=12)),
                                                                             column(width=6, sliderInput("sd_cor_lgd_labels","Legend Labels:", min=4, max=20, step=1, value=10))
                                                                           ),
                                                                           fluidRow(
                                                                             column(width=6, sliderInput("sd_cor_ann_title", "Annot. Titles:", min=4, max=20, step=1, value=10)),
                                                                             column(width=6, sliderInput("sd_cor_ann_labels","Annot. Labels:", min=4, max=20, step=1, value=10))
                                                                           ),
                                                                           tags$hr(),
                                                                           sliderInput("sd_cor_height", "Plot Height:", min=200, max=1500, step=50, value=800)
                                                          )
                                         ),
                                         conditionalPanel("input.groupplot_tabset=='AlignQC'",
                                                          tags$hr(),
                                                          conditionalPanel("input.AlignQC_tabset=='Plot Other Variables'",
                                                                           selectizeInput("alignQC_var", label="Select Numerica Variable to Plot", choices=NULL,multiple=FALSE)),
                                                          conditionalPanel("input.AlignQC_tabset!='Top Gene List'",
                                                                           sliderInput("alignQC_height", "Plot Height:", min = 200, max = 3000, step = 50, value = 700),
                                                                           sliderInput("alignQC_width", "Plot Width:", min = 200, max = 3000, step = 50, value = 1000)),
                                                          conditionalPanel("input.AlignQC_tabset=='Top Gene List'",
                                                                           sliderInput("alignQC_TG_Ngene", "Max # of Top Genes Per Sample:", min = 1, max = 100, step = 1, value = 30),
                                                                           sliderInput("alignQC_TG_Ntotal", "Max # of Genes to Show", min = 10, max = 200, step = 5, value = 90),
                                                                           checkboxInput("alignQC_TG_list", "Show Table Instead of Graph?",  FALSE, width="90%"),
                                                                           conditionalPanel("input.alignQC_TG_list==0",
                                                                                            sliderInput("alignQC_TG_height", "Gene Plot Height:", min = 200, max = 3000, step = 50, value = 800),
                                                                                            sliderInput("alignQC_TG_width", "Gene Plot Width:", min = 200, max = 3000, step = 50, value = 800))),
                                                          conditionalPanel("input.AlignQC_tabset!='Top Gene List' || input.alignQC_TG_list==0",
                                                                           sliderInput("alignQC_fontsize", "Label Font Size:", min = 1, max = 18, step = 1, value = 12)),
                                                          conditionalPanel("input.AlignQC_tabset=='Top Gene Ratio' || input.AlignQC_tabset=='Top Gene List'",
                                                                           radioButtons("convert_exp", label="Convert Expression Data", inline = TRUE, choices = c("None", "log to linear"), selected = "log to linear"),
                                                                           conditionalPanel("input.convert_exp=='log to linear'",
                                                                                            column(width=6,numericInput("convert_logbase", label= "log base",  value = 2, min=1, step=1)),
                                                                                            column(width=6,numericInput("convert_small", label= "small value",  value=1, min=0, step=0.1))),
                                                                           conditionalPanel("input.AlignQC_tabset!='Top Gene List' | input.alignQC_TG_list==0",
                                                                                            textInput("alignQC_TGR_Y", "Label for Top Gene %", value="% of Total TPM", width = "100%")))
                                         )
                                       )),
                                column(9,
                                       tabsetPanel(id="groupplot_tabset",
                                                   tabPanel(title="PCA Plot", actionButton("pcaplot", "Save to output"),
                                                            actionButton("plot_PCA", "Plot/Refresh", style="color: #0961E3; background-color: #F6E98C ; border-color: #2e6da4"),
                                                            plotOutput("pcaplot",height = 800)),
                                                   tabPanel(title="Covariates",
                                                            tabsetPanel(id="covartiate_tabset",
                                                                        tabPanel(title="Summary", actionButton("compute_PC", "Compute/Refresh", style="color: #0961E3; background-color: #F6E98C ; border-color: #2e6da4"),
                                                                                 textOutput("N_pairs"),tags$br(),dataTableOutput('covar_table')),
                                                                        tabPanel(title="Categorical Covariates",
                                                                                 h5("After changing parameters, please click Refresh button in the Summary panel to generate new plots."),
                                                                                 actionButton("covar_cat", "Save to output"),
                                                                                 sliderInput("covar_cat_height", "Plot Height:", min = 200, max = 3000, step = 50, value = 800),
                                                                                 tags$br(), uiOutput("plot.PC_covariatesC") ),
                                                                        tabPanel(title="Numerical Covariates",
                                                                                 h5("After changing parameters, please click Refresh button in the Summary panel to generate new plots."),
                                                                                 actionButton("covar_num", "Save to output"),
                                                                                 sliderInput("covar_num_height", "Plot Height:", min = 200, max = 3000, step = 50, value = 800),
                                                                                 tags$br(),uiOutput("plot.PC_covariatesN") )
                                                            )),
                                                   tabPanel(title="AlignQC",
                                                            tabsetPanel(id="AlignQC_tabset",
                                                                        tabPanel(title="Top Gene Ratio",
                                                                                 actionButton("alignQC_TGR", "Save to output"),
                                                                                 uiOutput("alignQC_TGR_plot")),
                                                                        tabPanel(title="Top Gene List",
                                                                                 actionButton("alignQC_TGL", "Save to output"),
                                                                                 uiOutput("alignQC_TG_plot")),
                                                                        tabPanel(title="Mapped Read Allocation",
                                                                                 actionButton("alignQC_RA", "Save to output"),
                                                                                 uiOutput("alignQC_RA_plot")),
                                                                        tabPanel(title="Plot Other Variables",
                                                                                 actionButton("alignQC_OV", "Save to output"),
                                                                                 uiOutput("alignQC_OV_plot"))
                                                            )),
                                                   tabPanel(title="Eigenvalues",  plotOutput("Eigenvalues",height = 650)),
                                                   tabPanel(title="PCA 3D Plot",  plotOutput("pca_legend",height = 100), rglwidgetOutput("plot3d",  width = 1000, height = 1000)),
                                                   tabPanel(title="PCA 3D Interactive", plotlyOutput("plotly3d",  width = 1000, height = 1000)),
                                                   tabPanel(title="Sample-sample Distance",
                                                            tabsetPanel(id="sd_tabset",
                                                                        tabPanel(title="Distance Heatmap",
                                                                                 actionButton("SampleDistance", "Save to output"),
                                                                                 uiOutput("sd_heatmap_plot_ui")
                                                                        ),
                                                                        tabPanel(title="Correlation Heatmap",
                                                                                 actionButton("SampleCorrelation", "Save to output"),
                                                                                 uiOutput("sd_cor_plot_ui")
                                                                        ),
                                                                        tabPanel(title="Correlation Data Table",
                                                                                 DT::dataTableOutput("pheatmap_table")
                                                                        )
                                                            )
                                                   ),
                                                   tabPanel(title="Dendrograms",actionButton("Dendrograms", "Save to output"), plotOutput("Dendrograms",height = 800)),
                                                   tabPanel(title="Box Plot", actionButton("QCboxplot", "Save to output"), plotOutput("QCboxplot",height = 800)),
                                                   tabPanel(title="CV Distribution", actionButton("histplot", "Save to output"), plotOutput("histplot",height = 800)),
                                                   tabPanel(title="Help", htmlOutput('help_QC'))
                                       )
                                )
                              )
                     ),
                     
                     ##########################################################################################################
                     ## Gene/Protein Expression Plot
                     ##########################################################################################################
                     tabPanel("Expression Plot", value = 'Exp_Plot',
                              fluidRow(
                                column(3,
                                       wellPanel(
                                         conditionalPanel("input.expression_tabset=='Searched Expression Data' | input.expression_tabset=='Browsing' ",
                                                          # selectizeInput("sel_group", label="Select Groups", choices=NULL, multiple=TRUE),
                                                          column(width=12,uiOutput("selectGroupSampleExpression"))
                                         ),
                                         conditionalPanel("input.expression_tabset=='Searched Expression Data' || input.expression_tabset=='Rank Abundance Curve'",
                                                          radioButtons("exp_subset",label="Genes Used in Plot", choices=c("Select", "Upload Genes", "Geneset"),inline = TRUE, selected="Select"),
                                                          conditionalPanel("input.exp_subset=='Upload Genes'",
                                                                           textAreaInput("exp_list", "Enter Gene List", "", cols = 5, rows=6)),
                                                          conditionalPanel("input.exp_subset=='Select'",
                                                                           radioButtons("exp_label",label="Select Gene Label",inline = TRUE, choices=c("UniqueID", "Gene.Name"), selected="Gene.Name"),
                                                                           selectizeInput("sel_gene",	label="Gene Name (Select 1 or more)",	choices = NULL,	multiple=TRUE, options = list(placeholder =	'Type to search'))),
                                                          conditionalPanel("input.exp_subset=='Geneset'",   uiOutput("html_geneset_exp") ) ),
                                         conditionalPanel("input.expression_tabset=='Searched Expression Data'",
                                                          radioButtons("SeparateOnePlot", label="Separate or One Plot", inline = TRUE, choices = c("Separate" = "Separate", "OnePlot" = "OnePlot"))),
                                         conditionalPanel("input.expression_tabset=='Browsing'",
                                                          column(width=6,numericInput("expression_fccut", label= "Choose Fold Change Threshold",  value = 1.2, min=1, step=0.1)),
                                                          column(width=6,numericInput("expression_pvalcut", label= "Choose P-value Threshold",  value=0.01, min=0, step=0.001)),
                                                          radioButtons("expression_psel", label= "P value or P.adj Value?",
                                                                       choices= c("Pval"="Pval","Padj"="Padj"),inline = TRUE),
                                                          selectInput("expression_test", label="Select Test", choices=NULL),
                                                          textOutput("expfilteredgene"),
                                                          tags$head(tags$style("#expfilteredgene{color: red; font-size: 20px; font-style: italic;}")),
                                                          column(width=6,selectInput("sel_page",	label="Select Page",	choices = 1,	selected=1)),
                                                          column(width=6,selectInput("numperpage", label= "Plot Number per Page", choices= c("4"=4,"6"=6,"9"=9), selected=6)),
                                                          radioButtons("browsing_gene_order", label="Order Genes by", inline = TRUE, choices = c("abs(Fold Change)","P value"))),
                                         conditionalPanel("input.expression_tabset=='Searched Expression Data' | input.expression_tabset=='Browsing' ",
                                                          #  selectizeInput("sel_group", label="Select Groups", choices=NULL, multiple=TRUE),
                                                          # column(width=12,textOutput("selectGroupSampleExpression")),
                                                          radioButtons("sel_geneid",label="Select Gene Label",inline = TRUE, choices=""),
                                                          radioButtons("plotformat", label="Select Plot Format", inline = TRUE, choices = c("Box Plot" = "boxplot","Bar Plot" = "barplot", "violin" = "violin","line" = "line")),
                                                          conditionalPanel("input.SeparateOnePlot=='OnePlot' & input.expression_tabset=='Searched Expression Data'",
                                                                           h5("OnePlot only supports Bar and Line plots")),
                                                          radioButtons("IndividualPoint", label="Show Individual Point?", inline = TRUE, choices = c("YES" = "YES","NO" = "NO")),
                                                          selectizeInput("plotx", label="X Axis", choices=NULL, multiple = FALSE),
                                                          conditionalPanel("input.SeparateOnePlot!='OnePlot' | input.expression_tabset=='Browsing'",
                                                                           selectizeInput("colorby", label="Color By:", choices=NULL, multiple = FALSE)),
                                                          conditionalPanel("input.colorby=='None'",
                                                                           colourInput("barcol", "Select colour", "#1E90FF", palette = "limited")),
                                                          conditionalPanel("input.colorby!='None'",
                                                                           selectInput("colpalette", label= "Select palette", choices=c("Accent(8)"="Accent","Dark2(8)"="Dark2","Paired(12)"="Paired","Pastel1(9)"="Pastel1",
                                                                                                                                        "Pastel2(8)"="Pastel2","Set1(9)"="Set1","Set2(8)"="Set2","Set3(12)"="Set3", "NPG(10)"="npg", "AAAS(10)"="aaas", "NEJM(8)"="nejm", 
                                                                                                                                        "Lancet(9)"="lancet", "JAM(7)"="jama", "D3(10)"="d3", "Uchiago(9)"="uchicago"), selected="Dark2")),
                                                          #radioButtons("ColPattern", label="Bar Colors", inline = TRUE, choices = c("Palette" = "Palette", "Single" = "Single")),
                                                          sliderInput("expression_axisfontsize", "Axis Font Size:", min = 10, max = 28, step = 1, value = 16),
                                                          sliderInput("expression_titlefontsize", "Title Font Size:", min = 12, max = 28, step = 1, value = 16),
                                                          sliderInput("exp_plot_ncol", label= "Column Number", min = 1, max = 6, step = 1, value = 3),
                                                          textInput("Ylab", "Y label", width = "100%"),
                                                          textInput("Xlab", "X label", width = "100%"),
                                                          sliderInput("Xangle", label= "X Angle", min = 0, max = 90, step = 5, value = 90),
                                                          radioButtons("exp_plot_Y_scale", label="Y Axis Scale", inline = TRUE, choices = c("Log","Linear"), selected = "Log"),
                                                          conditionalPanel("input.exp_plot_Y_scale=='Linear'",
                                                                           tags$p("Linear values are computed using the log base and samll value from the log based expression unit, e.g. log2(TPM+0.25). Please update the values below and the Y Label above as needed."),
                                                                           column(width=6,numericInput("linear_base", label= "log Base",  value = 2)),
                                                                           column(width=6,numericInput("linear_small_value", label= "Small Value",  value=0.25))),
                                                          radioButtons("exp_plot_Y_range", label="Y Axis Range", inline = TRUE, choices = c("Auto","Manual"), selected = "Auto"),
                                                          conditionalPanel("input.exp_plot_Y_range=='Manual'",
                                                                           column(width=6,numericInput("exp_plot_Ymin", label= "Y Min",  value = 0, step=0.1)),
                                                                           column(width=6,numericInput("exp_plot_Ymax", label= "Y Max",  value=5, step=0.1))),
                                                          h5("After changing parameters, please click Plot/Refresh button in the plot panel to generate expression plot.")),
                                         conditionalPanel("input.expression_tabset=='Rank Abundance Curve'",
                                                          sliderInput("scurve_axisfontsize", "Axis Font Size:", min = 12, max = 28, step = 4, value = 16),
                                                          sliderInput("scurve_labelfontsize", "Label Font Size:", min = 2, max = 12, step = 1, value = 6),
                                                          textInput("scurveYlab", "Y label", value="Abundance", width = "100%"),
                                                          textInput("scurveXlab", "X label", value="Rank", width = "100%"),
                                                          sliderInput("scurveXangle", label= "X Angle", min = 0, max = 90, step = 15, value = 45),
                                                          radioButtons("scurveright", label="Density or histogram on Right:", inline = TRUE, choices = c("densigram" = "densigram", "density" = "density","histogram" = "histogram","boxplot" = "boxplot", "violin"= "violin"))),
                                         conditionalPanel("input.expression_tabset=='expression_plot_data' || input.expression_tabset=='Result Table' ",
                                                          h5("Enter some genes in Search Expression Data tab, then come here for data table."),
                                                          radioButtons("exp_table_format", label="Table Format", inline = TRUE,
                                                                       choices = c("Wide Format" = "wide", "Long Format" = "long"),
                                                                       selected = "wide")
                                         )                                       
                                       )
                                ),
                                column(9,
                                       tabsetPanel(id="expression_tabset",
                                                   tabPanel(title="Browsing",actionButton("browsing", "Save to output"),
                                                            actionButton("plot_browsing", "Plot/Refresh", style="color: #0961E3; background-color: #F6E98C ; border-color: #2e6da4"),
                                                            plotOutput("browsing", height=800)),
                                                   tabPanel(title="Searched Expression Data",actionButton("boxplot", "Save to output"),
                                                            verbatimTextOutput("geneSearchInfo"),
                                                            actionButton("plot_exp", "Plot/Refresh", style="color: #0961E3; background-color: #F6E98C ; border-color: #2e6da4"),
                                                            uiOutput("plot.exp")),
                                                   tabPanel(title="Data Table",	value="expression_plot_data", DT::dataTableOutput("dat_dotplot")),
                                                   tabPanel(title="Result Table",	DT::dataTableOutput("res_dotplot")),
                                                   tabPanel(title="Rank Abundance Curve",plotOutput("SCurve", height=800)),
                                                   tabPanel(title="Help", htmlOutput('help_expression'))
                                       )
                                )
                              )
                     ),
                     
                     ##########################################################################################################
                     ## Heatmap
                     ##########################################################################################################
                     tabPanel("Heatmap", value = 'Heatmap',
                              fluidRow(
                                column(3,
                                       wellPanel(
                                         # in the Heatmap tab sidebar, uncomment/restore:
                                         #actionButton("action_heatmaps","Generate Interactive Heatmap"),
                                         # actionButton("action_heatmaps", "Plot/Refresh", style="color: #0961E3; background-color: #F6E98C ; border-color: #2e6da4"),
                                         #selectizeInput("heatmap_groups", label="Select Groups", choices=NULL, multiple=TRUE),
                                         #checkboxInput("HM_comp2sample", "Use Samples in Subset/Comparison?",  FALSE, width="90%"),
                                         #conditionalPanel("input.HM_comp2sample==1",
                                         #                 uiOutput("HM_samples_from_comp")	),
                                         #selectizeInput("heatmap_samples", label="Select Samples", choices=NULL, multiple=TRUE),
                                         column(width=12,uiOutput("selectGroupSampleHeatmap")),
                                         radioButtons("heatmap_subset",label="Genes used for heatmap", choices=c("All","Subset","Upload Genes", "Geneset"),inline = TRUE, selected="All"),
                                         conditionalPanel("input.heatmap_subset=='Upload Genes'",
                                                          radioButtons("heatmap_upload_type", label="Select upload type", inline = TRUE, choices = c("Gene List","Annotated Gene File"), selected = "Gene List"),
                                                          conditionalPanel("input.heatmap_upload_type=='Gene List'",
                                                                           textAreaInput("heatmap_list", "Enter Gene List", "", cols = 5, rows=6)),
                                                          conditionalPanel("input.heatmap_upload_type=='Annotated Gene File'",
                                                                           uiOutput("gene_annot_file"))),
                                         conditionalPanel("input.heatmap_subset=='Geneset'",   uiOutput("html_geneset_hm") ),
                                         conditionalPanel("input.heatmap_subset=='All'",
                                                          radioButtons("heatmap_submethod", label= "Plot Random Genes or Variable Genes", choices= c("Random"="Random","Variable"="Variable"),inline = TRUE),
                                                          numericInput("maxgenes",label="Set Number of Genes to Plot", min=1, max= 5000, value=100, step=1)),
                                         conditionalPanel("input.heatmap_subset=='Subset'",
                                                          selectInput("heatmap_test", label="Select Genes from Test:", choices=NULL),
                                                          column(width=6,numericInput("heatmap_fccut", label= "Fold Change Cutoff", value = 1.2, min=1, step=0.1)),
                                                          column(width=6,numericInput("heatmap_pvalcut", label= "P Value Cutoff", value=0.01, min=0, step=0.001)),
                                                          radioButtons("heatmap_psel", label= "P value or P.adj Value?", choices= c("Pval"="Pval","Padj"="Padj"),inline = TRUE),
                                                          column(width=12,textOutput("heatmapfilteredgene"),
                                                                 tags$head(tags$style("#heatmapfilteredgene{color: red; font-size: 16px; font-style: italic; }"))),
                                                          uiOutput("Test_to_sample"),
                                                          tags$hr()
                                         ),
                                         conditionalPanel( "input.heatmap_tabset=='Static Heatmap Layout 1'",
                                                           selectizeInput("heatmap_annot", label="Annotate Samples", choices=NULL, multiple = TRUE),
                                                           conditionalPanel("input.heatmap_annot && input.heatmap_annot.length > 0",
                                                                            radioButtons("heatmap_annot_reorder",
                                                                                         label = "Reorder Samples by Annotation Attributes?",
                                                                                         choices = c("No", "Yes"), selected = "Yes", inline = TRUE),
                                                                            conditionalPanel("input.heatmap_annot_reorder=='Yes'",
                                                                                             h5("Samples will be sorted in a nested/stratified order using the attributes above (first attribute = outermost grouping). Column clustering will be disabled.",
                                                                                                style = "color:red; font-size:13px; font-family:arial; font-style:italic")
                                                                            )
                                                           )),
                                         # fluidRow (not bare column()s) so this row's floats are cleared via
                                         # Bootstrap's .row clearfix -- otherwise, whenever "Apply Clustering:"
                                         # happens to wrap onto two lines while "Apply Scaling:" stays on one,
                                         # the uncleared height mismatch corrupts the float layout of every
                                         # sibling that follows (including other fluidRows below, since a row's
                                         # own clearfix only clears ITS children, not an uncleared row before it).
                                         fluidRow(
                                           column(width=5,selectInput("dendrogram", "Apply Clustering:", c("both" ,"none", "row", "column"))),
                                           column(width=5,selectInput("scale", "Apply Scaling:", c("none","row", "column"),selected="row"))
                                         ),
                                         conditionalPanel( "input.heatmap_tabset=='Static Heatmap Layout 2'",
                                                           fluidRow(
                                                             column(width=5,selectInput("key", "Color Key:", c("TRUE", "FALSE"))),
                                                             column(width=5,selectInput("srtCol", "angle of label", c("45", "60","90"))),
                                                             column(width=5,sliderInput("hxfontsize", "Column Font Size:", min = 0, max = 3, step = 0.5, value = 1)),
                                                             column(width=5,sliderInput("hyfontsize", "Row Font Size:", min = 0, max = 3, step = 0.5, value = 1)),
                                                             column(width=5,sliderInput("right", "Set Margin Width", min = 0, max = 20, value = 5)),
                                                             column(width=5,sliderInput("bottom", "Set Margin Height", min = 0, max = 20, value = 5))
                                                           )
                                         ),
                                         # conditionalPanel( "input.heatmap_tabset=='Interactive Heatmap'",
                                         #                   column(width=5,sliderInput("hxfontsizei", "Column Font Size:", min = 0, max = 3, step = 0.5, value = 1)),
                                         #                   column(width=5,sliderInput("hyfontsizei", "Row Font Size:", min = 0, max = 3, step = 0.5, value = 1))
                                         # ),
                                         conditionalPanel( "input.heatmap_tabset=='Static Heatmap Layout 1'",
                                                           fluidRow(
                                                             column(width=5,sliderInput("hxfontsizep", "Column Font Size:", min = 0, max = 20, step = 1, value = 10)),
                                                             column(width=5,sliderInput("hyfontsizep", "Row Font Size:", min = 0, max = 20, step = 1, value = 7))
                                                           ),
                                                           radioButtons("heatmap_label",label="Gene Label",inline = TRUE, choices=""),
                                                           sliderInput("heatmap_N_genes", "Max Number of Genes to Label:", min = 0, max = 500, step = 10, value = 100),
                                                           h5("After changing parameters, please click Plot/Refresh button in the plot panel to generate heatmap."),
                                                           radioButtons("heatmap_more_options", label="Show More Options", inline = TRUE, choices = c("Yes","No"), selected = "No"), 
                                                           conditionalPanel("input.heatmap_more_options=='Yes'",
                                                                            radioButtons("heatmap_annot_color",label="Color Setting for Annotations", choices=c("Auto-Set by Rand. Seed","Select Palette","Upload Colors"), selected="Auto-Set by Rand. Seed"),
                                                                            conditionalPanel("input.heatmap_annot_color=='Auto-Set by Rand. Seed'",
                                                                                             numericInput("hm_seed",label="Random Seed for Color Palettes for Categories", min=1, max= 5000, value=123, step=1)),
                                                                            conditionalPanel("input.heatmap_annot_color=='Select Palette'",
                                                                                             selectizeInput("heatmap_cat_pal", label="Color Palettes (One per Category Annotation)", 
                                                                                                            choices=c("Dark2(8)"="Dark2", "Accent(8)"="Accent",  "Set1(9)"="Set1","Set2(8)"="Set2","Set3(12)"="Set3", "NPG(10)"="npg",
                                                                                                                      "NEJM(8)"="nejm", "Lancet(9)"="lancet", "JAM(7)"="jama", "D3(10)"="d3", "Uchiago(9)"="uchicago"), multiple = TRUE),
                                                                                             selectizeInput("heatmap_num_pal", label="Set Max Colors for Numerical Annotations",
                                                                                                            choices=c("Dark2(8)"="Dark2", "Accent(8)"="Accent",  "Set1(9)"="Set1","Set2(8)"="Set2","Set3(12)"="Set3", "NPG(10)"="npg",
                                                                                                                      "NEJM(8)"="nejm", "Lancet(9)"="lancet", "JAM(7)"="jama", "D3(10)"="d3", "Uchiago(9)"="uchicago"), selected="Set1", multiple = FALSE)),
                                                                            conditionalPanel("input.heatmap_annot_color=='Upload Colors'",    
                                                                                             uiOutput("annot_color_file")),
                                                                            tags$hr(),
                                                                            radioButtons("heatmap_highlight", label="Highlight Subset of Genes:", inline = TRUE, choices = c("Yes","No"), selected = "No"),
                                                                            conditionalPanel("input.heatmap_highlight=='Yes'", 
                                                                                             uiOutput("gene_highlight_file"),
                                                                                             sliderInput("hl_font_size", "Font Size:", min = 0, max = 20, step = 1, value = 9)
                                                                            ),
                                                                            sliderInput("heatmap_height", "Heatmap Height:", min = 200, max = 3000, step = 50, value = 800),
                                                                            radioButtons("heatmap_row_dend", label="Show row dendrogram", inline = TRUE, choices = c("No" = FALSE,"Yes" = TRUE), selected = TRUE),
                                                                            radioButtons("heatmap_col_dend", label="Show column dendrogram", inline = TRUE, choices = c("No" = FALSE,"Yes" = TRUE), selected = TRUE),
                                                                            fluidRow(
                                                                              column(width=3,colourInput("lowColor", "Low", "blue")),
                                                                              column(width=3,colourInput("midColor", "Mid", "white")),
                                                                              column(width=3,colourInput("highColor", "High", "red"))
                                                                            ),
                                                                            fluidRow(
                                                                              column(width=5,selectInput("distanceMethod", "Distance Metric:", c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"))),
                                                                              column(width=5,selectInput("agglomerationMethod", "Linkage Algorithm:", c("complete", "single", "average", "centroid", "median", "mcquitty", "ward.D", "ward.D2"))),
                                                                              column(width=5,sliderInput("cutreerows", "cutree_rows:", min = 0, max = 8, step = 1, value = 0)),
                                                                              column(width=5,sliderInput("cutreecols", "cutree_cols:", min = 0, max = 8, step = 1, value = 0))
                                                                            ),
                                                                            fluidRow(
                                                                              textInput("heatmap_row_title", "Row Title", width = "100%"),
                                                                              sliderInput("heatmap_row_title_font_size", "Row Title Font Size:", min = 0, max = 30, step = 1, value = 16),
                                                                              textInput("heatmap_column_title", "Column Title", width = "100%"),
                                                                              sliderInput("heatmap_column_title_font_size", "Column Title Font Size:", min = 0, max = 30, step = 1, value = 16)                                           
                                                                            )
                                                           )
                                         )
                                       )
                                ),
                                column(9,
                                       tabsetPanel(id="heatmap_tabset",
                                                   tabPanel(title="Static Heatmap Layout 1",
                                                            br(),
                                                            actionButton("pheatmap2", "Save to output"),
                                                            actionButton("heatmap_gct", "Save GCT data file to output"),
                                                            actionButton("plot_heatmap", "Plot/Refresh", style="color: #0961E3; background-color: #F6E98C ; border-color: #2e6da4"),
                                                            uiOutput("plot.heatmap")),
                                                   tabPanel(title="Static Heatmap Layout 2", 
                                                            br(),
                                                            actionButton("staticheatmap", "Save to output"), 
                                                            plotOutput("staticheatmap", height = 800)),
                                                   # tabPanel(title="Interactive Heatmap",textOutput("text"), p(), plotlyOutput("interactiveheatmap", height = 800)),
                                                   tabPanel(title="Interactive Heatmap",
                                                            # actionButton("action_heatmaps", "Plot/Refresh", style="color: #0961E3; background-color: #F6E98C ; border-color: #2e6da4"),
                                                            actionButton("action_heatmaps", "Plot/Refresh",
                                                                         style="color: #0961E3 !important; background-color: #F6E98C !important; border-color: #2e6da4 !important;"),  
                                                            p(),
                                                            tags$details(
                                                              tags$summary(
                                                                tags$strong("Click to view Interactive Heatmap Guide",
                                                                            style = "color: #28a745; cursor: pointer;")
                                                              ),
                                                              tags$div(
                                                                style = "padding: 15px; background-color: #f4faf6; border-left: 5px solid #28a745; margin-top: 10px;",
                                                                tags$p("This interactive heatmap is powered by ", tags$strong("Morpheus"), ", a fast, canvas-based heatmap viewer built for exploring large expression matrices — zoom, pan, reorder, and inspect values directly in the browser."),
                                                                tags$p(
                                                                  "For a full walkthrough of Morpheus's features (row/column reordering, search, filtering, annotation tracks, and more), see the official tutorial: ",
                                                                  tags$a(href = "https://software.broadinstitute.org/morpheus/tutorial.html", target = "_blank",
                                                                         "Tutorial to explore Morpheus functions")
                                                                )
                                                              )
                                                            ),
                                                            p(),
                                                            # uiOutput("interactiveheatmap_ui")),
                                                            morpheus::morpheusOutput("interactiveheatmap", width="1500px", height="1200px")),
                                                   tabPanel(title="Help", htmlOutput('help_heatmap'))
                                       )
                                )
                              )
                     ),
                     
                     
                     ##########################################################################################################
                     ## Volcano Plot
                     ##########################################################################################################
                     tabPanel("DEGs",
                              fluidRow(
                                column(3,
                                       wellPanel(
                                         column(width=12,uiOutput("selectGroupSampleDEG")),
                                         conditionalPanel( "input.volcano_tabset!='DEGs in Two Comparisons' && input.volcano_tabset!='DEG Counts'",
                                                           selectInput("volcano_test", label="Select Comparison Groups for Volcano Plot", choices=NULL)),
                                         conditionalPanel( "input.volcano_tabset=='DEGs in Two Comparisons'",
                                                           selectInput("volcano_test1", label="1st Comparison (X-axis)", choices=NULL),
                                                           selectInput("volcano_test2", label="2nd Comparison (Y-aixs)", choices=NULL)),
                                         numericInput("volcano_FCcut", label= "Choose Fold Change Cutoff",  value = 1.2, min=1, step=0.1),
                                         numericInput("volcano_pvalcut", label= "Choose P Value Cutoff", value=0.01, min=0, step=0.001),
                                         radioButtons("volcano_psel", label= "P value or P.adj Value?", choices= c("Pval"="Pval","Padj"="Padj"),inline = TRUE),
                                         conditionalPanel( "input.volcano_tabset!='DEGs in Two Comparisons' && input.volcano_tabset!='DEG Counts'",
                                                           textOutput("volcano_filteredgene"),
                                                           tags$head(tags$style("#volcano_filteredgene{color: red; font-size: 17px; font-style: italic; }" )),
                                                           textOutput("volcano_filteredgene2"),
                                                           tags$head(tags$style("#volcano_filteredgene2{color: red; font-size: 15px; font-style: italic; }" ))),
                                         conditionalPanel( "input.volcano_tabset!='Volcano Plot (Interactive)' && input.volcano_tabset!='DEG Counts'",
                                                           radioButtons("volcano_genelabel",label="Select Gene Label",inline = TRUE, choices=""),
                                                           radioButtons("volcano_label", label="Label Genes:", inline = TRUE, choices = c("DEGs","None", "Upload", "Geneset"), selected = "DEGs"),
                                                           conditionalPanel("input.volcano_label!='None'",
                                                                            sliderInput("Ngenes", "# of Genes to Label", min = 10, max = 200, step = 5, value = 50),
                                                                            colourInput("volcano_subset_color", "Gene Label Color:", "blue")
                                                           ),
                                                           conditionalPanel("input.volcano_label=='Upload'",
                                                                            textAreaInput("volcano_gene_list", "List of genes to label\n(UniqueID, Gene.Name or Protein.ID)", "", cols = 5, rows=6)
                                                           ),
                                                           conditionalPanel("input.volcano_label=='Geneset'",
                                                                            uiOutput("html_geneset")
                                                           ),
                                                           tags$hr(),
                                                           radioButtons("volcano_subset_highlight", label="Highlight a Set of Genes in a Different Color?",
                                                                        inline = TRUE, choices = c("No","Yes"), selected = "No"),
                                                           conditionalPanel("input.volcano_subset_highlight=='Yes'",
                                                                            textAreaInput("volcano_subset_gene_list",
                                                                                          "Enter Gene List to Highlight\n(UniqueID, Gene.Name or Protein.ID)",
                                                                                          "", cols = 5, rows=6),
                                                                            colourInput("volcano_subset_highlight_color", "Color of Highlight Gene List and Label:", "black"),
                                                                            textAreaInput("volcano_subset_label_list",
                                                                                          "Enter Genes From the Highlight List to Label\n(only genes also present in the list above will be labeled; genes entered here are always labeled, never dropped)",
                                                                                          "", cols = 5, rows=6),
                                                                            actionButton("apply_highlight", "Apply Highlight",
                                                                                         style="color: #0961E3; background-color: #F6E98C; border-color: #2e6da4")
                                                           )
                                         ),
                                         conditionalPanel("input.volcano_tabset!='DEG Counts'",
                                                          radioButtons("more_options", label="Show More Options", inline = TRUE, choices = c("Yes","No"), selected = "No"),
                                                          conditionalPanel("input.more_options=='Yes'",
                                                                           numericInput("Max_logFC", label= "Max abs(log2FC) in plot (use 0 for full range)", value=0, min=0),
                                                                           numericInput("Max_Pvalue", label= "Max -log10(Stat Value) in plot (use 0 for full range)", value=0, min=0),
                                                                           conditionalPanel( "input.volcano_tabset!='Volcano Plot (Interactive)'",
                                                                                             sliderInput("lfontsize", "Label Font Size:", min = 1, max = 10, step = 1, value = 4),
                                                                                             sliderInput("yfontsize", "Legend Font Size:", min = 8, max = 24, step = 1, value = 14),
                                                                                             radioButtons("vlegendpos", label="Legend position", inline = TRUE, choices = c("bottom","right"), selected = "bottom"),
                                                                                             radioButtons("rasterize", label="Rasterize plot", inline = TRUE, choices = c("Yes","No"), selected = "No")
                                                                           ),
                                                                           conditionalPanel( "input.volcano_tabset=='DEGs in Two Comparisons'",
                                                                                             radioButtons("DEG_comp_XY", label="Make X and Y scale the same?", inline = TRUE, choices = c("Yes","No"), selected = "Yes"),
                                                                                             radioButtons("DEG_comp_color", label="Color DEGs using dot color?", inline = TRUE, choices = c("Yes","No"), selected = "Yes")))
                                                          
                                         )
                                       )
                                ),
                                column(9,
                                       tabsetPanel(id="volcano_tabset",
                                                   tabPanel(title="DEG Counts",
                                                            tags$p("Click a comparison name to view volcano plot."), tags$hr(),
                                                            DT::dataTableOutput("deg_counts")),
                                                   tabPanel(title="Volcano Plot (Static)",actionButton("volcano", "Save to output"), plotOutput("volcanoplotstatic", height=800)),
                                                   tabPanel(title="Volcano Plot (Interactive)",
                                                            hr(),
                                                            p("Use the lasso or box-select tool in the plot's toolbar (top-right of the plot) to circle a region and see the genes it contains in the table below. Double-click the plot to clear your selection."),
                                                            plotlyOutput("volcanoplot", height = 800),
                                                            hr(),
                                                            fluidRow(
                                                              column(12, textOutput("volcano_selected_count"))
                                                            ),
                                                            fluidRow(
                                                              column(12, actionButton("volcano_selected_to_highlight", "Send Selected Genes to Highlight List",
                                                                                      style="color: #0961E3; background-color: #F6E98C; border-color: #2e6da4"))
                                                            ),
                                                            hr(),
                                                            DT::dataTableOutput("volcano_selected_table")
                                                   ),
                                                   tabPanel(title="DEGs in Two Comparisons",actionButton("DEG_comp", "Save to output"), plotOutput("DEG_Compare", height=800)),
                                                   tabPanel(title="Data Table", actionButton("DEG_data", "Save to output"), DT::dataTableOutput("volcanoData")),
                                                   tabPanel(title="Help", htmlOutput('help_volcano'))
                                       )
                                )
                              )
                     ),
                     
                     ##########################################################################################################
                     ## Gene Set Enrichment
                     ##########################################################################################################
                     #tabPanel("Gene Set Enrichment", geneset_ui(id = "GS")), 
                     
                     
                     ##########################################################################################################
                     ## Pattern Clustering
                     ##########################################################################################################
                     ## 10/08/2020
                     ## eidted by bgao, add gene upload box and more plotting options
                     tabPanel("Pattern Clustering", value = 'Pattern_Clustering',
                              fluidRow(
                                column(3,
                                       wellPanel(
                                         column(width=12,uiOutput("selectGroupSamplePattern")),
                                         radioButtons("pattern_subset",label="Use subset genes or upload your own subset?", choices=c("subset","upload genes"),inline = TRUE, selected="subset"),
                                         conditionalPanel("input.pattern_subset=='subset'",
                                                          selectInput("pattern_test", label="Select Genes from Test:", choices=NULL),
                                                          column(width=6,numericInput("pattern_fccut", label= "Choose Fold Change Threshold",value = 1.2, min=1, step=0.1)),
                                                          column(width=6,numericInput("pattern_pvalcut", label= "Choose P-value Threshold", value=0.01, min=0, step=0.001)),
                                                          radioButtons("pattern_psel", label= "P value or P.adj Value?", choices= c("Pval"="Pval","Padj"="Padj"),inline = TRUE),
                                                          textOutput("patternfilteredgene")
                                         ),
                                         conditionalPanel("input.pattern_subset=='upload genes'", textAreaInput("pattern_list", "list", "", cols = 5, rows=6)),
                                         
                                         tags$head(tags$style("#patternfilteredgene{color: red; font-size: 20px; font-style: italic; }")),
                                         tags$br(),
                                         selectInput("pattern_attr", label="Select Pattern Variable", choices=NULL),
                                         # selectizeInput("pattern_group", label="Select Groups (re-order under Groups and Samples tab)", choices=NULL, multiple=TRUE),
                                         radioButtons("ClusterMethod", label="Cluster Method", inline = FALSE, choices = c("Soft Clustering" = "mfuzz", "K-means" = "kmeans")),
                                         sliderInput("k", "Cluster Number:", min = 3, max = 12, step = 1, value = 6),
                                         conditionalPanel("input.ClusterMethod=='kmeans'",
                                                          sliderInput("pattern_font", "Font Size:", min = 12, max = 24, step = 1, value = 14),
                                                          sliderInput("pattern_Xangle", label= "X Angle", min = 0, max = 90, step = 15, value = 45)),
                                         sliderInput("pattern_ncol", label= "Column Number", min = 1, max = 6, step = 1, value = 3),
                                         conditionalPanel("input.Pattern_tabset=='Data Table'",
                                                          radioButtons("DataFormat", label="Data Output Format:", inline = TRUE, choices = c("Wide Format" = "wide","Long Format" = "long"))
                                         )
                                       )
                                ),
                                column(9,
                                       tabsetPanel(id="Pattern_tabset",
                                                   #tabPanel(title="Optimal Number of Clusters", plotOutput("nbclust", height=800)),
                                                   tabPanel(title="Clustering of Centroid Profiles",
                                                            uiOutput('ui_sel_order_group'),
                                                            uiOutput('reset_group'),
                                                            tags$br(),
                                                            actionButton("pattern_plot", "Plot", style="color: #0961E3; background-color: #F6E98C ; border-color: #2e6da4"),
                                                            tags$br(),tags$hr(),
                                                            actionButton("pattern", "Save to output"),
                                                            plotOutput("pattern", height=800)),
                                                   #tabPanel(title="Time Series Plot",actionButton("Pattern_data", "Save to output"), DT::dataTableOutput("dat_pattern")),
                                                   tabPanel(title="Data Table", DT::dataTableOutput("dat_pattern")),
                                                   tabPanel(title="Help", htmlOutput('help_pattern'))
                                       )
                                )
                              )
                     ),
                     
                     ##########################################################################################################
                     ## Correlation Network
                     ##########################################################################################################
                     tabPanel("Correlation Network", value = 'Correlation_Network',
                              fluidRow(
                                column(3,
                                       wellPanel(
                                         column(width=12,uiOutput("selectGroupSampleNetwork")),
                                         radioButtons("network_label",label="Select Gene Label",inline = TRUE, choices=c("UniqueID", "Gene.Name"), selected="Gene.Name"),
                                         selectizeInput("sel_net_gene",	label="Gene Name (Select 1 or more)",	choices = NULL,	multiple=TRUE, options = list(placeholder =	'Type to search')),
                                         sliderInput("network_rcut", label= "Choose r Cutoff",  min = 0, max = 1, value = 0.9, step=0.02),
                                         selectInput("network_pcut", label= "Choose P Value Cutoff", choices= c("0.0001"=0.0001,"0.001"=0.001,"0.01"=0.01,"0.05"=0.05),selected=0.01),
                                         textOutput("networkstat"),
                                         uiOutput("myTabUI")
                                       )
                                ),
                                column(9,
                                       tabsetPanel(id="Network_tabset",
                                                   tabPanel("visNetwork", visNetworkOutput("visnetwork", height="800px"), style = "background-color: #eeeeee;"),
                                                   tabPanel("networkD3", forceNetworkOutput("networkD3", height="800px"), style = "background-color: #eeeeee;"),
                                                   tabPanel(title="Data Table",	DT::dataTableOutput("dat_network")),
                                                   tabPanel(title="Help", htmlOutput('help_network'))
                                       )
                                )
                              )
                     ),
                     
                     
                     
                     ##########################################################################################################
                     ## Venn Diagram
                     ##########################################################################################################
                     tabPanel("Venn Diagram",
                              fluidRow(
                                column(3,
                                       wellPanel(
                                         conditionalPanel("input.venn_combined=='Venn Diagram from Current Project'",
                                                          numericInput("venn_fccut", label= "Choose Fold Change Cutoff", value = 1.2, min=1, step=0.1),
                                                          numericInput("venn_pvalcut", label= "Choose P value Cutoff", value=0.01, min=0, step=0.001),
                                                          radioButtons("venn_psel", label= "P value or P.adj Value?", choices= c("Pval"="Pval","Padj"="Padj"),inline = TRUE),
                                                          radioButtons("venn_updown", label= "All, Up or Down?", choices= c("All"="All","Up"="Up","Down"="Down"),inline = TRUE),
                                                          selectInput("venn_test1", label="Select List 1", choices=NULL),
                                                          conditionalPanel("input.venn_tabset=='Venn Diagram'",
                                                                           colourInput("col1", "Select colour", "#0000FF",palette = "limited")
                                                          ),
                                                          
                                                          selectInput("venn_test2", label="Select List 2", choices=NULL),
                                                          conditionalPanel("input.venn_tabset=='Venn Diagram'",
                                                                           colourInput("col2", "Select colour", "#FF7F00",palette = "limited")
                                                          ),
                                                          
                                                          selectInput("venn_test3", label="Select List 3", choices=NULL),
                                                          conditionalPanel("input.venn_tabset=='Venn Diagram'",
                                                                           colourInput("col3", "Select colour", "#00FF00",palette = "limited")
                                                          ),
                                                          
                                                          selectInput("venn_test4", label="Select List 4", choices=NULL),
                                                          conditionalPanel("input.venn_tabset=='Venn Diagram'",
                                                                           colourInput("col4", "Select colour", "#FF00FF",palette = "limited")
                                                          ),
                                                          
                                                          selectInput("venn_test5", label="Select List 5", choices=NULL),
                                                          conditionalPanel("input.venn_tabset=='Venn Diagram'",
                                                                           colourInput("col5", "Select colour", "#FFFF00", palette = "limited")
                                                          ),
                                                          conditionalPanel("input.venn_tabset=='Intersection Output'",
                                                                           radioButtons("vennlistname", label= "Label name", choices= c("Gene.Name"="Gene.Name","AC Number"="AC", "UniqueID"="UniqueID"),inline = TRUE, selected = "Gene.Name"))
                                         ),
                                         conditionalPanel("input.venn_combined=='Venn Diagram Across Projects'",
                                                          numericInput("vennP_fccut", label= "Choose Fold Change Cutoff", value = 1.2, min=1, step=0.1),
                                                          numericInput("vennP_pvalcut", label= "Choose P value Cutoff", value=0.01, min=0, step=0.001),
                                                          radioButtons("vennP_psel", label= "P value or P.adj Value?", choices= c("Pval"="Pval","Padj"="Padj"),inline = TRUE),
                                                          checkboxInput("upperSymbols", "Upper Case Gene Symbols (e.g. mouse vs human)?",  FALSE, width="90%"),
                                                          selectInput("dataset1", "Data set1", choices=NULL),
                                                          selectInput("vennP_test1", label="Select List 1", choices=NULL),
                                                          conditionalPanel("input.vennP_tabset=='Venn Diagram'",
                                                                           colourInput("vennPcol1", "Select colour", "#0000FF",palette = "limited")
                                                          ),
                                                          selectInput("dataset2", "Data set2", choices=NULL),
                                                          selectInput("vennP_test2", label="Select List 2", choices=NULL),
                                                          conditionalPanel("input.vennP_tabset=='Venn Diagram'",
                                                                           colourInput("vennPcol2", "Select colour", "#FF7F00",palette = "limited")
                                                          ),
                                                          selectInput("dataset3", "Data set3", choices=NULL),
                                                          selectInput("vennP_test3", label="Select List 3", choices=NULL),
                                                          conditionalPanel("input.vennP_tabset=='Venn Diagram'",
                                                                           colourInput("vennPcol3", "Select colour", "#00FF00",palette = "limited")
                                                          ),
                                                          selectInput("dataset4", "Data set4", choices=NULL),
                                                          selectInput("vennP_test4", label="Select List 4", choices=NULL),
                                                          conditionalPanel("input.vennP_tabset=='Venn Diagram'",
                                                                           colourInput("vennPcol4", "Select colour", "#FF00FF",palette = "limited")
                                                          ),
                                                          selectInput("dataset5", "Data set5", choices=NULL),
                                                          selectInput("vennP_test5", label="Select List 5", choices=NULL),
                                                          conditionalPanel("input.vennP_tabset=='Venn Diagram'",
                                                                           colourInput("vennPcol5", "Select colour", "#FFFF00", palette = "limited")
                                                          )
                                         )
                                       )
                                ),
                                column(9,
                                       tabsetPanel(id="venn_combined",
                                                   tabPanel(title="Venn Diagram from Current Project",
                                                            tabsetPanel(id="venn_tabset",
                                                                        tabPanel(title="Venn Diagram",
                                                                                 column(9,
                                                                                        actionButton("vennDiagram", "Save to output"),plotOutput("vennDiagram", height = 800)
                                                                                 ),
                                                                                 column(3,
                                                                                        textInput("title", "Title", width = "100%"),
                                                                                        sliderInput("maincex", "Title Size", min = 0, max = 6, value = 3, width = "100%"),
                                                                                        sliderInput("alpha", "Opacity", min = 0, max = 1, value = 0.4, width = "100%"),
                                                                                        sliderInput("lwd", "Line Width", min = 1, max = 4, value = 1, width = "100%"),
                                                                                        sliderInput("lty", "Line Type", min = 1, max = 6, value = 1, width = "100%"),
                                                                                        radioButtons("fontface","Number Font face",list("plain", "bold", "italic"),selected = "plain",	inline = TRUE),
                                                                                        sliderInput("cex", "Font size", min = 1, max = 4, value = 2, width = "100%"),
                                                                                        radioButtons("catfontface","Label Font face",list("plain", "bold", "italic"),	selected = "plain", inline = TRUE),
                                                                                        sliderInput("catcex", "Font size", min = 1, max = 2, step=0.1, value = 1.8, width = "100%"),
                                                                                        sliderInput("catwrap", "Wrap Labels At (characters)", min = 5, max = 40, step = 1, value = 18, width = "100%"),
                                                                                        sliderInput("margin", "Margin", min = 0, max = 1, step=0.05, value = 0.1, width = "100%")
                                                                                 )
                                                                        ),
                                                                        tabPanel(title="Venn Diagram (black & white)", plotOutput("SvennDiagram",height = 800, width = 800)),
                                                                        tabPanel(title="Intersection Output",
                                                                                 rclipboard::rclipboardSetup(),
                                                                                 br(),
                                                                                 helpText("Click one or more rows below to select a gene list."),
                                                                                 uiOutput("venn_copy_btn"),
                                                                                 br(),
                                                                                 DT::dataTableOutput("venn_intersect_table")
                                                                        ),
                                                                        tabPanel(title="DEG Table", actionButton("venn_DEG_data", "Save to output"), DT::dataTableOutput("venn_DEG_Data")),
                                                                        tabPanel(title="Help", htmlOutput('help_venn'))
                                                            )
                                                   ),
                                                   tabPanel(title="Venn Diagram Across Projects",
                                                            tabsetPanel(id="vennP_tabset",
                                                                        tabPanel(title="Venn Diagram",
                                                                                 column(9,
                                                                                        plotOutput("vennPDiagram", height = 800)
                                                                                 ),
                                                                                 column(3,
                                                                                        textInput("vennPtitle", "Title", width = "100%"),
                                                                                        sliderInput("vennPmaincex", "Title Size", min = 0, max = 6, value = 3, width = "100%"),
                                                                                        sliderInput("vennPalpha", "Opacity", min = 0, max = 1, value = 0.4, width = "100%"),
                                                                                        sliderInput("vennPlwd", "Line Width", min = 1, max = 4, value = 1, width = "100%"),
                                                                                        sliderInput("vennPlty", "Line Type", min = 1, max = 6, value = 1, width = "100%"),
                                                                                        radioButtons("vennPfontface",	"Number Font face", list("plain", "bold", "italic"),	selected = "plain",	inline = TRUE),
                                                                                        sliderInput("vennPcex", "Font size", min = 1, max = 4, value = 2, width = "100%"),
                                                                                        radioButtons("vennPcatfontface", "Label Font face", list("plain", "bold", "italic"),	selected = "plain",	inline = TRUE),
                                                                                        sliderInput("vennPcatcex", "Font size", min = 1, max = 2, step=0.1, value = 1.8, width = "100%"),
                                                                                        sliderInput("vennPcatwrap", "Wrap Labels At (characters)", min = 5, max = 40, step = 1, value = 18, width = "100%"),
                                                                                        sliderInput("vennPmargin", "Margin", min = 0, max = 1, step=0.05, value = 0.2, width = "100%")
                                                                                 )
                                                                        ),
                                                                        tabPanel(title="Venn Diagram (black & white)", plotOutput("SvennPDiagram",height = 800,width = 800)),
                                                                        tabPanel(title="Intersection Output",
                                                                                 rclipboard::rclipboardSetup(),
                                                                                 br(),
                                                                                 helpText("Click one or more rows below to select a gene list."),
                                                                                 uiOutput("vennP_copy_btn"),
                                                                                 br(),
                                                                                 DT::dataTableOutput("vennP_intersect_table")
                                                                        ),
                                                                        tabPanel(title="Help", htmlOutput('help_vennp'))
                                                                        
                                                            )
                                                   )
                                       )
                                )
                              )
                     ),
                     
                     
                     ##########################################################################################################
                     ## Output
                     ##########################################################################################################
                     
                     tabPanel("Output",
                              sliderInput("pdf_width", "Plot File Page Width", min = 3, max = 30, step = 1, value = 12),
                              sliderInput("pdf_height", "Plot File Page Height", min = 3, max = 50, step = 1, value = 8),
                              actionButton("clear_saved_plots", "Clear all saved plots"),
                              tags$br(),
                              #  htmlOutput("saved_plot_list"),
                              checkboxGroupInput("plots_checked", "Plots to Save", choices=NULL, selected=NULL),
                              downloadButton('downloadPDF', 'Download PDF'),
                              downloadButton('downloadSVG', 'Download SVG (for the first selected plot)'),
                              tags$br(),tags$hr(),
                              downloadButton('downloadXLSX', 'Download tables in .xlsx'),
                              tags$br(),tags$hr(),
                              checkboxGroupInput("GCT_table_checked", "GCT tables to Save", choices=NULL, selected=NULL),
                              downloadButton('downloadGCT', 'Download tables in .gct'),
                              tags$br(),tags$hr(),
                              h4("Save Session"),
                              tags$p("Save your current settings (filters, comparisons, plot options, across all tabs) as a file -- upload it later to restore this exact state."),
                              actionButton("save_session_btn", "Save Session", icon = icon("save", lib = "glyphicon")),
                              uiOutput("save_session_url"),
                              tags$hr(),
                              tags$p("Or restore settings previously saved to a file:"),
                              fileInput("upload_session_file", "Restore Session from File", accept = ".rds", width = "100%")
                     ),
                     
                     
                     ##########################################################################################################
                     ## footer
                     ##########################################################################################################
                     
                     footer= HTML(footer_text)
                     
                     
          )
  )
) #for tagList
} #for ui function(request)

server <- function(input, output, session) {
  # actionButton/downloadButton click counts get bookmarked like any other
  # input -- restoring a non-zero count silently re-triggers that button's
  # own observeEvent right after restore, since Shiny can't distinguish
  # "value changed because it was restored" from "value changed because it
  # was clicked". Confirmed directly: restoring save_session_btn's own
  # saved (non-zero) click count re-triggered ANOTHER session$doBookmark()
  # immediately after the first restore, with no click involved. None of
  # these buttons' click counts are meaningful state to restore -- their
  # actual effects are already captured via whatever reactiveValues/state
  # they produced, not via the click count itself -- so exclude them all.
  # Module-scoped buttons (correlation.R/wgcna.R/genesetmodule.R/
  # TimeSeries.R) are excluded separately, inside each module's own
  # moduleServer(), using that module's own bare (un-namespaced) IDs.
  session$setBookmarkExclude(c(
    "save_session_btn",
    # Every fileInput's value (name/size/type/datapath) gets bookmarked like
    # any other input. Restoring one is never safe: the saved datapath
    # always points at a per-session temp file (always contains "/"), and
    # Shiny's OWN restore code explicitly rejects any file input value whose
    # datapath contains "/" -- "Invalid '/' found in file input path" -- an
    # UNCAUGHT error that breaks the whole session (confirmed directly: this
    # is what actually produced a stuck/greyed-out "dead loop" after
    # restoring a session that had a file upload in it, not a reactive
    # loop in groupandsample.R). For upload_session_file specifically there
    # was a second-order issue too: a restored non-NULL value re-triggers
    # its own observeEvent below, generating yet another bookmark and
    # redirecting again -- same root cause as the actionButton
    # auto-retrigger issue elsewhere in this file, compounding the crash
    # into a genuine repeat-forever loop.
    "upload_session_file", "file1", "file2", "sd_cor_annot_color_file",
    "file_gene_highlight", "file_gene_annot", "annot_color_file",
    "F_sample", "F_exp", "F_comp", "F_annot",
    "DEG_comp", "DEG_data", "Dendrograms", "PCA_refresh_sample", "Pattern_data", "ProteinGeneName",
    "QCboxplot", "SampleCorrelation", "SampleDistance", "action_heatmaps",
    "alignQC_OV", "alignQC_RA", "alignQC_TGL", "alignQC_TGR", "apply_highlight",
    "boxplot", "browsing", "clear_saved_plots", "compute_PC", "covar_cat", "covar_num",
    "data_wide", "downloadGCT", "downloadPDF", "downloadSVG", "downloadXLSX",
    "heatmap_gct", "heatmap_test2sample", "histplot", "pattern", "pattern_plot",
    "pcaplot", "pheatmap2", "plot_PCA", "plot_browsing", "plot_exp", "plot_heatmap",
    "reset_all_types", "reset_group", "results", "sample", "staticheatmap",
    "vennDiagram", "venn_DEG_data", "volcano", "volcano_selected_to_highlight",
    "gennet", "customData", "uploadData"
  ))

  output$dynamic_sidebar_css <- renderUI({
    main_pct <- 100 - input$sidebar_width_pct
    tags$style(HTML(sprintf(
      ".col-sm-3 { width: %d%%; } .col-sm-9 { width: %d%%; }",
      input$sidebar_width_pct, main_pct
    )))
  })

  # Some restored inputs' *choices* are only ever populated lazily, when the
  # user actually visits the tab that computes them -- e.g. sel_net_gene on
  # Correlation Network: NetworkReactive() is gated on
  # input$menu=="Correlation_Network" and can be an expensive correlation
  # computation ("may take a few minutes" per its own progress message), so
  # it's deliberately NOT forced to run eagerly just because a session is
  # being restored. A fixed-duration poll (like the retry mechanism below)
  # is useless here since there's no bound on how long until the user visits
  # that tab. Instead, stash the restored value here; network.R's own
  # observeEvent consumes and clears it the first time it actually computes
  # choices, whenever that happens to be.
  restored_sel_net_gene <- reactiveVal(NULL)

  # Pattern Clustering's group_source/group_dest (Groups to Plot / Drag Here
  # to Remove) are shinyjqui::orderInput widgets that get completely
  # regenerated by pattern.R's own renderUI every time input$pattern_attr
  # changes -- same fixed IDs are reused for whichever attribute is
  # currently selected, so a raw restored value only makes sense once
  # pattern_attr has itself been restored back to the matching attribute.
  # Stashed here and consumed by pattern.R's own renderUI once that lines up.
  restored_pattern_group_snapshot <- reactiveVal(NULL)

  # State ID of the most recent session bookmark, so the download handler
  # below can serve that exact bookmark's saved input.rds as a file.
  last_bookmark_state_id <- reactiveVal(NULL)

  ##########################################################################################################
  ## Save Session (first pass): Shiny's built-in server-side bookmarking captures every plain input$-bound
  ## widget's current value and saves it server-side under a short state ID; the returned URL restores them
  ## all on load. Known gaps for a later pass: renderUI-created controls that don't exist yet at restore
  ## time, the custom drag-and-drop Groups-and-Samples UI (multidrag.js, not a standard Shiny input), and
  ## this session's DT/plotly selection state (lives in plain reactiveVals, not input$) -- none of those are
  ## covered by plain bookmarking and would need explicit onBookmark/onRestore handling if wanted later.
  ##########################################################################################################
  observeEvent(input$save_session_btn, {
    session$doBookmark()
  })

  onBookmarked(function(url) {
    # The state ID is the last query-string segment Shiny appended to the
    # current URL, e.g. "...?_state_id_=6ba545c11e3c6206".
    state_id <- sub(".*_state_id_=", "", url)
    last_bookmark_state_id(state_id)
    output$save_session_url <- renderUI({
      downloadButton("download_session_file", "Download Session File", class = "btn-default btn-sm")
    })
    showNotification("Session saved -- click the button below to download it.", type = "message")
  })

  # Serves the exact input.rds Shiny already wrote for the bookmark above --
  # same values, same exclusions (action buttons etc.), just handed to the
  # user as a portable file instead of (or alongside) the link. A link only
  # keeps working as long as this exact server retains that state ID's
  # folder on disk; the file survives a redeploy, a disk cleanup, or moving
  # to a different server, and can be emailed/archived like any other file.
  output$download_session_file <- downloadHandler(
    filename = function() {
      paste0("Quickomics_session_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds")
    },
    content = function(file) {
      state_id <- last_bookmark_state_id()
      req(state_id)
      src <- file.path("shiny_bookmarks", state_id, "input.rds")
      req(file.exists(src))
      file.copy(src, file, overwrite = TRUE)
    }
  )

  # Restoring from an uploaded file reuses the exact same restore pipeline
  # as the URL: write the uploaded input.rds into a freshly-generated
  # bookmark folder (same on-disk layout Shiny's own doBookmark() uses),
  # then navigate the browser to that state ID's URL. Every onRestored()
  # handler in this app (including the polling ones in modules) fires
  # exactly as it would for a normal saved-link restore.
  observeEvent(input$upload_session_file, {
    req(input$upload_session_file)
    uploaded <- tryCatch(readRDS(input$upload_session_file$datapath), error = function(e) NULL)
    if (is.null(uploaded) || !is.list(uploaded)) {
      showNotification("That file doesn't look like a Quickomics session file.", type = "error", duration = NULL)
      return()
    }
    # Defensively strip any fileInput value the uploaded file might already
    # contain (e.g. saved before this exclusion existed, or from a session
    # that had one of these set at save time) -- restoring any of them
    # crashes with "Invalid '/' found in file input path" (see the exclude
    # list above for why). Belt-and-suspenders: setBookmarkExclude only
    # stops *new* saves from including these; this stops a *previously*
    # saved file from ever reintroducing one.
    for (fk in c("upload_session_file", "file1", "file2", "sd_cor_annot_color_file",
                 "file_gene_highlight", "file_gene_annot", "annot_color_file",
                 "F_sample", "F_exp", "F_comp", "F_annot", "GS-custom_gmt_file")) {
      uploaded[[fk]] <- NULL
    }
    new_state_id <- paste(sample(c(letters[1:6], 0:9), 16, replace = TRUE), collapse = "")
    dest_dir <- file.path("shiny_bookmarks", new_state_id)
    dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(uploaded, file.path(dest_dir, "input.rds"))
    showNotification("Session file loaded -- restoring...", type = "message")
    shinyjs::runjs(sprintf(
      "window.location.href = window.location.pathname + '?_state_id_=%s';",
      new_state_id
    ))
  })

  onRestored(function(state) {
    showNotification("Session restored from saved link.", type = "message")

    # The active top-level tab (input$menu) needs an explicit updateTabsetPanel()
    # rather than relying on Shiny's normal input restoration. Confirmed directly:
    # 4 of these tabs (Gene Set Enrichment/gsea, WGCNA/wgcna, Correlation
    # Analysis/Correlation, Time Course Analysis/time_series) are added via
    # insertTab() at server runtime rather than being part of the static
    # ui(request) output -- they don't exist in the DOM yet when Shiny's normal
    # bookmark restoration tries to apply input$menu, so it silently falls back
    # to the first tab instead. Re-applying it here works because insertTab()
    # has already run (synchronously, earlier in this same server() call) by
    # the time onRestored() fires.
    if (!is.null(state$input$menu)) {
      updateTabsetPanel(session, "menu", selected = state$input$menu)
    }

    # sel_net_gene (Correlation Network) -- see restored_sel_net_gene's
    # definition above for why this is deferred rather than polled.
    if (!is.null(state$input$sel_net_gene)) {
      restored_sel_net_gene(state$input$sel_net_gene)
    }

    # group_source/group_dest (Pattern Clustering) -- see
    # restored_pattern_group_snapshot's definition above.
    if (!is.null(state$input$group_source) || !is.null(state$input$group_dest)) {
      restored_pattern_group_snapshot(list(
        attr   = state$input$pattern_attr,
        source = state$input$group_source,
        dest   = state$input$group_dest
      ))
    }

    # Expression Plot's colorby/plotx/sel_geneid/expression_test/sel_page/
    # sel_gene (barboxplot.R) get their *choices* populated by server-side
    # observe() blocks that only run once project data has loaded -- at
    # restore time those haven't run yet, so these inputs' own "preserve
    # previous selection" isolate(input$X) check reads NULL and falls back
    # to a hardcoded default, clobbering the restored value before we ever
    # see it. sel_page's observer also *reactively* depends on
    # expression_test/expression_fccut/expression_pvalcut/numperpage, so it
    # re-fires (recomputing its own choices from scratch) again after THOSE
    # get restored.
    #
    # A single post-flush reapply (session$onFlushed(fn, once=TRUE)) is not
    # good enough for any of these: onFlushed(once=TRUE) only runs on the
    # NEXT flush that happens to occur, and if nothing else in the app
    # triggers one, it may simply never fire at all (confirmed directly --
    # for sel_gene specifically, a bare onFlushed(once=TRUE) reapply never
    # ran, not even once, in an otherwise-idle session). Poll instead:
    # invalidateLater() *guarantees* its own recurring flush cycles, so the
    # reapply is retried every ~300ms (up to ~6s) until each restored value
    # actually sticks, regardless of what else is or isn't happening in the
    # reactive graph. In practice this resolves within 1-2 attempts.
    restored <- state$input
    pending <- list(
      colorby         = list(value = restored$colorby,         apply = function(v) updateSelectInput(session, "colorby", selected = v)),
      plotx           = list(value = restored$plotx,           apply = function(v) updateSelectInput(session, "plotx", selected = v)),
      sel_geneid      = list(value = restored$sel_geneid,       apply = function(v) updateRadioButtons(session, "sel_geneid", selected = v)),
      expression_test = list(value = restored$expression_test, apply = function(v) updateSelectizeInput(session, "expression_test", selected = v)),
      # Must re-supply the full page choices on every apply -- sel_page's own
      # populate observer (barboxplot.R) computes choices from the current
      # gene-count/numperpage and re-fires whenever expression_test etc.
      # change (including when restored above), so a bare selected= outside
      # its still-default choices=1 would silently fail to apply.
      sel_page = list(value = restored$sel_page, apply = function(v) {
        req(DataQCReactive())
        results_long <- DataQCReactive()$tmp_results_long
        req(results_long)
        expression_test <- isolate(input$expression_test)
        expression_fccut <- log2(as.numeric(isolate(input$expression_fccut)))
        expression_pvalcut <- as.numeric(isolate(input$expression_pvalcut))
        numperpage <- as.numeric(isolate(input$numperpage))
        req(expression_test, numperpage)
        if (identical(isolate(input$expression_psel), "Padj")) {
          filteredgene <- results_long %>% dplyr::filter(abs(logFC) > expression_fccut & Adj.P.Value < expression_pvalcut) %>% dplyr::filter(test == expression_test)
        } else {
          filteredgene <- results_long %>% dplyr::filter(abs(logFC) > expression_fccut & P.Value < expression_pvalcut) %>% dplyr::filter(test == expression_test)
        }
        page_choices <- seq_len(ceiling(nrow(filteredgene) / numperpage))
        updateSelectInput(session, "sel_page", choices = page_choices, selected = v)
      }),
      # Must re-supply the full choices list (DataIngenesReactive(), defined
      # in barboxplot.R) on every apply -- a server=TRUE selectize update
      # with no choices registers an empty searchable dataset, so the
      # restored gene could never be found/rendered even though `selected`
      # itself was set correctly.
      sel_gene        = list(value = restored$sel_gene,         apply = function(v) updateSelectizeInput(session, "sel_gene", choices = isolate(DataIngenesReactive()), selected = v, server = TRUE)),
      # Ylab/linear_base/linear_small_value (barboxplot.R's "linear value
      # parameters" observe(), ~line 103) get unconditionally overwritten
      # with a computed default every time exp_plot_Y_scale changes -- with
      # no "preserve previous value" guard at all (unlike colorby/plotx
      # etc.). Restoring exp_plot_Y_scale itself (a plain radio button,
      # which restores natively -- no retry needed for it specifically)
      # re-triggers that observer, clobbering a just-restored custom Ylab.
      Ylab               = list(value = restored$Ylab,               apply = function(v) updateTextInput(session, "Ylab", value = v)),
      linear_base        = list(value = restored$linear_base,        apply = function(v) updateTextInput(session, "linear_base", value = v)),
      linear_small_value = list(value = restored$linear_small_value, apply = function(v) updateTextInput(session, "linear_small_value", value = v)),

      # QC Plots (qcplot.R) -- PCAcolorby/PCAshapeby/PCAsizeby/PCA_label and
      # the Covariates tab's sd_cor_annotate_by/sd_cor_label_by all get their
      # choices+selected forced to a hardcoded default by the same
      # observeEvent(all_metadata(), ...) block, which also fires again at
      # restore time.
      PCAcolorby = list(value = restored$PCAcolorby, apply = function(v) {
        attrs <- sort(setdiff(colnames(all_metadata()), c("sampleid", "Order", "ComparePairs")))
        updateSelectInput(session, "PCAcolorby", choices = attrs, selected = v)
      }),
      PCAshapeby = list(value = restored$PCAshapeby, apply = function(v) {
        attrs <- sort(setdiff(colnames(all_metadata()), c("sampleid", "Order", "ComparePairs")))
        updateSelectInput(session, "PCAshapeby", choices = c("none", attrs), selected = v)
      }),
      PCAsizeby = list(value = restored$PCAsizeby, apply = function(v) {
        attrs <- sort(setdiff(colnames(all_metadata()), c("sampleid", "Order", "ComparePairs")))
        updateSelectInput(session, "PCAsizeby", choices = c("none", attrs), selected = v)
      }),
      PCA_label = list(value = restored$PCA_label, apply = function(v) {
        attrs <- sort(setdiff(colnames(all_metadata()), c("Order", "ComparePairs")))
        updateRadioButtons(session, "PCA_label", inline = TRUE, choices = attrs, selected = v)
      }),
      sd_cor_annotate_by = list(value = restored$sd_cor_annotate_by, apply = function(v) {
        attrs <- sort(setdiff(colnames(all_metadata()), c("sampleid", "Order", "ComparePairs")))
        updateSelectizeInput(session, "sd_cor_annotate_by", choices = attrs, selected = v)
      }),
      sd_cor_label_by = list(value = restored$sd_cor_label_by, apply = function(v) {
        attrs <- sort(setdiff(colnames(all_metadata()), c("Order", "ComparePairs")))
        updateSelectInput(session, "sd_cor_label_by", choices = attrs, selected = v)
      }),

      # Heatmap (heatmap.R) -- same "observe(){...} fires on every
      # all_metadata()/test_order() change, including at restore" pattern.
      heatmap_test = list(value = restored$heatmap_test, apply = function(v) {
        updateSelectizeInput(session, "heatmap_test", choices = test_order(), selected = v)
      }),
      heatmap_label = list(value = restored$heatmap_label, apply = function(v) {
        updateRadioButtons(session, "heatmap_label", inline = TRUE, choices = ProteinGeneNameHeader()[-1], selected = v)
      }),
      heatmap_annot = list(value = restored$heatmap_annot, apply = function(v) {
        attrs <- sort(setdiff(colnames(all_metadata()), c("sampleid", "Order", "ComparePairs")))
        updateSelectInput(session, "heatmap_annot", choices = attrs, selected = v)
      }),

      # Volcano Plot (volcano.R) -- same pattern, keyed off test_order().
      volcano_test = list(value = restored$volcano_test, apply = function(v) {
        updateSelectizeInput(session, "volcano_test", choices = test_order(), selected = v)
      }),
      volcano_test1 = list(value = restored$volcano_test1, apply = function(v) {
        updateSelectizeInput(session, "volcano_test1", choices = test_order(), selected = v)
      }),
      volcano_test2 = list(value = restored$volcano_test2, apply = function(v) {
        updateSelectizeInput(session, "volcano_test2", choices = test_order(), selected = v)
      }),
      volcano_genelabel = list(value = restored$volcano_genelabel, apply = function(v) {
        updateRadioButtons(session, "volcano_genelabel", inline = TRUE, choices = ProteinGeneNameHeader()[-1], selected = v)
      }),

      # Pattern Clustering (pattern.R) -- pattern_attr forced to "group" by
      # observeEvent(DataQCReactive(), ...); group_source/group_dest are
      # handled separately via restored_pattern_group_snapshot since they're
      # regenerated (not update*Input-able) and keyed off pattern_attr.
      pattern_attr = list(value = restored$pattern_attr, apply = function(v) {
        req(DataQCReactive())
        attrs <- sort(setdiff(colnames(DataQCReactive()$MetaData), c("sampleid", "Order", "ComparePairs")))
        updateSelectInput(session, "pattern_attr", choices = attrs, selected = v)
      }),
      pattern_test = list(value = restored$pattern_test, apply = function(v) {
        updateSelectizeInput(session, "pattern_test", choices = c("ALL", test_order()), selected = v)
      }),

      # AlignQC (alignQC.R) -- alignQC_var forced to a hardcoded default by
      # observe(){...} whenever DataQCReactive() changes, including restore.
      alignQC_var = list(value = restored$alignQC_var, apply = function(v) {
        req(DataQCReactive())
        MetaData <- DataQCReactive()$MetaData
        num_col <- colnames(dplyr::select_if(MetaData, is.numeric))
        req(length(num_col) > 0)
        updateSelectizeInput(session, "alignQC_var", choices = num_col, selected = v)
      }),

      # QC Plots (qcplot.R) -- PCA_list (List of Samples to Label) gets
      # unconditionally overwritten with the full current sample list by
      # observe(){ updateTextAreaInput(...) } every time sample_order()
      # changes, including at restore.
      PCA_list = list(value = restored$PCA_list, apply = function(v) {
        updateTextAreaInput(session, "PCA_list", value = v)
      }),

      # Volcano Plot (volcano.R) -- volcano_gene_list (List of genes to
      # label, Upload mode) gets unconditionally overwritten with an
      # auto-suggested DEG sample by observe(){...} every time volcano_test/
      # volcano_FCcut/volcano_pvalcut/Ngenes/DataQCReactive() change,
      # including at restore.
      volcano_gene_list = list(value = restored$volcano_gene_list, apply = function(v) {
        updateTextAreaInput(session, "volcano_gene_list", value = v)
      }),

      # Venn Diagram (venn.R) -- venn_test1..5, same "forced default on every
      # project load" pattern, sharing venn_test_choices() with the populate
      # observer.
      venn_test1 = list(value = restored$venn_test1, apply = function(v) {
        updateSelectizeInput(session, "venn_test1", choices = venn_test_choices(), selected = v)
      }),
      venn_test2 = list(value = restored$venn_test2, apply = function(v) {
        updateSelectizeInput(session, "venn_test2", choices = venn_test_choices(), selected = v)
      }),
      venn_test3 = list(value = restored$venn_test3, apply = function(v) {
        updateSelectizeInput(session, "venn_test3", choices = venn_test_choices(), selected = v)
      }),
      venn_test4 = list(value = restored$venn_test4, apply = function(v) {
        updateSelectizeInput(session, "venn_test4", choices = venn_test_choices(), selected = v)
      }),
      venn_test5 = list(value = restored$venn_test5, apply = function(v) {
        updateSelectizeInput(session, "venn_test5", choices = venn_test_choices(), selected = v)
      }),

      # Venn Diagram Across Projects (vennprojects.R) -- dataset1..5 pick
      # which project each slot compares, and vennP_test1..5's choices
      # depend on whichever project the matching datasetN currently points
      # to, so must be read fresh (isolate) on every reapply -- same
      # cascading pattern as sel_group/sel_attribute in correlation.R.
      dataset1 = list(value = restored$dataset1, apply = function(v) {
        updateSelectizeInput(session, "dataset1", choices = vennP_dataset_choices(), selected = v)
      }),
      dataset2 = list(value = restored$dataset2, apply = function(v) {
        updateSelectizeInput(session, "dataset2", choices = vennP_dataset_choices(), selected = v)
      }),
      dataset3 = list(value = restored$dataset3, apply = function(v) {
        updateSelectizeInput(session, "dataset3", choices = vennP_dataset_choices(), selected = v)
      }),
      dataset4 = list(value = restored$dataset4, apply = function(v) {
        updateSelectizeInput(session, "dataset4", choices = vennP_dataset_choices(), selected = v)
      }),
      dataset5 = list(value = restored$dataset5, apply = function(v) {
        updateSelectizeInput(session, "dataset5", choices = vennP_dataset_choices(), selected = v)
      }),
      vennP_test1 = list(value = restored$vennP_test1, apply = function(v) {
        tests <- vennP_test_choices_for(isolate(input$dataset1))
        req(tests)
        updateSelectizeInput(session, "vennP_test1", choices = tests, selected = v)
      }),
      vennP_test2 = list(value = restored$vennP_test2, apply = function(v) {
        tests <- vennP_test_choices_for(isolate(input$dataset2))
        req(tests)
        updateSelectizeInput(session, "vennP_test2", choices = tests, selected = v)
      }),
      vennP_test3 = list(value = restored$vennP_test3, apply = function(v) {
        tests <- vennP_test_choices_for(isolate(input$dataset3))
        req(tests)
        updateSelectizeInput(session, "vennP_test3", choices = tests, selected = v)
      }),
      vennP_test4 = list(value = restored$vennP_test4, apply = function(v) {
        tests <- vennP_test_choices_for(isolate(input$dataset4))
        req(tests)
        updateSelectizeInput(session, "vennP_test4", choices = tests, selected = v)
      }),
      vennP_test5 = list(value = restored$vennP_test5, apply = function(v) {
        tests <- vennP_test_choices_for(isolate(input$dataset5))
        req(tests)
        updateSelectizeInput(session, "vennP_test5", choices = tests, selected = v)
      })
    )
    pending <- Filter(function(p) !is.null(p$value), pending)

    if (length(pending) > 0) {
      # Always reapply every pending input on every tick -- deliberately NOT
      # stopping as soon as isolate(input[[nm]]) appears to match. Confirmed
      # directly (Ylab/exp_plot_Y_scale case) that this "stop on first
      # match" shortcut is unsafe: a clobbering observer's update message
      # can already have changed the value client-side (visually) before
      # that change round-trips back to update input[[nm]] server-side, so
      # a check right after applying our fix can read the OLD, still-
      # correct server value, declare victory, and stop -- while the client
      # has already moved on to the wrong one, with nothing left to correct
      # it afterward. Unconditionally reapplying for a fixed ~3s window
      # instead guarantees our value is re-asserted after any such stray
      # clobber, at the trivial cost of a few redundant no-op messages.
      attempts_left <- 10  # ~3s at 300ms
      restore_observer <- NULL
      restore_observer <- observe({
        invalidateLater(300, session)
        isolate({
          attempts_left <<- attempts_left - 1
          for (nm in names(pending)) {
            p <- pending[[nm]]
            # tryCatch so one input's apply() failing on an early attempt
            # (e.g. sel_gene's DataIngenesReactive() req()-ing out before
            # data has loaded yet) can't abort the rest of this tick's loop
            # -- it just gets retried again next tick like normal.
            tryCatch(p$apply(p$value), error = function(e) NULL)
          }
          if (attempts_left <= 0) {
            restore_observer$destroy()
          }
        })
      })
    }
  })


  source("inputdata.R",local = TRUE)
  source("process_uploaded_files.R",local = TRUE)
  source("groupandsample.R",local=TRUE)
  source("qcplot.R",local=TRUE)
  source("alignQC.R",local=TRUE)
  source("volcano.R",local = TRUE)
  source("heatmap.R",local = TRUE)
  source("barboxplot.R",local = TRUE)
  source("venn.R",local = TRUE)
  source("genesetmodule.R",local = TRUE)
  insertTab(session=session,  inputId = "menu", target = "DEGs",  position = "after",
            tabPanel("Gene Set Enrichment", value = "gsea", geneset_ui(id = "GS")) )
  geneset_server(id = "GS")
  source("wgcna.R",local = TRUE)
  insertTab(session=session,  inputId = "menu", target = "gsea",  position = "after",
            tabPanel("WGCNA", value = "wgcna", wgcna_ui(id = "wgcna")) )
  wgcna_server(id = "wgcna", parent_session = session)
  source("correlation.R",local = TRUE)
  insertTab(session=session,  inputId = "menu", target = "wgcna",  position = "after",
            tabPanel("Correlation Analysis", value = "Correlation", correlation_ui(id = "Corr")) )
  correlation_server(id = "Corr", parent_session = session)
  source("pattern.R",local = TRUE) 
  source("TimeSeries.R",local = TRUE)
  insertTab(session=session,  inputId = "menu", target = "Pattern_Clustering",  position = "after",
            tabPanel("Time Course Analysis", value = "time_series", TimeSeries_ui(id = "time_series")) )
  TimeSeries_server(id = "time_series", parent_session = session)
  source("vennprojects.R",local = TRUE)
  source("network.R",local = TRUE)
  source("help.R",local = TRUE)
  source("output.R",local = TRUE)
  source("scurve.R",local = TRUE)  
}

shinyApp(ui, server, enableBookmarking = "server")