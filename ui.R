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
tagList(tags$head(tags$style(type = 'text/css','.navbar-brand{display:none;}')),
fluidPage(theme = shinytheme("cerulean"),
          windowTitle = "Quickomics",
          titlePanel(
            fluidRow(
              column(4, img(height =75 , src = "Quickomics.png")), 
              column(8,  h2(strong(textOutput('project')), align = 'left'))
            ),
              windowTitle = "Quickomics" ),
            
          
navbarPage(title ="", 
	
##########################################################################################################
## Select Dataset
##########################################################################################################
tabPanel("Select Dataset",
         fluidRow(
           column(3,
                  wellPanel(
                    radioButtons("select_dataset",label="Select data set", choices=c("Saved Projects","Upload RData File"),inline = F, selected="Saved Projects"),
                    conditionalPanel("input.select_dataset=='Saved Projects'",
                                     selectInput("sel_project", label="Available Dataset", 
                                                 choices=c("", projects), selected=NULL)),
                    conditionalPanel("input.select_dataset=='Upload RData File'",
                                     fileInput("file1", "Choose data file"),
                                     fileInput("file2", "Choose network file"),                                     
                                     uiOutput('ui.action') )
                  )
           ),
           column(9,
                  tabsetPanel(id="Tables",
                              #tabPanel(title="Introduction",htmlOutput('intro')),
                              #tabPanel(title="Project Table", DT::dataTableOutput('projecttable')),
                              tabPanel(title="Sample Table", actionButton("sample", "Save to output"), dataTableOutput('sample')),
                              tabPanel(title="Project Overview", htmlOutput("summary"), tableOutput('group_table')),                              
                              tabPanel(title="Result Table", actionButton("results", "Save to output"), dataTableOutput('results')),
                              tabPanel(title="Data Table", actionButton("data_wide", "Save to output"), dataTableOutput('data_wide')),
                              tabPanel(title="Protein Gene Names", actionButton("ProteinGeneName", "Save to output"), dataTableOutput('ProteinGeneName')),
                              tabPanel(title="Help", htmlOutput('help_input'))
                  )
           )
         )
), 

##########################################################################################################
## QC Plots
##########################################################################################################
tabPanel("QC Plots",
	fluidRow(
		column(3,
			wellPanel(
			  tags$style(mycss),
				selectizeInput("QC_groups", label="Select Groups", choices=NULL, multiple=TRUE),
				selectizeInput("QC_samples", label="Select Samples", choices=NULL,multiple=TRUE),
				selectInput("PCAcolorby", label="Color By", choices=NULL),
				conditionalPanel("input.groupplot_tabset=='PCA Plot' || input.groupplot_tabset=='PCA 3D Interactive'",
				                 selectInput("PCAshapeby", label="Shape By", choices=NULL)),
				conditionalPanel("input.groupplot_tabset=='PCA Plot'",
				  selectInput("PCAsizeby", label="Size By", choices=NULL),
					selectizeInput("pcnum",label="Select Principal Components", choices=1:10, multiple=TRUE, selected=1:2, options = list(maxItems = 2)),
					radioButtons("ellipsoid", label="Plot Ellipsoid (>3 per Group)", inline = TRUE, choices = c("No" = FALSE,"Yes" = TRUE)),
					radioButtons("mean_point", label="Show Mean Point", inline = TRUE, choices = c("No" = FALSE,"Yes" = TRUE)),
					radioButtons("rug", label="Show Marginal Rugs", inline = TRUE, choices = c("No" = FALSE,"Yes" = TRUE)),
					selectInput("PCAcolpalette", label= "Select palette", choices=c("Accent"="Accent","Dark2"="Dark2","Paired"="Paired","Pastel1"="Pastel1","Pastel2"="Pastel2","Set1"="Set1","Set2"="Set2","Set3"="Set3"), selected="Dark2"),
					sliderInput("PCAdotsize", "Dot Size (when size by not used):", min = 1, max = 20, step = 1, value = 4),
					radioButtons("PCA_subsample", label="Label Samples:", inline = TRUE, choices = c("All","None", "Subset"), selected = "All"),
					conditionalPanel("input.PCA_subsample!='None'",
					               sliderInput("PCAfontsize", "Label Font Size:", min = 1, max = 20, step = 1, value = 10),
					               radioButtons("PCA_label",label="Select Sample Label",inline = TRUE, choices="")),
					conditionalPanel("input.PCA_subsample=='Subset'",
					                 actionButton("PCA_refresh_sample", "Reload Sample IDs"),
					                 textAreaInput("PCA_list", "List of Samples to Label", "", cols = 5, rows=6))
				),
				conditionalPanel("input.groupplot_tabset=='PCA 3D Plot'",
				                 radioButtons("ellipsoid3d", label="Plot Ellipsoid (>3 per Group)", inline = TRUE, choices = c("No" = "No","Yes" = "Yes")),
				                 radioButtons("dotlabel", label="Dot Label", inline = TRUE, choices =  c("No" = "No","Yes" = "Yes"))
				                 
				),
				conditionalPanel("input.groupplot_tabset=='Dendrograms'",
					sliderInput("DendroCut", label="tree cut number:", min = 2, max = 10, step = 1, value = 4),
					sliderInput("DendroFont", label= "Label Font Size:", min = 0.5, max = 4, step = 0.5, value = 1),
					radioButtons("dendroformat", label="Select Plot Format", inline = TRUE, choices = c("tree" = "tree","horizontal" = "horiz", "circular" = "circular"), selected="circular")
				)
				
			)
		),
		column(9,
			tabsetPanel(id="groupplot_tabset",
				tabPanel(title="PCA Plot", actionButton("pcaplot", "Save to output"), plotOutput("pcaplot",height = 800)),
				tabPanel(title="Eigenvalues",  plotOutput("Eigenvalues",height = 650)),
				tabPanel(title="PCA 3D Plot",  plotOutput("pca_legend",height = 100), rglwidgetOutput("plot3d",  width = 1000, height = 1000)),
				tabPanel(title="PCA 3D Interactive", plotlyOutput("plotly3d",  width = 1000, height = 1000)),
				tabPanel(title="Sample-sample Distance",actionButton("SampleDistance", "Save to output"), plotOutput("pheatmap",height = 800)),
				tabPanel(title="Dendrograms",actionButton("Dendrograms", "Save to output"), plotOutput("Dendrograms",height = 800)),
				tabPanel(title="Box Plot", actionButton("QCboxplot", "Save to output"), plotOutput("QCboxplot",height = 800)),
				tabPanel(title="CV Distribution", actionButton("histplot", "Save to output"), plotOutput("histplot",height = 800)),
				tabPanel(title="Order Groups", actionButton("reset_group", "Reset Groups"), uiOutput('reorder_group')),
				tabPanel(title="Help", htmlOutput('help_QC'))
			)
		)
	)
), 

##########################################################################################################
## Volcano Plot
##########################################################################################################
tabPanel("Volcano Plot",
	fluidRow(
		column(3,
			wellPanel(
			  conditionalPanel( "input.valcano_tabset!='DEGs in Two Comparisons'",
				selectInput("valcano_test", label="Select Comparison Groups for Volcano Plot", choices=NULL)),
				conditionalPanel( "input.valcano_tabset=='DEGs in Two Comparisons'",
				selectInput("valcano_test1", label="1st Comparison (X-axis)", choices=NULL),
				selectInput("valcano_test2", label="2nd Comparison (Y-aixs)", choices=NULL)),
				numericInput("valcano_FCcut", label= "Choose Fold Change Cutoff",  value = 1.2, min=1, step=0.1),
				numericInput("valcano_pvalcut", label= "Choose P Value Cutoff", value=0.01, min=0, step=0.001),
  			radioButtons("valcano_psel", label= "P value or P.adj Value?", choices= c("Pval"="Pval","Padj"="Padj"),inline = TRUE),
				conditionalPanel( "input.valcano_tabset!='DEGs in Two Comparisons'",		
				  textOutput("valcano_filteredgene"),
				  tags$head(tags$style("#valcano_filteredgene{color: red; font-size: 20px; font-style: italic; }" ))),
				conditionalPanel( "input.valcano_tabset!='Volcano Plot (Interactive)'",
					radioButtons("valcano_label",label="Select Gene Label",inline = TRUE, choices=""),
					radioButtons("volcano_label", label="Label Genes:", inline = TRUE, choices = c("DEGs","None", "Upload", "Geneset"), selected = "DEGs"),
					conditionalPanel("input.volcano_label!='None'",
					    sliderInput("Ngenes", "# of Genes to Label", min = 10, max = 200, step = 5, value = 50)),
					conditionalPanel("input.volcano_label=='Upload'",
					    textAreaInput("volcano_gene_list", "List of genes to label\n(UniqueID, Gene.Name or Protein.ID)", "", cols = 5, rows=6)
						),
          conditionalPanel("input.volcano_label=='Geneset'",
              uiOutput("html_geneset") )	
					),
				radioButtons("more_options", label="Show More Options", inline = TRUE, choices = c("Yes","No"), selected = "No"),
				conditionalPanel("input.more_options=='Yes'",
				  numericInput("Max_logFC", label= "Max abs(logFC) in plot (use 0 for full range)", value=0, min=0),
				  conditionalPanel( "input.valcano_tabset!='DEGs in Two Comparisons'",
				                    numericInput("Max_Pvalue", label= "Max -log10(Stat Value) in plot (use 0 for full range)", value=0, min=0) ),
				  conditionalPanel( "input.valcano_tabset!='Volcano Plot (Interactive)'",
				                    sliderInput("lfontsize", "Label Font Size:", min = 1, max = 10, step = 1, value = 4),
				                    sliderInput("yfontsize", "Legend Font Size:", min = 8, max = 24, step = 1, value = 14)),
				  conditionalPanel( "input.valcano_tabset=='DEGs in Two Comparisons'",
				                    radioButtons("DEG_comp_XY", label="Make X and Y scale the same?", inline = TRUE, choices = c("Yes","No"), selected = "Yes"))
				                    
				)
			)
		),
		column(9,
			tabsetPanel(id="valcano_tabset",
				tabPanel(title="Volcano Plot (Static)",actionButton("valcano", "Save to output"), plotOutput("volcanoplotstatic", height=800)),
				tabPanel(title="Volcano Plot (Interactive)", plotlyOutput("volcanoplot", height=800)),
			  tabPanel(title="DEGs in Two Comparisons",actionButton("DEG_comp", "Save to output"), plotOutput("DEG_Compare", height=800)),
				tabPanel(title="Data Table", actionButton("DEG_data", "Save to output"), DT::dataTableOutput("valcanoData")),
				tabPanel(title="Help", htmlOutput('help_volcano'))
			)
		)
	)
), 

##########################################################################################################
## Heat Map
##########################################################################################################
tabPanel("Heat Map",
	fluidRow(
		column(3,
			wellPanel(
				#actionButton("action_heatmaps","Generate Interactive Heatmap"),
				selectizeInput("heatmap_groups", label="Select Groups", choices=NULL, multiple=TRUE),
				selectizeInput("heatmap_samples", label="Select Samples", choices=NULL, multiple=TRUE),
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
				  numericInput("maxgenes",label="Choose Gene Number", min=1, max= 5000, value=100, step=1)), 
				conditionalPanel("input.heatmap_subset=='Subset'",
					selectInput("heatmap_test", label="Select Test to Select genes", choices=NULL),
					column(width=6,numericInput("heatmap_fccut", label= "Fold Change Cutoff", value = 1.2, min=1, step=0.1)),
					column(width=6,numericInput("heatmap_pvalcut", label= "P Value Cutoff", value=0.01, min=0, step=0.001)),
					radioButtons("heatmap_psel", label= "P value or P.adj Value?", choices= c("Pval"="Pval","Padj"="Padj"),inline = TRUE),
					column(width=12,textOutput("heatmapfilteredgene"),
					tags$head(tags$style("#heatmapfilteredgene{color: red; font-size: 20px; font-style: italic; }")))
					), 
				conditionalPanel( "input.heatmap_tabset=='Static Heatmap Layout 1'",
				    selectizeInput("heatmap_annot", label="Annotate Samples", choices=NULL, multiple = TRUE)),			                  
				    column(width=5,selectInput("dendrogram", "Apply Clustering:", c("both" ,"none", "row", "column"))),
				    column(width=5,selectInput("scale", "Apply Scaling:", c("none","row", "column"),selected="row")),
				conditionalPanel( "input.heatmap_tabset=='Static Heatmap Layout 2'",
					column(width=5,selectInput("key", "Color Key:", c("TRUE", "FALSE"))),
					column(width=5,selectInput("srtCol", "angle of label", c("45", "60","90"))),
					column(width=5,sliderInput("hxfontsize", "Column Font Size:", min = 0, max = 3, step = 0.5, value = 1)),
					column(width=5,sliderInput("hyfontsize", "Row Font Size:", min = 0, max = 3, step = 0.5, value = 1)),
					column(width=5,sliderInput("right", "Set Margin Width", min = 0, max = 20, value = 5)),
					column(width=5,sliderInput("bottom", "Set Margin Height", min = 0, max = 20, value = 5))
				),
				conditionalPanel( "input.heatmap_tabset=='Interactive Heatmap'",
					column(width=5,selectInput("key", "Color Key:", c("TRUE", "FALSE"))),
					column(width=5,selectInput("srtCol", "angle of label", c("45", "60","90"))),
					column(width=5,sliderInput("hxfontsizei", "Column Font Size:", min = 0, max = 3, step = 0.5, value = 1)),
					column(width=5,sliderInput("hyfontsizei", "Row Font Size:", min = 0, max = 3, step = 0.5, value = 1)),
					column(width=5,sliderInput("l", "Set Margin Width", min = 0, max = 200, value = 120)),
					column(width=5,	sliderInput("b", "Set Margin Height", min = 0, max = 200, value = 120))
				),
				conditionalPanel( "input.heatmap_tabset=='Static Heatmap Layout 1'",
				  column(width=5,sliderInput("hxfontsizep", "Column Font Size:", min = 0, max = 20, step = 1, value = 10)),
					column(width=5,sliderInput("hyfontsizep", "Row Font Size:", min = 0, max = 20, step = 1, value = 7)),
				  radioButtons("heatmap_label",label="Gene Label",inline = TRUE, choices=""),
				  sliderInput("heatmap_N_genes", "Max Number of Genes to Label:", min = 0, max = 500, step = 10, value = 100),
  				h5("After changing parameters, please click Plot/Refresh button in the plot panel to generate heatmap."),
					radioButtons("heatmap_highlight", label="Highlight Subset of Genes:", inline = TRUE, choices = c("Yes","No"), selected = "No"),
					conditionalPanel("input.heatmap_highlight=='Yes'",
					               uiOutput("gene_highlight_file"),
					              sliderInput("hl_font_size", "Font Size:", min = 0, max = 20, step = 1, value = 9)),
	  			radioButtons("heatmap_more_options", label="Show More Options", inline = TRUE, choices = c("Yes","No"), selected = "No"),
		  		conditionalPanel("input.heatmap_more_options=='Yes'",
				                 sliderInput("heatmap_height", "Heatmap Height:", min = 200, max = 3000, step = 50, value = 800),
				                 radioButtons("heatmap_row_dend", label="Show row dendrogram", inline = TRUE, choices = c("No" = FALSE,"Yes" = TRUE), selected = TRUE),
				                 radioButtons("heatmap_col_dend", label="Show column dendrogram", inline = TRUE, choices = c("No" = FALSE,"Yes" = TRUE), selected = TRUE),
				                 column(width=3,colourInput("lowColor", "Low", "blue")),
				                 column(width=3,colourInput("midColor", "Mid", "white")),
				                 column(width=3,colourInput("highColor", "High", "red")),
				                 column(width=5,selectInput("distanceMethod", "Distance Metric:", c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"))),
				                 column(width=5,selectInput("agglomerationMethod", "Linkage Algorithm:", c("complete", "single", "average", "centroid", "median", "mcquitty", "ward.D", "ward.D2"))),
				                 column(width=5,sliderInput("cutreerows", "cutree_rows:", min = 0, max = 8, step = 1, value = 0)),
				                 column(width=5,sliderInput("cutreecols", "cutree_cols:", min = 0, max = 8, step = 1, value = 0)),
				                 radioButtons("custom_color", label="Upload Colors for Annotations", inline = TRUE, choices = c("Yes","No"), selected = "No"),
				                 conditionalPanel("input.custom_color=='Yes'",
				                                  uiOutput("annot_color_file"))
				                                  
				                 
				)
				)			
			)
		),
		column(9,
			tabsetPanel(id="heatmap_tabset",
				tabPanel(title="Static Heatmap Layout 1",actionButton("pheatmap2", "Save to output"),  
				         actionButton("plot_heatmap", "Plot/Refresh", style="color: #0961E3; background-color: #F6E98C ; border-color: #2e6da4"), 
				         uiOutput("plot.heatmap")),
				tabPanel(title="Static Heatmap Layout 2", actionButton("staticheatmap", "Save to output"), plotOutput("staticheatmap", height = 800)),
				tabPanel(title="Interactive Heatmap",textOutput("text"), p(), plotlyOutput("interactiveheatmap", height = 800)),
				tabPanel(title="Help", htmlOutput('help_heatmap'))
			)
		)
	)
),


##########################################################################################################
## Protein Expression Plot
##########################################################################################################
tabPanel("Expression Plot",
         fluidRow(
           column(3,
                  wellPanel(
                    conditionalPanel("input.expression_tabset=='Searched Expression Data' || input.expression_tabset=='Data Output' ",
                                    radioButtons("exp_subset",label="Genes used in plot", choices=c("Select", "Upload Genes", "Geneset"),inline = TRUE, selected="Select"),
                                    conditionalPanel("input.exp_subset=='Upload Genes'", 
                                        textAreaInput("exp_list", "Enter Gene List", "", cols = 5, rows=6)),
                                    conditionalPanel("input.exp_subset=='Select'", 
                                        selectizeInput("sel_gene",	label="Gene Name (Select 1 or more)",	choices = NULL,	multiple=TRUE, options = list(placeholder =	'Type to search'))),
                                    conditionalPanel("input.exp_subset=='Geneset'",   uiOutput("html_geneset_exp") ),  
                                    radioButtons("SeparateOnePlot", label="Separate or One Plot", inline = TRUE, choices = c("Separate" = "Separate", "OnePlot" = "OnePlot"))),
                    conditionalPanel("input.expression_tabset=='Browsing'",
                                     column(width=6,numericInput("expression_fccut", label= "Choose Fold Change Threshold",  value = 1.2, min=1, step=0.1)),
                                     column(width=6,numericInput("expression_pvalcut", label= "Choose P-value Threshold",  value=0.01, min=0, step=0.001)),
                                     radioButtons("expression_psel", label= "P value or P.adj Value?", choices= c("Pval"="Pval","Padj"="Padj"),inline = TRUE),
                                     selectInput("expression_test", label="Select Test", choices=NULL),
                                     textOutput("expfilteredgene"),
                                     tags$head(tags$style("#expfilteredgene{color: red; font-size: 20px; font-style: italic;}")),
                                     column(width=6,selectInput("sel_page",	label="Select Page",	choices = NULL,	selected=1)),
                                     column(width=6,selectInput("numperpage", label= "Plot Number per Page", choices= c("4"=4,"6"=6,"9"=9), selected=6))
                    ),
                    selectizeInput("sel_group", label="Select Groups (u can re-order)", choices=NULL, multiple=TRUE),
                    radioButtons("sel_geneid",label="Select Gene Label",inline = TRUE, choices=""),
                    radioButtons("plotformat", label="Select Plot Format", inline = TRUE, choices = c("Box Plot" = "boxplot","Bar Plot" = "barplot", "violin" = "violin","line" = "line")),
                    radioButtons("IndividualPoint", label="Show Individual Point?", inline = TRUE, choices = c("YES" = "YES","NO" = "NO")),
                    radioButtons("ColPattern", label="Bar Colors", inline = TRUE, choices = c("Palette" = "Palette", "Single" = "Single")),
                    conditionalPanel("input.ColPattern=='Single'",
                                     colourInput("barcol", "Select colour", "#1E90FF", palette = "limited")),
                    conditionalPanel("input.ColPattern=='Palette'",
                                     selectInput("colpalette", label= "Select palette", choices=c("Accent"="Accent","Dark2"="Dark2","Paired"="Paired","Pastel1"="Pastel1","Pastel2"="Pastel2","Set1"="Set1","Set2"="Set2","Set3"="Set3"), selected="Dark2")),
                    sliderInput("expression_axisfontsize", "Axis Font Size:", min = 10, max = 28, step = 1, value = 16),
                    sliderInput("expression_titlefontsize", "Title Font Size:", min = 12, max = 28, step = 1, value = 16),
                    textInput("Ylab", "Y label", width = "100%"),
                    textInput("Xlab", "X label", width = "100%"),
                    sliderInput("Xangle", label= "X Angle", min = 0, max = 90, step = 15, value = 45),
                    radioButtons("exp_plot_Y_scale", label="Y Axis Scale", inline = TRUE, choices = c("Auto","Manual"), selected = "Auto"),
                    conditionalPanel("input.exp_plot_Y_scale=='Manual'",
                      column(width=6,numericInput("exp_plot_Ymin", label= "Y Min",  value = 0, step=0.1)),
                      column(width=6,numericInput("exp_plot_Ymax", label= "Y Max",  value=5, step=0.1))),
                    conditionalPanel("input.expression_tabset=='Searched Expression Data'",
                      h5("After changing parameters, please click Plot/Refresh button in the plot panel to generate expression plot."))
                    #sliderInput("hjust", label= "hjust", min = 0, max = 1, step = 0.1, value = 1),
                    #sliderInput("vjust", label= "vjust", min = 0, max = 1, step = 0.1, value = 0.5)
                  )
           ),
           column(9,
                  tabsetPanel(id="expression_tabset",
                              tabPanel(title="Browsing",actionButton("browsing", "Save to output"),plotOutput("browsing", height=800)),
                              tabPanel(title="Searched Expression Data",actionButton("boxplot", "Save to output"),
                                       actionButton("plot_exp", "Plot/Refresh", style="color: #0961E3; background-color: #F6E98C ; border-color: #2e6da4"), 
                                       uiOutput("plot.exp")),
                              tabPanel(title="Data Table",	DT::dataTableOutput("dat_dotplot")),
                              tabPanel(title="Result Table",	DT::dataTableOutput("res_dotplot")),
                              tabPanel(title="Help", htmlOutput('help_expression'))
                  )
           )
         )
),

##########################################################################################################
## Gene Set Enrichment
##########################################################################################################
tabPanel("Gene Set Enrichment",
	fluidRow(
		column(3,
			wellPanel(
				selectInput("geneset_test", label="Select Comparison for Gene Set analysis", choices=NULL),
				conditionalPanel("input.geneset_tabset!='KEGG Pathway View'",
				  column(width=6,numericInput("geneset_FCcut", label= "Fold Change Cutoff", value = 1.2, min=1, step=0.1)),
				  column(width=6,numericInput("geneset_pvalcut", label= "P Value Cutoff", value=0.01, min=0, step=0.001)),
				  radioButtons("geneset_psel", label= "P value or P.adj Value?", choices= c("Pval"="Pval","Padj"="Padj"),inline = TRUE),
				  radioButtons("geneset_direction", label= "Up- or Down-Regulated Genes?", choices= c("Both", "Up", "Down"),inline = TRUE),
				  textOutput("geneset_filteredgene"),
				  tags$head(tags$style("#geneset_filteredgene{color: red; font-size: 20px; font-style: italic;}")),
				),
				conditionalPanel("input.geneset_tabset=='KEGG Pathway View'",
				  radioButtons("kegg_more_tests", label= "Add more comparisons?", choices= c("Yes", "No"),selected="No", inline = TRUE),
				  conditionalPanel("input.kegg_more_tests=='Yes'",
				  selectInput("geneset_test2", label="2nd Comparison", choices=NULL),
				  selectInput("geneset_test3", label="3rd Comparison", choices=NULL),
				  selectInput("geneset_test4", label="4th Comparison", choices=NULL),
				  selectInput("geneset_test5", label="5th Comparison", choices=NULL)),
				),				
				conditionalPanel("input.geneset_tabset=='Gene Set Heat Map'",
				radioButtons("gs_heatmap_label",label="Gene Label",inline = TRUE, choices=c("UniqueID", "Gene.Name"), selected="Gene.Name")),
				conditionalPanel("input.geneset_tabset=='Gene Set Enrichment'",
					radioButtons("MSigDB", label= "MSigDB Collections",
						choices= c("KEGG Pathway" = "KEGG",
							"c2.cgp (chemical and genetic perturbations)" = "c2.cgp.v6.1",
							"c2.cp.biocarta" = "c2.cp.biocarta.v6.1",
							"c2.cp.reactome" =  "c2.cp.reactome.v6.1",
							"c2.cp (Canonical pathways)" = "c2.cp.v6.1",
							"c3.all (motif gene sets)" = "c3.all.v6.1",
							"c3.tft (transcription factor targets)" = "c3.tft.v6.1",
							"c5.bp (GO biological process)" = "c5.bp.v6.1",
							"c5.cc (GO cellular component)" = "c5.cc.v6.1",
							"c5.mf (GO molecular function)" = "c5.mf.v6.1",
							"c6.all (oncogenic signatures)" = "c6.all.v6.1",
							"c7.all (immunologic signatures)" = "c7.all.v6.1",
						"h.all.v6.1 (hallmark gene sets)" = "h.all.v6.1"), selected = "KEGG")
					),
				conditionalPanel("input.geneset_tabset=='KEGG Pathway View'",
				      selectInput("kegg_logFC", label= "KEGG view logFC Range:", choices= c(1, 2, 3), selected=1),
				      radioButtons("kegg_mapsample", label= "Map Symbols to KEGG Nodes?", choices= c("YEs"=TRUE, "No"=FALSE),inline = TRUE)),
				)
			),
			column(9,
				tabsetPanel(id="geneset_tabset",
					tabPanel(title="Gene Set Enrichment", DT::dataTableOutput("MSigDB")),
					tabPanel(title="Gene Expression",textInput('x1', 'Row ID'),  DT::dataTableOutput("Expression")),
					tabPanel(title="Gene Set Heat Map", textInput('x2', 'Row ID'), plotOutput('SetHeatMap',height="auto", width = 900)),
					tabPanel(title="KEGG Pathway View", textInput('x3', 'Row ID'), plotOutput('keggView')),
					tabPanel(title="Help", htmlOutput('help_geneset'))
			)
			)
		)
	), 
	
##########################################################################################################
## Pattern Clustering
##########################################################################################################
## 10/08/2020
## eidted by bgao, add gene upload box and more plotting options	   
tabPanel("Pattern Clustering",
	fluidRow(
		column(3,
			wellPanel(

			  radioButtons("pattern_subset",label="Use subset genes or upload your own subset?", choices=c("subset","upload genes"),inline = TRUE, selected="subset"),
				conditionalPanel("input.pattern_subset=='subset'",
				                 column(width=6,numericInput("pattern_fccut", label= "Choose Fold Change Threshold",value = 1.2, min=1, step=0.1)),
				                 column(width=6,numericInput("pattern_pvalcut", label= "Choose P-value Threshold", value=0.01, min=0, step=0.001)),
				                 radioButtons("pattern_psel", label= "P value or P.adj Value?", choices= c("Pval"="Pval","Padj"="Padj"),inline = TRUE),
				                 textOutput("patternfilteredgene")
				),
				conditionalPanel("input.pattern_subset=='upload genes'", textAreaInput("pattern_list", "list", "", cols = 5, rows=6)),
				
				tags$head(tags$style("#patternfilteredgene{color: red; font-size: 20px; font-style: italic; }")),
				selectizeInput("pattern_group", label="Select Groups (re-order under QC tab)", choices=NULL, multiple=TRUE),
				radioButtons("ClusterMehtod", label="Cluster Method", inline = FALSE, choices = c("Soft Clustering" = "mfuzz", "K-means" = "kmeans", "Partitioning Around Medoids (disabled)" = "pam")),
				sliderInput("k", "Cluster Number:", min = 3, max = 12, step = 1, value = 6),
				conditionalPanel("input.ClusterMehtod=='kmeans'",
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
				tabPanel(title="Clustering of Centroid Profiles",actionButton("pattern", "Save to output"),plotOutput("pattern", height=800)),
				#tabPanel(title="Data Table",actionButton("Pattern_data", "Save to output"), DT::dataTableOutput("dat_pattern")),
				tabPanel(title="Data Table", DT::dataTableOutput("dat_pattern")),
				tabPanel(title="Help", htmlOutput('help_pattern'))
			)
		)
	)
),

##########################################################################################################
## Correlation Network
##########################################################################################################
tabPanel("Correlation Network",
	fluidRow(
		column(3,
			wellPanel(
			  radioButtons("network_label",label="Select Gene Label",inline = TRUE, choices=c("UniqueID", "Gene.Name"), selected="Gene.Name"),
				selectizeInput("sel_net_gene",	label="Gene Name (Select 1 or more)",	choices = NULL,	multiple=TRUE, options = list(placeholder =	'Type to search')),
				sliderInput("network_rcut", label= "Choose r Cutoff",  min = 0.7, max = 1, value = 0.9, step=0.02),
	 			selectInput("network_pcut", label= "Choose P Value Cutoff", choices= c("0.0001"=0.0001,"0.001"=0.001,"0.01"=0.01,"0.05"=0.05),selected=0.01),
       	textOutput("networkstat"),
        uiOutput("myTabUI")
			)
		),
		column(9,
			tabsetPanel(id="Network_tabset",
			tabPanel("visNetwork", visNetworkOutput("visnetwork", height="800px"), style = "background-color: #eeeeee;"),
			tabPanel("networkD3(disabled)"), #forceNetworkOutput("networkD3", height="800px"), style = "background-color: #eeeeee;"),
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
		column(2,
			wellPanel(
				column(width=6,numericInput("venn_fccut", label= "Choose Fold Change Cutoff", value = 1.2, min=1, step=0.1)),
				column(width=6,numericInput("venn_pvalcut", label= "Choose P value Cutoff", value=0.01, min=0, step=0.001)),
				radioButtons("venn_psel", label= "P value or P.adj Value?", choices= c("Pval"="Pval","Padj"="Padj"),inline = TRUE),
				radioButtons("venn_updown", label= "All, Up or Down?", choices= c("All"="All","Up"="Up","Down"="Down"),inline = TRUE),
				selectInput("venn_test1", label="Select List 1", choices=NULL),
				conditionalPanel("input.venn_tabset=='vennDiagram'",
				colourInput("col1", "Select colour", "#0000FF",palette = "limited")
				),

				selectInput("venn_test2", label="Select List 2", choices=NULL),
				conditionalPanel("input.venn_tabset=='vennDiagram'",
				colourInput("col2", "Select colour", "#FF7F00",palette = "limited")
				),

				selectInput("venn_test3", label="Select List 3", choices=NULL),
				conditionalPanel("input.venn_tabset=='vennDiagram'",
				colourInput("col3", "Select colour", "#00FF00",palette = "limited")
				),

				selectInput("venn_test4", label="Select List 4", choices=NULL),
				conditionalPanel("input.venn_tabset=='vennDiagram'",
				colourInput("col4", "Select colour", "#FF00FF",palette = "limited")
				),

				selectInput("venn_test5", label="Select List 5", choices=NULL),
				conditionalPanel("input.venn_tabset=='vennDiagram'",
				colourInput("col5", "Select colour", "#FFFF00", palette = "limited")
				),
				conditionalPanel("input.venn_tabset=='Intersection Output'",
				radioButtons("vennlistname", label= "Label name", choices= c("Gene"="Gene","AC Number"="AC", "UniqueID"="UniqueID"),inline = TRUE, selected = "Gene"))
				) 
			),
			column(10,
				tabsetPanel(id="venn_tabset",
					tabPanel(title="vennDiagram",
						column(9,
							actionButton("vennDiagram", "Save to output"),plotOutput("vennDiagram", height = 800)
						),
						column(3,
							textInput("title", "Title", width = "100%"),
							sliderInput("maincex", "Title Size", min = 0, max = 6, value = 3, width = "100%"),
							sliderInput("alpha", "Opacity", min = 0, max = 1, value = 0.4, width = "100%"),
							sliderInput("lwd", "line thick", min = 1, max = 4, value = 1, width = "100%"),
							sliderInput("lty", "line type", min = 1, max = 6, value = 1, width = "100%"),
							radioButtons("fontface","Number Font face",list("plain", "bold", "italic"),selected = "plain",	inline = TRUE),
							sliderInput("cex", "Font size", min = 1, max = 4, value = 2, width = "100%"),
							radioButtons("catfontface","Label Font face",list("plain", "bold", "italic"),	selected = "plain", inline = TRUE),
							sliderInput("catcex", "Font size", min = 1, max = 2, step=0.1, value = 1.8, width = "100%"),
							sliderInput("margin", "Margin", min = 0, max = 1, step=0.05, value = 0.1, width = "100%")
						)
					),
					tabPanel(title="VennDiagram(black & white)", plotOutput("SvennDiagram",height = 800, width = 800)),
					tabPanel(title="Intersection Output", htmlOutput("vennHTML")),
					tabPanel(title="DEG Table", actionButton("venn_DEG_data", "Save to output"), DT::dataTableOutput("venn_DEG_Data")),
					tabPanel(title="Help", htmlOutput('help_venn'))
				)
			)
		)
),


##########################################################################################################
## Venn Across Projects
##########################################################################################################
tabPanel("Venn Across Projects",
	fluidRow(
		column(3,
			wellPanel(
			  numericInput("vennP_fccut", label= "Choose Fold Change Cutoff", value = 1.2, min=1, step=0.1),
				radioButtons("vennP_psel", label= "P value or P.adj Value?", choices= c("Pval"="Pval","Padj"="Padj"),inline = TRUE),
				numericInput("vennP_pvalcut", label= "Choose P value Cutoff", value=0.01, min=0, step=0.001),

				selectInput("dataset1", "Data set1", choices=NULL),
				selectInput("vennP_test1", label="Select List 1", choices=NULL),
				conditionalPanel("input.vennP_tabset=='vennDiagram'",
				colourInput("vennPcol1", "Select colour", "#0000FF",palette = "limited")
				),
				selectInput("dataset2", "Data set2", choices=NULL),
				selectInput("vennP_test2", label="Select List 2", choices=NULL),
				conditionalPanel("input.vennP_tabset=='vennDiagram'",
				colourInput("vennPcol2", "Select colour", "#FF7F00",palette = "limited")
				),
				selectInput("dataset3", "Data set3", choices=NULL),
				selectInput("vennP_test3", label="Select List 3", choices=NULL),
				conditionalPanel("input.vennP_tabset=='vennDiagram'",
				colourInput("vennPcol3", "Select colour", "#00FF00",palette = "limited")
				),
				selectInput("dataset4", "Data set4", choices=NULL),
				selectInput("vennP_test4", label="Select List 4", choices=NULL),
				conditionalPanel("input.vennP_tabset=='vennDiagram'",
					colourInput("vennPcol4", "Select colour", "#FF00FF",palette = "limited")
				),
				selectInput("dataset5", "Data set5", choices=NULL),
				selectInput("vennP_test5", label="Select List 5", choices=NULL),
				conditionalPanel("input.vennP_tabset=='vennDiagram'",
				colourInput("vennPcol5", "Select colour", "#FFFF00", palette = "limited")
				)
			) 
		), 
		column(9,
			tabsetPanel(id="vennP_tabset",
				tabPanel(title="vennDiagram",
					column(9,
						plotOutput("vennPDiagram", height = 800)
					),
					column(3,
						textInput("vennPtitle", "Title", width = "100%"),
						sliderInput("vennPmaincex", "Title Size", min = 0, max = 6, value = 3, width = "100%"),
						sliderInput("vennPalpha", "Opacity", min = 0, max = 1, value = 0.4, width = "100%"),
						sliderInput("vennPlwd", "line thick", min = 1, max = 4, value = 1, width = "100%"),
						sliderInput("vennPlty", "line type", min = 1, max = 6, value = 1, width = "100%"),
						radioButtons("vennPfontface",	"Number Font face", list("plain", "bold", "italic"),	selected = "plain",	inline = TRUE),
						sliderInput("vennPcex", "Font size", min = 1, max = 4, value = 2, width = "100%"),
						radioButtons("vennPcatfontface", "Label Font face", list("plain", "bold", "italic"),	selected = "plain",	inline = TRUE),
						sliderInput("vennPcatcex", "Font size", min = 1, max = 2, step=0.1, value = 1.8, width = "100%"),
						sliderInput("vennPmargin", "Margin", min = 0, max = 1, step=0.05, value = 0.2, width = "100%")
					)
				),
				tabPanel(title="VennDiagram(black & white)", plotOutput("SvennPDiagram",height = 800,width = 800)),
				tabPanel(title="Intersection Output", htmlOutput("vennPHTML")),
				tabPanel(title="Help", htmlOutput('help_vennp'))
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
	downloadButton('downloadXLSX', 'Download tables in .xlsx')
	),


##########################################################################################################
## footer
##########################################################################################################

footer= HTML(footer_text)


)
)
) #for tagList
