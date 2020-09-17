###########################################################################################################
## Proteomics Visualization R Shiny App
##
##This software belongs to Biogen Inc. All right reserved.
##
##@file: input.R
##@Developer : Benbo Gao (benbo.gao@Biogen.com)
##@Date : 5/31/2019
##@version 1.0
###########################################################################################################
saved_plots <- reactiveValues()  
saved_table <- reactiveValues() 
group_order <- reactiveVal()
#saved_palette <- reactiveVal()

ProjectInfo<-reactive({
  req(input$sel_project)
  ProjectID=input$sel_project
  Name=saved_projects$Name[saved_projects$ProjectID==ProjectID]
  Species=saved_projects$Species[saved_projects$ProjectID==ProjectID]
  return(list(ProjectID=ProjectID, Name=Name, Species=Species))
  }) #later on can use customer uploaded data

output$project <- renderText({
  if (input$sel_project==""){"Please select or upload a date set"} else {ProjectInfo()$Name}
  })


output$ui.action <- renderUI({
  if (is.null(input$file1)) return()
  tagList(
  textInput("projet_name", label="Rename Project", value=input$file1$name),
  radioButtons("species",label="Select species", choices=c("human","mouse", "rat"), inline = F, selected="human"),
  actionButton("customData", "Submit Data")
  )
})


DataReactive <- reactive({
  withProgress(message = 'Fetching data.',
               detail = 'This may take a while...',
               value = 0,
               {
                 Pinfo=ProjectInfo()
                 #setwd('H:/Rcode/ptxvisv3')
               RDataFile <- paste("data/",  Pinfo$ProjectID, ".RData", sep = "")
               #  RDataFile <- ("D:/Test/temp/Mouse_microglia_RNA/data/Mouse_microglia_RNA-Seq.RData")   
               #  RDataFile <- ("D:/Test/temp/Mouse_microglia_RNA/data/2019_XH_OGA_iPSC_Neuron.RData")   
                 
                 load(RDataFile)
                 
                 results_long <-
                   results_long %>% mutate_if(is.factor, as.character)  %>% left_join(ProteinGeneName, ., by = "UniqueID")
                 data_long <-
                   data_long %>% mutate_if(is.factor, as.character)  %>% left_join(ProteinGeneName, ., by = "UniqueID")
                 
                 group_names <- as.character(MetaData$Order[MetaData$Order != ""])
                 if (length(group_names) == 0) {
                   group_names <- as.character(unique(MetaData$group))
                 }
                 tests  <-
                   as.character(MetaData$ComparePairs[MetaData$ComparePairs != ""])
                 tests <-  gsub("-", "vs", tests)
                 if (length(tests) == 0) {
                   tests = unique(as.character(results_long$test))
                 }
                 group_order(group_names)
                 return(
                   list(
                     "groups" = group_names,
                     "MetaData" = MetaData,
                     "results_long" = results_long,
                     "data_long" = data_long,
                     "ProteinGeneName" = ProteinGeneName,
                     "data_wide" = data_wide,
                     "data_results" = data_results,
                     "tests" = tests
                   )
                 )
               })
  
})


DataNetworkReactive <- reactive({
  DataIn = DataReactive()
  ProteinGeneName <- DataIn$ProteinGeneName
  
  query <- parseQueryString(session$clientData$url_search)

  
  Pinfo=ProjectInfo()
  CorResFile <- paste("networkdata/", Pinfo$ProjectID, ".RData", sep = "")
  if (file.exists(CorResFile)) {
    load(CorResFile)
  } else {
    data_wide <- DataIn$data_wide
    cor_res <- Hmisc::rcorr(as.matrix(t(data_wide)))
    cormat <- cor_res$r
    pmat <- cor_res$P
    ut <- upper.tri(cormat)
    network <- tibble (
      from = rownames(cormat)[row(cormat)[ut]],
      to = rownames(cormat)[col(cormat)[ut]],
      cor  = signif(cormat[ut], 2),
      p = signif(pmat[ut], 2),
      direction = as.integer(sign(cormat[ut]))
    )
    network <- network %>% mutate_if(is.factor, as.character) %>%
      dplyr::filter(!is.na(cor) & abs(cor) > 0.7 & p < 0.05)
    save(network,
         file =  paste("networkdata/", ProjectID, ".RData", sep = ""))
  }
  
  sel_gene = input$sel_net_gene
  tmpids = ProteinGeneName[unique(na.omit(c(
    apply(ProteinGeneName, 2, function(k)
      match(sel_gene, k))
  ))), ]
  
  edges.sel <-
    network %>% filter((from %in% tmpids$UniqueID) |
                         (to %in% tmpids$UniqueID))
  rcutoff <- as.numeric(input$network_rcut)
  pvalcutoff <- as.numeric(as.character(input$network_pcut))
  edges <-
    dplyr::filter(edges.sel, abs(cor) > rcutoff & p < pvalcutoff)
  networks_ids <-
    unique(c(as.character(edges$from), as.character(edges$to)))
  nodes <-
    ProteinGeneName %>% dplyr::filter(UniqueID %in% networks_ids) %>%
    dplyr::select(UniqueID, Gene.Name) %>%
    dplyr::rename(id = UniqueID, label = Gene.Name)
  net <- list("nodes" = nodes, "edges" = edges)
  return(net)
})

output$results <- DT::renderDataTable({
	DataIn <- DataReactive()
	results <- DataIn$data_results %>%
	dplyr::select(-one_of(c("Fasta.headers","UniqueID","id")))
	results[,sapply(results,is.numeric)] <- signif(results[,sapply(results,is.numeric)],3)
	DT::datatable(results, 
  options = list(
  	pageLength = 15
  ),rownames= FALSE)
})

output$sample <- DT::renderDataTable({
	DT::datatable(DataReactive()$MetaData, options = list(pageLength = 15))
	
})

output$data_wide <- DT::renderDataTable({
	DT::datatable(DataReactive()$data_wide, 
	              extensions = 'FixedColumns',
  options = list(
  	pageLength = 15,
    dom = 't',
    scrollX = TRUE,
    fixedColumns = list(leftColumns = 1)
  ))
})

output$ProteinGeneName <- DT::renderDataTable({
	DT::datatable(DataReactive()$ProteinGeneName, options = list(pageLength = 15),rownames= FALSE)
})

observeEvent(input$results, {
	DataIn <- DataReactive()
	results = DataIn$data_results 
	results[,sapply(results,is.numeric)] <- signif(results[,sapply(results,is.numeric)],3)
	saved_table$results <- results
})

observeEvent(input$sample, {
	saved_table$sample <- DataReactive()$MetaData
})

observeEvent(input$data_wide, {
	saved_table$data <- DataReactive()$data_wide
})

observeEvent(input$ProteinGeneName, {
	saved_table$ProteinGeneName <- DataReactive()$ProteinGeneName
})



