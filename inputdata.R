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
  ShortName=saved_projects$ShortNames[saved_projects$ProjectID==ProjectID]
  return(list(ProjectID=ProjectID, Name=Name, Species=Species, ShortName=ShortName))
  }) #later on can use customer uploaded data

output$project <- renderText({
  if (input$sel_project==""){"Please select or upload a date set"} else {ProjectInfo()$Name}
  })

html_geneset<-reactive({
  req(ProjectInfo())
  Species=ProjectInfo()$Species
  string=str_replace(html_geneset0, "human", Species)
 # cat(string, "\n") #debug
  return(string)
})
output$html_geneset=renderUI({
  HTML(html_geneset())
})

html_geneset_hm<-reactive({
  req(ProjectInfo())
  Species=ProjectInfo()$Species
  string=str_replace(html_geneset_hm0, "human", Species)
 # cat(string, "\n") #debug
  return(string)
})
output$html_geneset_hm=renderUI({
  HTML(html_geneset_hm())
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
               RDataFile <- paste("data/",  Pinfo$ProjectID, ".RData", sep = "")
 
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

project_summary<-reactive({
  req(DataReactive())
  DataIn = DataReactive()
  groups=DataIn$groups
  tests=DataIn$tests
  summary=str_c('<style type="text/css">
.disc {
 list-style-type: disc;
}
.square {
 list-style-type: square;
 margin-left: -2em;
 font-size: small
}
</style>',
"<h2>Project ", ProjectInfo()$ShortName, "</h2><br>",
    '<ul class="disc"><li>Species: ', ProjectInfo()$Species, "</li>",
    "<li>Number of Samples: ", nrow(DataIn$MetaData), "</li>",
    "<li>Number of Groups: ", length(groups), " (please see group table below)</li>",  
"<li>Number of Comparison Tests: ", length(tests), "</li>",
'<ul class="square">', paste(str_c("<li>", tests, "</li>"), collapse=""), "</ul></li></ul><br><hr>",
"<h4>Number of Samples in Each Group</h4>")
})
output$summary=renderText(project_summary())

group_info<-reactive({
  DataIn <- DataReactive()
  group_info<-DataIn$MetaData%>%group_by(group)%>%count()
  return(t(group_info))
})
output$group_table=renderTable(group_info(), colnames=F)



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
	DT::datatable(results,  extensions = 'Buttons',
  options = list(
    dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
  	pageLength = 15
  ),rownames= T)
})

output$sample <- DT::renderDataTable({
  meta<-DataReactive()$MetaData%>%dplyr::select(-Order, -ComparePairs)
	DT::datatable(meta,  extensions = 'Buttons',  options = list(
	  dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel', 'pdf', 'print'), pageLength = 15))
	
})

output$data_wide <- DT::renderDataTable({
  data_w<-DataReactive()$data_wide
  data_w=round(data_w*1000)/1000
	DT::datatable(data_w, extensions = c('FixedColumns', 'Buttons'),
  options = list(
  	pageLength = 15,
  	dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
    scrollX = TRUE,
    fixedColumns = list(leftColumns = 1)
  ))
})

output$ProteinGeneName <- DT::renderDataTable({
	DT::datatable(DataReactive()$ProteinGeneName, extensions = 'Buttons', options = list(
	  dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
	  pageLength = 15),rownames= FALSE)
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



