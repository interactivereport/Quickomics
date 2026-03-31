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
library(data.table)

#global reactive values
saved_plots <- reactiveValues()  
saved_table <- reactiveValues() 
samples_excludeM<-reactiveVal() #manually excluded samples
samples_excludeF<-reactiveVal() #samples excluded from filtering on sample attributes
samples_excludeM(""); samples_excludeF("")
attribute_filters<-reactiveVal()
attribute_filters(NULL)
resetComp2Sample<-reactiveVal(); resetComp2Sample(FALSE) #control when to reset the tool to Get Samples from Comparison.,
all_groups <-reactiveVal()
group_order <- reactiveVal()
all_samples <-reactiveVal()
sample_order <- reactiveVal()
all_tests<-reactiveVal()
test_order<-reactiveVal()
all_metadata<-reactiveVal()
MetaData_long <-reactiveVal()
# Initialize the reactiveVal to store the IDs for each test
test_id_lookup <- reactiveVal(list())
upload_message <- reactiveVal()
ProteinGeneNameHeader<- reactiveVal()
exp_unit<-reactiveVal()
#saved_palette <- reactiveVal()
ProjectInfo<-reactiveValues(ProjectID=NULL, Name=NULL, Species=NULL, ShortName=NULL, file1=NULL, file2=NULL, file3=NULL, Path=NULL)
showAlert<-reactiveVal()
plot_pca_control<-reactiveVal(0)
plot_heatmap_control<-reactiveVal(0)
plot_exp_control<-reactiveVal(0)
gsea_control<-reactiveVal(0)
ora_control<-reactiveVal(0)


observeEvent(input$exp_unit, {
  Eu=input$exp_unit; exp_unit(Eu)
  })

observe({
query <- parseQueryString(session$clientData$url_search)
if (!is.null(query[['project']])) {
  ProjectID = query[['project']]
  validate(need(ProjectID %in% saved_projects$ProjectID , message = "Please pass a valid ProjectID from URL."))
  ProjectInfo$ProjectID=ProjectID
  ProjectInfo$Name=saved_projects$Name[saved_projects$ProjectID==ProjectID]
  ProjectInfo$Species=saved_projects$Species[saved_projects$ProjectID==ProjectID]
  ProjectInfo$ShortName=saved_projects$ShortNames[saved_projects$ProjectID==ProjectID]
  ProjectInfo$file1= paste("data/",  ProjectID, ".RData", sep = "")  #data file
  ProjectInfo$file2= paste("networkdata/", ProjectID, "_network.RData", sep = "") #Correlation results
  ProjectInfo$file3= paste("data/wgcna_data/wgcna_", ProjectID, ".RData", sep = "") #wgcna results
}
if (!is.null(query[['unlisted']])) {
  ProjectID = query[['unlisted']]
  validate(need(file.exists(str_c("unlisted/",  ProjectID, ".csv")), 
                message = "Please pass a valid ProjectID from URL. Files must be located in unlisted folder" ))
  unlisted_project=read.csv(str_c("unlisted/", ProjectID, ".csv"))
  ProjectInfo$ProjectID=ProjectID
  ProjectInfo$Name=unlisted_project$Name
  ProjectInfo$Species=unlisted_project$Species
  ProjectInfo$ShortName=unlisted_project$ShortName
  ProjectInfo$file1= paste("unlisted/",  ProjectID, ".RData", sep = "")  #data file
  ProjectInfo$file2= paste("unlisted/", ProjectID, "_network.RData", sep = "") #Correlation results
  ProjectInfo$file3= paste("unlisted/wgcna_", ProjectID, ".RData", sep = "") #wgcna results
  if ("Path" %in% names(unlisted_project)) {ProjectInfo$Path=unlisted_project$Path} 
  if ("ExpressionUnit" %in% names(unlisted_project)) {updateTextInput(session, "exp_unit", value=unlisted_project$ExpressionUnit[1]) }
}
if (!is.null(query[['serverfile']])) {
  ProjectID = query[['serverfile']]
  if (!is.null(server_dir)) {
    validate(need(file.exists(str_c(server_dir, "/",  ProjectID, ".csv")), 
                  message = "Please pass a valid ProjectID from URL. Files must be located in server file folder" ))
    unlisted_project=read.csv(str_c(server_dir, "/",  ProjectID, ".csv"))
    ProjectInfo$ProjectID=ProjectID
    ProjectInfo$Name=unlisted_project$Name
    ProjectInfo$Species=unlisted_project$Species
    ProjectInfo$ShortName=unlisted_project$ShortName
    ProjectInfo$file1= paste(server_dir, "/",   ProjectID, ".RData", sep = "")  #data file
    ProjectInfo$file2= paste(server_dir, "/",  ProjectID, "_network.RData", sep = "") #Correlation results
    ProjectInfo$file3= paste(server_dir, "/",  "wgcna_", ProjectID, ".RData", sep = "") #wgcna results
    if ("Path" %in% names(unlisted_project)) {ProjectInfo$Path=unlisted_project$Path} 
    if ("ExpressionUnit" %in% names(unlisted_project)) {updateTextInput(session, "exp_unit", value=unlisted_project$ExpressionUnit[1]) }
  }
}
if (!is.null(query[['testfile']])) {
  ProjectID = query[['testfile']]
  if (!is.null(test_dir)) {
    validate(need(file.exists(str_c(test_dir, "/",  ProjectID, ".csv")), 
                  message = "Please pass a valid ProjectID from URL. Files must be located in test file folder" ))
    unlisted_project=read.csv(str_c(test_dir, "/",  ProjectID, ".csv"))
    ProjectInfo$ProjectID=ProjectID
    ProjectInfo$Name=unlisted_project$Name
    ProjectInfo$Species=unlisted_project$Species
    ProjectInfo$ShortName=unlisted_project$ShortName
    ProjectInfo$file1= paste(test_dir, "/",   ProjectID, ".RData", sep = "")  #data file
    ProjectInfo$file2= paste(test_dir, "/",  ProjectID, "_network.RData", sep = "") #Correlation results
    ProjectInfo$file3= paste(test_dir, "/",  "wgcna_", ProjectID, ".RData", sep = "") #wgcna results
    if ("Path" %in% names(unlisted_project)) {ProjectInfo$Path=unlisted_project$Path} 
    if ("ExpressionUnit" %in% names(unlisted_project)) {updateTextInput(session, "exp_unit", value=unlisted_project$ExpressionUnit[1]) }
  }
}
})

observe({
if (input$sel_project!="") {
  ProjectID=input$sel_project
  ProjectInfo$ProjectID=ProjectID
  ProjectInfo$Name=saved_projects$Name[saved_projects$ProjectID==ProjectID]
  ProjectInfo$Species=saved_projects$Species[saved_projects$ProjectID==ProjectID]
  ProjectInfo$ShortName=saved_projects$ShortNames[saved_projects$ProjectID==ProjectID]
  ProjectInfo$file1= paste("data/",  ProjectID, ".RData", sep = "")  #data file
  ProjectInfo$file2= paste("networkdata/", ProjectID, "_network.RData", sep = "") #Correlation results
  ProjectInfo$file3= paste("data/wgcna_data/wgcna_", ProjectID, ".RData", sep = "") #wgcna results
  # updateTabsetPanel(session, "Tables", selected = "Sample Table")
}
})

observeEvent(ProjectInfo$ProjectID, {
  #cat("load file UI for", ProjectInfo$ProjectID, "\n")
  updateRadioButtons(session, "heatmap_subset",  selected="All")
  output$gene_highlight_file=renderUI({
    tagList(fileInput("file_gene_highlight", "Highlight Genes (csv with headers like Genes, Pathways, Color)"))
  })
  updateRadioButtons(session, "heatmap_highlight",  selected="No")
  output$gene_annot_file=renderUI({
    tagList(fileInput("file_gene_annot", "Choose gene annotation file (csv with headers like Genes, Pathways, Color)"))
  })
  updateRadioButtons(session, "custom_color",  selected="No")
  output$annot_color_file=renderUI({
    tagList(fileInput("annot_color_file", "Upload annotation Colors (csv with 3 headers: Attribute, Value and Color)"))
  })
  updateTabsetPanel(session, "Tables", selected = "Project Overview")
})

output$project <- renderText({
  if (is.null(ProjectInfo$Name)){"Please select or upload a date set"} else {ProjectInfo$Name}
})

html_geneset<-reactive({
  req(ProjectInfo)
  Species=ProjectInfo$Species
  string=str_replace(html_geneset0, "human", Species)
 #cat(string, "\n") #debug
  return(string)
})
output$html_geneset=renderUI({
  HTML(html_geneset())
})

html_geneset_hm<-reactive({
  req(ProjectInfo)
  Species=ProjectInfo$Species
  string=str_replace(html_geneset_hm0, "human", Species)
  #cat(string, "\n") #debug
  return(string)
})
output$html_geneset_hm=renderUI({
  HTML(html_geneset_hm())
})

html_geneset_exp<-reactive({
  req(ProjectInfo)
  Species=ProjectInfo$Species
  string=str_replace(html_geneset_exp0, "human", Species)
  return(string)
})
output$html_geneset_exp=renderUI({
  HTML(html_geneset_exp())
})


output$ui.action <- renderUI({
  if (is.null(input$file1) ) return()
  tagList(
  textInput("project_name", label="Rename Project", value=input$file1$name),
  radioButtons("species",label="Select species", choices=c("human","mouse", "rat"), inline = F, selected="human"),
  actionButton("customData", "Submit Data")
  )
})


observeEvent(input$customData, {  
  ProjectInfo$ProjectID=str_replace(input$file1$name,  regex(".RData", ignore_case = TRUE), "")
  ProjectInfo$Name=input$project_name
  ProjectInfo$Species=input$species
  ProjectInfo$ShortName=input$project_name
  ProjectInfo$file1=input$file1$datapath; ProjectInfo$file2=input$file2$datapath
  #browser() #debug
})


observe({
  RDataFile <- ProjectInfo$file1
  req(RDataFile)
  objs <- load(RDataFile)
  
  comp_only <- identical(objs, c("results_long", "ProteinGeneName"))
  
  all_tabs <- c("Groups and Samples", "QC_Plots", "Heatmap", "Exp_Plot", "Pattern_Clustering", "time_series", "Correlation_Network", "Correlation", 'wgcna')
  all_data_tables <- c("sample_table", "Result Table", "data_table")
  
  if (comp_only) {
    # hide everything
    lapply(all_tabs, function(t) hideTab("menu", t))
    lapply(all_data_tables, function(t) hideTab("Tables", t))
  } else {
    lapply(all_tabs, function(t) showTab("menu", t))
    lapply(all_data_tables, function(t) showTab("Tables", t))
  }
})

DataReactive <- reactive({
  req(ProjectInfo$ProjectID)
  withProgress(message = 'Fetching data.',
               detail = 'This may take a while...',
               value = 0,
               {
                 RDataFile <- ProjectInfo$file1
                 comp_info <- NULL 
                 objs <- load(RDataFile)
                 if (identical(objs, c("results_long", "ProteinGeneName"))) {
                   tests <- unique(as.character(results_long$test))                     
                   group_names <- NULL
                   MetaData <- NULL
                   data_long <- NULL
                   # ProteinGeneName <- GetProteinGeneNames(ProjectInfo$species)
                   results_long <- results_long %>% 
                     mutate_if(is.factor, as.character) %>% 
                     dplyr::select(UniqueID, test, logFC, P.Value, Adj.P.Value) %>%
                     # mutate(UniqueID = stringr::str_replace(UniqueID, "\\.\\d+$", "")) %>%
                     left_join(ProteinGeneName, by = "UniqueID")
                   data_wide <- NULL
                   data_results <- NULL
                   comp_info <- NULL
                   sel_comp <- NULL
                   all_tests(tests)
                   test_order(tests)
                   ProteinGeneNameHeader(colnames(ProteinGeneName))
                 } else {
                   if (!is.data.frame(data_wide)) {data_wide=data.frame(data_wide, check.names = FALSE)}  #change data_wide to data frame from numeric matrix if needed
                   if (!"Protein.ID" %in% names(ProteinGeneName)) {ProteinGeneName$Protein.ID=NA} #Add Protein.ID column as it is required for certain tools.
                   
                   MetaData_long <- MetaData %>%
                     dplyr::select(-any_of(c("Order", "ComparePairs", "Treatments"))) %>%
                     #dplyr::mutate_if(is.numeric, as.character) %>% #this will fail when there are columns in Time format
                     dplyr::mutate_all(as.character) %>%
                     tidyr::pivot_longer(cols = -sampleid,  names_to = "type",values_to = "group")
                   
                   results_long <-
                     results_long %>% mutate_if(is.factor, as.character)  %>% dplyr::select(UniqueID, test, logFC, P.Value, Adj.P.Value) %>% 
                     # mutate(UniqueID = stringr::str_replace(UniqueID, "\\.\\d+$", "")) %>%
                     left_join(ProteinGeneName, by = "UniqueID")
                   data_long <-
                     data_long %>% mutate_if(is.factor, as.character)  %>% left_join(ProteinGeneName, by = "UniqueID") %>%
                     left_join(MetaData %>% dplyr::select(-any_of(c('group', "Order", "ComparePairs", "Treatments"))), by = "sampleid")
                   
                   group_names <- as.character(unique((MetaData$Order[MetaData$Order != "" & !is.na(MetaData$Order)])))
                   if (length(group_names) == 0) {
                     group_names <- as.character(unique(MetaData$group))
                   }
                   tests  <-
                     as.character(MetaData$ComparePairs[MetaData$ComparePairs != ""])
                   tests<-unique(tests[!is.na(tests)])
                   comp_tests=as.character(unique(results_long$test))
                   if (!all(tests %in% comp_tests) ) { tests <-  gsub("-", "vs", tests) } #for projects where - used in MetaData, "vs" used in results_long
                   if (length(tests) == 0) {
                     tests = unique(as.character(results_long$test))
                   }
                   samples <- as.character( MetaData$sampleid[order(match(MetaData$group,group_names))])
                   group_order(group_names)
                   sample_order(samples)
                   all_samples(samples)
                   all_groups(group_names)
                   all_metadata(MetaData)
                   MetaData_long(MetaData_long)
                   all_tests(tests)
                   test_order(tests)
                   ProteinGeneNameHeader(colnames(ProteinGeneName))
                   sel_comp=NULL
                   # browser() #debug
                   if (!is.null(comp_info)) {
                     # 1. Convert to data.table locally for speed
                     # keep.rownames = "Comparison" preserves your test names
                     sel_comp_dt <- data.table::as.data.table(comp_info, keep.rownames = "Comparison")
                     sel_comp_dt <- sel_comp_dt[!is.na(Group_name) & Group_name != ""]
                     
                     if (nrow(sel_comp_dt) > 0) {
                       # Local data.table version of MetaData to keep the original as a data.frame
                       md_dt <- data.table::as.data.table(MetaData)
                       
                       # Helper function for vectorized processing per comparison row
                       get_test_data <- function(sg_string, group_col, g_test, g_ctrl, md) {
                         
                         # Determine if a subsetting filter actually exists
                         has_subset <- !is.na(sg_string) && sg_string != "" && grepl(":", sg_string)
                         
                         # Logical vector for rows to keep
                         keep <- rep(TRUE, nrow(md))
                         
                         # Step A: Process Subsetting_group if it exists
                         if (has_subset) {
                           conds <- strsplit(sg_string, ";")[[1]]
                           for (cond in conds) {
                             parts <- strsplit(cond, ":")[[1]]
                             # Efficiently update the logical mask
                             keep <- keep & (md[[parts[1]]] == parts[2])
                           }
                           # Assign the comma-separated string of IDs that passed the subset
                           sub_list_val <- paste(md$sampleid[keep], collapse = ",")
                         } else {
                           # PER YOUR REQUEST: If no subsetting string, subset_list is NA
                           sub_list_val <- NA_character_
                         }
                         
                         # Step B: Further filter for Group_test and Group_ctrl for the DEG analysis
                         keep_final <- keep & (md[[group_col]] %in% c(g_test, g_ctrl))
                         final_ids <- md$sampleid[keep_final]
                         
                         return(list(
                           N = length(final_ids),
                           s_list = paste(final_ids, collapse = ","),
                           sub_list = sub_list_val,
                           ids_vec = final_ids 
                         ))
                       }
                       
                       # 2. Execute row-wise using mapply
                       res_list <- mapply(get_test_data, 
                                          sel_comp_dt$Subsetting_group, 
                                          sel_comp_dt$Group_name, 
                                          sel_comp_dt$Group_test, 
                                          sel_comp_dt$Group_ctrl, 
                                          MoreArgs = list(md = md_dt), 
                                          SIMPLIFY = FALSE)
                       
                       # 3. Bind results back to the table
                       # We only bind the first 3 elements (N, s_list, sub_list) to the table
                       res_summary <- data.table::rbindlist(lapply(res_list, function(x) x[1:3]))
                       
                       sel_comp_dt[, `:=`(
                         N_samples = res_summary$N,
                         sample_list = res_summary$s_list,
                         subset_list = res_summary$sub_list
                       )]
                       
                       # 4. Update the reactiveVal for direct ID access in other parts of the app
                       id_mapping <- lapply(res_list, `[[`, "ids_vec")
                       names(id_mapping) <- sel_comp_dt$Comparison
                       test_id_lookup(id_mapping)
                       
                       # 5. Filter out empty tests and convert back to data.frame for the rest of the app
                       sel_comp <- as.data.frame(sel_comp_dt[N_samples > 0])
                       if (nrow(sel_comp) == 0) sel_comp <- NULL
                     }
                   } else if (length(tests) > 0) {
                     # 1. Create the base data.frame using data.table for speed
                     # We use tstrsplit to split "GroupA vs GroupB" into two columns at once
                     sel_comp_dt <- data.table::data.table(Comparison = tests)
                     sel_comp_dt[, c("Group_test", "Group_ctrl") := data.table::tstrsplit(Comparison, "vs", fixed = TRUE)]
                     
                     # 2. Set Row Names (converting to data.frame at the end)
                     # 3. Prepare for Sample ID lookup
                     # Ensure we have a local data.table of MetaData for fast subsetting
                     md_dt <- data.table::as.data.table(MetaData)
                     
                     # 4. Define the lookup and list-building function
                     # This uses the 'group' column from MetaData (adjust if your column name differs)
                     get_samples_for_test <- function(g_test, g_ctrl, md) {
                       # Get vector of IDs for each side
                       ids_test <- md[group == g_test, as.character(sampleid)]
                       ids_ctrl <- md[group == g_ctrl, as.character(sampleid)]
                       
                       combined_vec <- c(ids_test, ids_ctrl)
                       
                       return(list(
                         n = length(combined_vec),
                         s_list = paste(combined_vec, collapse = ","),
                         ids_vec = combined_vec
                       ))
                     }
                     
                     # 5. Execute lookup
                     res_list <- mapply(get_samples_for_test, 
                                        sel_comp_dt$Group_test, 
                                        sel_comp_dt$Group_ctrl, 
                                        MoreArgs = list(md = md_dt), 
                                        SIMPLIFY = FALSE)
                     
                     # 6. Bind results back to the table
                     res_summary <- data.table::rbindlist(lapply(res_list, function(x) x[1:2]))
                     sel_comp_dt[, `:=`(
                       N_samples = res_summary$n,
                       sample_list = res_summary$s_list,
                       Group_name = "", # Defaulting to 'group' column
                       subset_list = NA_character_
                     )]
                     
                     # 7. Update the reactiveVal with raw vectors for instant access
                     id_mapping <- lapply(res_list, `[[`, "ids_vec")
                     names(id_mapping) <- sel_comp_dt$Comparison
                     
                     test_id_lookup(id_mapping)
                     
                     # 8. Final conversion to data.frame with Row Names
                     sel_comp <- as.data.frame(sel_comp_dt)
                     rownames(sel_comp) <- sel_comp$Comparison
                   }
                 }
                 return(
                   list(
                     "groups" = group_names,
                     "MetaData" = MetaData,
                     "results_long" = results_long,
                     "data_long" = data_long,
                     "ProteinGeneName" = ProteinGeneName,
                     "data_wide" = data_wide,
                     "data_results" = data_results,
                     "tests" = tests,
                     "comp_info"=comp_info,
                     "sel_comp"=sel_comp
                   )
                 )
               })
  
})

observeEvent(DataReactive(), {
  req(DataReactive()$MetaData)
  plot_pca_control(plot_pca_control()+1)
  plot_heatmap_control( plot_heatmap_control()+1)
  plot_exp_control(plot_exp_control()+1)
})

observeEvent(DataReactive(), {
  req(DataReactive()$results_long)
  gsea_control(gsea_control()+1)
  ora_control(ora_control()+1)
})

project_summary<-reactive({
  req(DataReactive())
  DataIn = DataReactive()
  restricted_msg <- if (!public_dataset) {
    "<h4>This project is loaded in restricted mode. Individual sample data have been masked for privacy.</h4><br>"
  } else {
    ""
  }
  
  if (is.null(DataIn$MetaData)) {
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
                  "<h2>Project ", ProjectInfo$ShortName, "</h2><br>",
                  restricted_msg,
                  '<ul class="disc"><li>Species: ', ProjectInfo$Species, "</li>",
                  "<li>Description: ", ProjectInfo$Name, "</li>",
                  "<li>Data Path: ", ProjectInfo$Path, "</li>",
                  "<li>This is a project only has comparison results</li>",
                  "<li>Number of Comparison Tests: ", length(tests), "</li>",
                  '<ul class="square">', paste(str_c("<li>", tests, "</li>"), collapse=""), "</ul></li></ul><br><hr>")
  } else {
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
                  "<h2>Project ", ProjectInfo$ShortName, "</h2><br>",
                  restricted_msg,
                  '<ul class="disc"><li>Species: ', ProjectInfo$Species, "</li>",
                  "<li>Description: ", ProjectInfo$Name, "</li>",
                  "<li>Data Path: ", ProjectInfo$Path, "</li>",
                  "<li>Number of Samples: ", nrow(DataIn$MetaData), "</li>",
                  "<li>Number of Groups: ", length(groups), " (please see group table below)</li>",  
                  "<li>Number of Genes/Proteins: ", nrow(DataIn$data_wide), "</li>",
                  "<li>Number of Comparison Tests: ", length(tests), "</li>",
                  '<ul class="square">', paste(str_c("<li>", tests, "</li>"), collapse=""), "</ul></li></ul><br><hr>",
                  "<h4>Number of Samples in Each Group</h4>")
    
  }
})
    
    output$summary=renderText(project_summary())
    
    group_info<-reactive({
      req(DataReactive()$MetaData)
      DataIn <- DataReactive()
      group_info<-DataIn$MetaData%>%group_by(group)%>%dplyr::count()
      #browser() #bebug
      return(t(group_info))
    })
    
    output$group_table <- renderTable({
      req(group_info())   # stops if NULL
      group_info()
    }, colnames = FALSE)
    
    
    output$results <- DT::renderDataTable({
      DataIn <- DataReactive()
      req(DataIn$data_results) 
      results <- DataIn$data_results %>%
        dplyr::select(-one_of(c("Fasta.headers","UniqueID","id")))
      results[,sapply(results,is.numeric)] <- signif(results[,sapply(results,is.numeric)],3)
      DT::datatable(results,  extensions = 'Buttons',
                    options = list(
                      dom = 'lBfrtip', buttons = c('csv', 'excel', 'print'),
                      pageLength = 15
                    ),rownames= T)
    })
    
    output$sample <-  DT::renderDT(server=FALSE,{
      req(public_dataset)
      req(DataReactive()$MetaData)
      meta<-DataReactive()$MetaData%>%dplyr::select(-Order, -ComparePairs)
      DT::datatable(meta,  extensions = 'Buttons',  options = list(
        dom = 'lBfrtip', pageLength = 15,
        buttons = list(
          list(extend = "csv", text = "Download Page", filename = "Page_Samples",
               exportOptions = list(modifier = list(page = "current"))),
          list(extend = "csv", text = "Download All", filename = "All_Samples",
               exportOptions = list(modifier = list(page = "all")))
        )
      ), rownames= F)
    })
    
    observe({
      if (public_dataset) {
        showTab(inputId = "Tables", target = "sample_table")
      } else {
        hideTab(inputId = "Tables", target = "sample_table")
      }
    })
    
    output$comp_info <- renderUI ({
      if (is.null(DataReactive()$comp_info)) return()
      output$comparison <-  DT::renderDT(server=FALSE,{
        DT::datatable(DataReactive()$comp_info,  extensions = 'Buttons',  options = list(
          dom = 'lBfrtip', pageLength = 15,
          buttons = list(
            list(extend = "csv", text = "Download Page", filename = "Page_results",
                 exportOptions = list(modifier = list(page = "current"))),
            list(extend = "csv", text = "Download All", filename = "All_Results",
                 exportOptions = list(modifier = list(page = "all")))
          )
        )
        )
      })
      tagList(
        h4("Comparison Table (shown only when RData file contains comp_info)"),
        dataTableOutput('comparison')
      )
    })
    
    output$data_wide <- DT::renderDataTable({
      req(public_dataset)
      data_w<-DataReactive()$data_wide
      req(data_w)
      data_w=round(data_w*1000)/1000
      DT::datatable(data_w, extensions = c('FixedColumns', 'Buttons'),
                    options = list(
                      pageLength = 15,
                      dom = 'lBfrtip', buttons = c('csv', 'excel', 'print'),
                      scrollX = TRUE,
                      fixedColumns = list(leftColumns = 1)
                    ))
    })
    
    observe({
      if (public_dataset) {
        showTab(inputId = "Tables", target = "data_table")
      } else {
        hideTab(inputId = "Tables", target = "data_table")
      }
    })
    
    output$ProteinGeneName <- DT::renderDataTable({
      if (is.null(DataReactive()$ProteinGeneName)) return()
      DT::datatable(DataReactive()$ProteinGeneName, extensions = 'Buttons', options = list(
        dom = 'lBfrtip', buttons = c('csv', 'excel', 'print'),
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
    
    
    
    
    
    
    
    
    
    
    
    