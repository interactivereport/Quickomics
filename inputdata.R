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
saved_gcts <- reactiveValues() 

samples_excludeM<-reactiveVal() #manually excluded samples
samples_excludeF<-reactiveVal() #samples excluded from filtering on sample attributes
samples_excludeM(""); samples_excludeF("")
attribute_filters<- reactiveVal()
attribute_filters(NULL)
numerical_attributes <- reactiveVal()
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
ProjectInfo<-reactiveValues(ProjectID=NULL, Name=NULL, Species=NULL, ShortName=NULL, file1=NULL, file2=NULL, file3=NULL, Path=NULL, source=NULL, SourceLabel=NULL)
# Wall-clock bookkeeping for the "[timing]" console messages: `start` is set
# the moment a project ID is input (URL parameter, project dropdown, or
# upload) and consumed once when DataReactive() finishes for that project.
load_timing <- new.env(parent = emptyenv())
load_timing$start <- NULL
load_timing$id <- NULL
mark_project_start <- function(ProjectID) {
  load_timing$start <- Sys.time()
  load_timing$id <- ProjectID
}
showAlert<-reactiveVal()
plot_pca_control<-reactiveVal(0)
plot_heatmap_control<-reactiveVal(0)
plot_exp_control<-reactiveVal(0)
gsea_control<-reactiveVal(0)
ora_control<-reactiveVal(0)

observeEvent(input$exp_unit, {
  Eu=input$exp_unit; exp_unit(Eu)
})

# Loads a project from one of the raw URL query-string parameters below.
# Factored out (rather than left inline in the observe() below) so
# app.R's onRestored() can also call it directly -- project/unlisted/
# serverfile/testfile/api_serverfile/api_testfile/datalake/auto select the project straight from
# session$clientData$url_search rather than from a Shiny input, so they're
# never captured by Shiny's bookmarking and are gone from the URL once
# Save Session's doBookmark() replaces it with "?_state_id_=...". Without
# replaying this same lookup at restore time from a value captured at save
# time, a session whose project was opened via one of these parameters
# restores to no project at all.
#
# Uses showNotification()+return() rather than validate()/need() for
# invalid input, on purpose: this function is called both from a plain
# observe() below (which has no render output of its own, so validate()
# would fail completely silently there) AND from onRestored() in app.R, a
# non-reactive context where an uncaught validate() condition surfaces as a
# raw error instead of a clean message. showNotification() works correctly
# in both places.
#
# With auto = TRUE (used only by the ?auto= key) the unlisted/serverfile/testfile/
# api_* /datalake loaders stay silent on "not found" and API errors and just
# return FALSE, so the next source in auto_source_priority can be tried.
load_project_from_query <- function(key, ProjectID, auto = FALSE) {
  if (!auto) mark_project_start(ProjectID)
  log_debug("load_project_from_query: key=", key, " ProjectID=", ProjectID, if (auto) " (auto)")
  if (key == 'project') {
    if (!(ProjectID %in% saved_projects$ProjectID)) {
      log_warn("Project load failed (key=project): '", ProjectID, "' not found in saved_projects.")
      showNotification(
        paste0("Project '", ProjectID, "' could not be found. Please check the project ID in the URL."),
        type = "error", duration = NULL
      )
      return()
    }
    ProjectInfo$ProjectID=ProjectID
    ProjectInfo$Name=saved_projects$Name[saved_projects$ProjectID==ProjectID]
    ProjectInfo$Species=saved_projects$Species[saved_projects$ProjectID==ProjectID]
    ProjectInfo$ShortName=saved_projects$ShortNames[saved_projects$ProjectID==ProjectID]
    ProjectInfo$file1= paste("data/",  ProjectID, ".RData", sep = "")  #data file
    ProjectInfo$file2= paste("networkdata/", ProjectID, "_network.RData", sep = "") #Correlation results
    ProjectInfo$file3= paste("data/wgcna_data/wgcna_", ProjectID, ".RData", sep = "") #wgcna results
    ProjectInfo$source=NULL  #local data/ folder project, not API-backed
    project_query_source(list(key = key, value = ProjectID))
    set_project_source(key)
    log_info("Project loaded: ", ProjectID, " (source=project)")
  } else if (key == 'unlisted') {
    if (!file.exists(str_c("unlisted/",  ProjectID, ".csv"))) {
      log_warn("Project load failed (key=unlisted): '", ProjectID, "' -- no unlisted/", ProjectID, ".csv")
      if (!auto) showNotification(
        paste0("Project '", ProjectID, "' could not be found in the unlisted folder. Please check the project ID in the URL."),
        type = "error", duration = NULL
      )
      return(FALSE)
    }
    unlisted_project=read.csv(str_c("unlisted/", ProjectID, ".csv"))
    ProjectInfo$ProjectID=ProjectID
    ProjectInfo$Name=unlisted_project$Name
    ProjectInfo$Species=unlisted_project$Species
    ProjectInfo$ShortName=unlisted_project$ShortName
    ProjectInfo$file1= paste("unlisted/",  ProjectID, ".RData", sep = "")  #data file
    ProjectInfo$file2= paste("unlisted/", ProjectID, "_network.RData", sep = "") #Correlation results
    ProjectInfo$file3= paste("unlisted/wgcna_", ProjectID, ".RData", sep = "") #wgcna results
    ProjectInfo$source=NULL  #unlisted/ folder project, not API-backed
    if ("Path" %in% names(unlisted_project)) {ProjectInfo$Path=unlisted_project$Path}
    if ("ExpressionUnit" %in% names(unlisted_project)) {updateTextInput(session, "exp_unit", value=unlisted_project$ExpressionUnit[1]) }
    project_query_source(list(key = key, value = ProjectID))
    if (!auto) set_project_source(key)
    log_info("Project loaded: ", ProjectID, " (source=unlisted)")
    return(TRUE)
  } else if (key == 'serverfile') {
    # Reads {server_dir}/{ProjectID}.RData directly with base::load().
    load_project_from_dir(key, ProjectID, server_dir, "server_dir", auto)
  } else if (key == 'testfile') {
    # Reads {test_dir}/{ProjectID}.RData directly with base::load().
    load_project_from_dir(key, ProjectID, test_dir, "test_dir", auto)
  } else if (key == 'api_serverfile') {
    # Same published projects as serverfile, fetched through the Quickomics API.
    load_project_from_api(key, ProjectID, "published", auto)
  } else if (key == 'api_testfile') {
    # Same unpublished projects as testfile, fetched through the Quickomics API.
    load_project_from_api(key, ProjectID, "unpublished", auto)
  } else if (key == 'datalake') {
    # Datalake projects are only reachable through the Quickomics API.
    load_project_from_api(key, ProjectID, "datalake", auto)
  } else if (key == 'auto') {
    # Tries each source in auto_source_priority and loads from the first one
    # that has the project. Save Session records the key that worked (via
    # project_query_source() in the loader), not 'auto', so a restored
    # session reopens from the same source.
    for (k in auto_source_priority) {
      if (isTRUE(load_project_from_query(k, ProjectID, auto = TRUE))) {
        set_project_source(k, via_auto = TRUE)
        return(invisible(TRUE))
      }
    }
    log_warn("Project load failed (key=auto): '", ProjectID, "' not found in any source (tried ",
             paste(auto_source_priority, collapse = ", "), ").")
    showNotification(
      paste0("Project '", ProjectID, "' could not be found in any data source (tried ",
             paste(auto_source_priority, collapse = ", "), "). Please check the project ID in the URL."),
      type = "error", duration = NULL
    )
    return(invisible(FALSE))
  }
}

# Order the ?auto= key tries sources in, highest priority first.
auto_source_priority <- c('datalake', 'api_serverfile', 'serverfile', 'api_testfile', 'testfile', 'unlisted')

# Where a project's data came from, shown on the Project Overview page
# (ProjectInfo$SourceLabel) and in a notification when the project loads.
project_source_labels <- c(
  project        = "Local data folder (project)",
  unlisted       = "Unlisted folder (unlisted)",
  serverfile     = "Server folder, read directly (serverfile)",
  testfile       = "Test folder, read directly (testfile)",
  api_serverfile = "QuickOmics API, published projects (api_serverfile)",
  api_testfile   = "QuickOmics API, unpublished projects (api_testfile)",
  datalake       = "Data lake, through the QuickOmics API (datalake)"
)

set_project_source <- function(key, via_auto = FALSE) {
  label <- unname(project_source_labels[key])
  if (via_auto) label <- paste0(label, ", chosen automatically")
  ProjectInfo$SourceLabel <- label
  showNotification(paste0("Project data source: ", label), type = "message", duration = 10)
}

# serverfile/testfile: project registered in a server-side folder as
# {dir}/{ProjectID}.csv, with its data in {dir}/{ProjectID}.RData.
load_project_from_dir <- function(key, ProjectID, dir, dir_name, auto = FALSE) {
  if (is.null(dir)) {
    log_warn("Project load failed (key=", key, "): ", dir_name, " is not configured.")
    if (!auto) showNotification(
      paste0("Cannot open '", ProjectID, "': ", dir_name, " is not configured on this server."),
      type = "error", duration = NULL
    )
    return(FALSE)
  }
  if (!file.exists(str_c(dir, "/",  ProjectID, ".csv"))) {
    log_warn("Project load failed (key=", key, "): '", ProjectID, "' -- no ", dir, "/", ProjectID, ".csv")
    if (!auto) showNotification(
      paste0("Project '", ProjectID, "' could not be found. Please check the project ID in the URL."),
      type = "error", duration = NULL
    )
    return(FALSE)
  }
  unlisted_project=read.csv(str_c(dir, "/",  ProjectID, ".csv"))
  ProjectInfo$ProjectID=ProjectID
  ProjectInfo$Name=unlisted_project$Name
  ProjectInfo$Species=unlisted_project$Species
  ProjectInfo$ShortName=unlisted_project$ShortName
  ProjectInfo$file1= paste(dir, "/",   ProjectID, ".RData", sep = "")  #data file
  ProjectInfo$file2= paste(dir, "/",  ProjectID, "_network.RData", sep = "") #Correlation results
  ProjectInfo$file3= paste(dir, "/",  "wgcna_", ProjectID, ".RData", sep = "") #wgcna results
  ProjectInfo$source=NULL  #read directly from disk, not API-backed
  if ("Path" %in% names(unlisted_project)) {ProjectInfo$Path=unlisted_project$Path}
  if ("ExpressionUnit" %in% names(unlisted_project)) {updateTextInput(session, "exp_unit", value=unlisted_project$ExpressionUnit[1]) }
  project_query_source(list(key = key, value = ProjectID))
  if (!auto) set_project_source(key)
  log_info("Project loaded: ", ProjectID, " (source=", key, ")")
  TRUE
}

# api_serverfile/api_testfile/datalake: project fetched through the
# Quickomics API, where `source` is the API's {source} path segment
# ("published", "unpublished" or "datalake").
load_project_from_api <- function(key, ProjectID, source, auto = FALSE) {
  api_result <- tryCatch(
    list(ok = TRUE, meta = api_get_project_meta(ProjectID, source)),
    error = function(e) list(ok = FALSE, error = conditionMessage(e))
  )
  if (!api_result$ok) {
    log_warn("Project load failed (key=", key, "): Quickomics API error for '", ProjectID, "': ", api_result$error)
    if (!auto) showNotification(
      paste("Quickomics API is not working:", api_result$error),
      type = "error", duration = NULL
    )
    return(FALSE)
  }
  meta <- api_result$meta
  if (is.null(meta)) {
    log_warn("Project load failed (key=", key, "): '", ProjectID, "' not found in the Quickomics API (source=", source, ").")
    if (!auto) showNotification(
      paste0("Project '", ProjectID, "' could not be found. Please check the project ID in the URL."),
      type = "error", duration = NULL
    )
    return(FALSE)
  }
  ProjectInfo$ProjectID=ProjectID
  ProjectInfo$Name=meta$Name
  ProjectInfo$Species=meta$Species
  ProjectInfo$ShortName=meta$ShortName
  # file1/file2/file3 are placeholders, not real paths, for API-sourced
  # projects -- actual data loading goes through api_load()/
  # api_load_network()/api_load_wgcna() (see below, network.R, wgcna.R),
  # none of which dereference these strings as real file paths. They're
  # kept non-NULL only because other code gates on `req(ProjectInfo$file1)`.
  ProjectInfo$file1= paste0("api://", source, "/", ProjectID, ".RData")
  ProjectInfo$file2= paste0("api://", source, "/", ProjectID, "_network.RData")
  ProjectInfo$file3= paste0("api://", source, "/wgcna_", ProjectID, ".RData")
  ProjectInfo$source=source  #routes file1/network/wgcna loading through the Quickomics API
  if (!is.null(meta$Path)) {ProjectInfo$Path=meta$Path}
  if (!is.null(meta$ExpressionUnit)) {updateTextInput(session, "exp_unit", value=meta$ExpressionUnit)}
  project_query_source(list(key = key, value = ProjectID))
  if (!auto) set_project_source(key)
  log_info("Project loaded: ", ProjectID, " (source=", key, ")")
  TRUE
}

# URL query keys that select a project, in the order the observe() below
# checks them.
project_query_keys <- c('project', 'unlisted', 'serverfile', 'testfile',
                        'api_serverfile', 'api_testfile', 'datalake', 'auto')

observe({
  query <- parseQueryString(session$clientData$url_search)
  # No early exit on a match -- matches the original behavior where all
  # blocks ran unconditionally, so if more than one of these parameters
  # were somehow present at once, the last one in this vector (auto) wins.
  for (key in project_query_keys) {
    if (!is.null(query[[key]])) {
      load_project_from_query(key, query[[key]])
    }
  }
})

observe({
  if (input$sel_project!="") {
    ProjectID=input$sel_project
    mark_project_start(ProjectID)
    ProjectInfo$ProjectID=ProjectID
    ProjectInfo$Name=saved_projects$Name[saved_projects$ProjectID==ProjectID]
    ProjectInfo$Species=saved_projects$Species[saved_projects$ProjectID==ProjectID]
    ProjectInfo$ShortName=saved_projects$ShortNames[saved_projects$ProjectID==ProjectID]
    ProjectInfo$file1= paste("data/",  ProjectID, ".RData", sep = "")  #data file
    ProjectInfo$file2= paste("networkdata/", ProjectID, "_network.RData", sep = "") #Correlation results
    ProjectInfo$file3= paste("data/wgcna_data/wgcna_", ProjectID, ".RData", sep = "") #wgcna results
    ProjectInfo$source=NULL  #local data/ folder project, not API-backed
    ProjectInfo$SourceLabel="Local data folder (project list)"
    # updateTabsetPanel(session, "Tables", selected = "Sample Table")
    log_info("Project loaded: ", ProjectID, " (source=sel_project dropdown)")
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
  mark_project_start(ProjectInfo$ProjectID)
  ProjectInfo$Name=input$project_name
  ProjectInfo$Species=input$species
  ProjectInfo$ShortName=input$project_name
  ProjectInfo$file1=input$file1$datapath; ProjectInfo$file2=input$file2$datapath
  ProjectInfo$source=NULL  #user-uploaded file, not API-backed
  ProjectInfo$SourceLabel="Uploaded file"
  #browser() #debug
})


# One API download per project load: the tab-visibility observer below and
# DataReactive() both need the project's bundle, so whichever runs first
# downloads it and the other reuses it. It is released once both have used
# it, so the session doesn't keep a second copy of the data around.
api_bundle_holder <- new.env(parent = emptyenv())
get_api_bundle <- function(user) {
  key <- paste0(ProjectInfo$source, "/", ProjectInfo$ProjectID)
  h <- api_bundle_holder
  if (!identical(h$key, key) || is.null(h$bundle)) {
    h$key <- key
    h$bundle <- api_fetch_bundle(ProjectInfo$ProjectID, ProjectInfo$source)
    h$users <- character(0)
  }
  bundle <- h$bundle
  h$users <- union(h$users, user)
  if (all(c("tabs", "data") %in% h$users)) h$bundle <- NULL
  bundle
}

observe({
  RDataFile <- ProjectInfo$file1
  req(RDataFile)
  # ProjectInfo$source is set ("published"/"unpublished"/"datalake") only for
  # projects loaded via ?serverfile=/?testfile=/?datalake=; everything else
  # (local data/ folder, unlisted/ folder, file upload) still uses
  # base::load() unchanged.
  observe_start <- Sys.time()
  if (!is.null(ProjectInfo$source)) {
    objs <- api_load(ProjectInfo$ProjectID, ProjectInfo$source, bundle = get_api_bundle("tabs"))
  } else {
    objs <- load(RDataFile)
  }
  message(sprintf(
    "[timing] tab-visibility observe load for '%s' (source=%s) took %.2f sec",
    ProjectInfo$ProjectID,
    if (is.null(ProjectInfo$source)) "local" else ProjectInfo$source,
    as.numeric(difftime(Sys.time(), observe_start, units = "secs"))
  ))
  
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
                 reactive_entry <- Sys.time()
                 RDataFile <- ProjectInfo$file1
                 comp_info <- NULL
                 # ProjectInfo$source is set ("published"/"unpublished"/
                 # "datalake") only for projects loaded via
                 # ?serverfile=/?testfile=/?datalake=; every other source
                 # (local data/ folder, unlisted/ folder, file upload) still
                 # uses base::load() completely unchanged.
                 fetch_start <- Sys.time()
                 if (!is.null(ProjectInfo$source)) {
                   objs <- api_load(ProjectInfo$ProjectID, ProjectInfo$source, bundle = get_api_bundle("data"))
                 } else {
                   objs <- load(RDataFile)
                 }
                 fetch_secs <- as.numeric(difftime(Sys.time(), fetch_start, units = "secs"))
                 message(sprintf(
                   "[timing] DataReactive load for '%s' (source=%s) took %.2f sec",
                   ProjectInfo$ProjectID,
                   if (is.null(ProjectInfo$source)) "local" else ProjectInfo$source,
                   fetch_secs
                 ))
                 if (identical(objs, c("results_long", "ProteinGeneName"))) {
                   tests <- unique(as.character(results_long$test))                     
                   group_names <- NULL
                   MetaData <- NULL
                   data_long <- NULL
                   # ProteinGeneName <- GetProteinGeneNames(ProjectInfo$species)
                   # Vectorized data.table alternative
                   results_long <- data.table::as.data.table(results_long)
                   # Convert factors to character
                   for (col in names(results_long)) {
                     if (is.factor(results_long[[col]])) {
                       results_long[[col]] <- as.character(results_long[[col]])
                     }
                   }
                   # Select specific columns
                   results_long <- results_long[, .(UniqueID, test, logFC, P.Value, Adj.P.Value)]
                   # Left join with ProteinGeneName using data.table
                   pgn_dt <- data.table::as.data.table(ProteinGeneName)
                   data.table::setkey(pgn_dt, UniqueID)
                   results_long <- pgn_dt[results_long]
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
                   if (!is.character(MetaData$sampleid)) {
                     MetaData$sampleid <- as.character(MetaData$sampleid)
                   } 
                   if (!is.character(data_long$sampleid)) {
                     data_long$sampleid <- as.character(data_long$sampleid)
                   }
                   if (!is.character(data_results$id)) {
                     data_results$id <- as.character(data_results$id)
                   }
                   # Vectorized data.table alternative to dplyr
                   DT_MetaData_long <- data.table::as.data.table(MetaData)
                   cols_to_drop <- c("Order", "ComparePairs", "Treatments")
                   cols_to_drop <- cols_to_drop[cols_to_drop %in% names(DT_MetaData_long)]
                   if (length(cols_to_drop) > 0) {
                     DT_MetaData_long[, (cols_to_drop) := NULL]
                   }
                   # Convert all columns to character
                   DT_MetaData_long <- DT_MetaData_long[, lapply(.SD, as.character)]
                   # Pivot longer (long format)
                   DT_MetaData_long <- data.table::melt(DT_MetaData_long, 
                                                        id.vars = "sampleid",
                                                        variable.name = "type",
                                                        value.name = "group")
                   # Missing values show up as an "NA" tile in Groups and Samples,
                   # and the tile hands back the string "NA"; store them as "NA"
                   # so that tile matches its samples (true NA never matches).
                   DT_MetaData_long[is.na(group), group := "NA"]
                   
                   # Returns names of columns where is.numeric is TRUE
                   num_cols <- names(MetaData)[sapply(MetaData, is.numeric)]
                   
                   # Vectorized data.table alternative
                   results_long <- data.table::as.data.table(results_long)
                   # Convert factors to character
                   for (col in names(results_long)) {
                     if (is.factor(results_long[[col]])) {
                       results_long[[col]] <- as.character(results_long[[col]])
                     }
                   }
                   # Select specific columns
                   results_long <- results_long[, .(UniqueID, test, logFC, P.Value, Adj.P.Value)]
                   # Left join with ProteinGeneName using data.table
                   pgn_dt <- data.table::as.data.table(ProteinGeneName)
                   data.table::setkey(pgn_dt, UniqueID)
                   results_long <- pgn_dt[results_long]
                   # Vectorized data.table alternative
                   data_long <- data.table::as.data.table(data_long)
                   # Convert factors to character
                   for (col in names(data_long)) {
                     if (is.factor(data_long[[col]])) {
                       data_long[[col]] <- as.character(data_long[[col]])
                     }
                   }
                   # Ensure sampleid is character for join compatibility
                   if (is.numeric(data_long$sampleid)) {
                     data_long[, sampleid := as.character(sampleid)]
                   }
                   
                   # Left join with ProteinGeneName
                   # pgn_dt <- data.table::as.data.table(ProteinGeneName)
                   # data.table::setkey(pgn_dt, UniqueID)
                   data_long <- pgn_dt[data_long]
                   
                   # Left join with MetaData (drop unwanted columns first)
                   md_dt <- data.table::as.data.table(MetaData)
                   # Ensure sampleid is character for join compatibility
                   if (is.numeric(md_dt$sampleid)) {
                     md_dt[, sampleid := as.character(sampleid)]
                   }
                   # Drop unwanted columns (only if they exist)
                   cols_to_drop_md <- c("group", "Order", "ComparePairs", "Treatments")
                   cols_to_drop_md <- cols_to_drop_md[cols_to_drop_md %in% names(md_dt)]
                   if (length(cols_to_drop_md) > 0) {
                     md_dt[, (cols_to_drop_md) := NULL]
                   }
                   data.table::setkey(md_dt, sampleid)
                   data_long <- md_dt[data_long, on = .(sampleid)]
                   
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
                   MetaData_long(DT_MetaData_long)
                   numerical_attributes(num_cols)
                   all_tests(tests)
                   test_order(tests)
                   ProteinGeneNameHeader(colnames(ProteinGeneName))
                   sel_comp=NULL
                   # browser() #debug
                   # Start every project load with an empty lookup so a previous
                   # project's comparisons can't leak into this one.
                   test_id_lookup(list())
                   # comp_info is only usable when at least one comparison has
                   # Group_name, Group_test and Group_ctrl filled in. Otherwise (e.g.
                   # all NA, as in TargetALS_Plasma) fall back to splitting the
                   # comparison name at "vs", same as when comp_info is missing.
                   filled <- function(x) !is.na(x) & trimws(as.character(x)) != ""
                   comp_usable <- !is.null(comp_info) &&
                     all(c("Group_name", "Group_test", "Group_ctrl") %in% colnames(comp_info)) &&
                     any(filled(comp_info$Group_name) & filled(comp_info$Group_test) & filled(comp_info$Group_ctrl))
                   if (!is.null(comp_info) && !comp_usable) {
                     log_info("comp_info has no usable Group_name/Group_test/Group_ctrl values; matching comparisons to samples by splitting test names at 'vs'.")
                   }
                   if (comp_usable) {
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
                     # 1. Split each test name into its test and control group, e.g.
                     # "A vs B", "A_vs_B", "A.vs.B". Separators are tried from most to
                     # least specific; a split counts only when both halves are values
                     # of MetaData$group. Tests that can't be split get NA groups (no
                     # samples, no lookup entry), so Groups and Samples keeps them.
                     known_groups <- unique(as.character(MetaData$group))
                     split_test_name <- function(test_name) {
                       for (sep in c("_vs_", ".vs.", " vs ", "-vs-", "vs")) {
                         parts <- strsplit(test_name, sep, fixed = TRUE)[[1]]
                         if (length(parts) != 2) next
                         parts <- gsub("^[[:space:]_.-]+|[[:space:]_.-]+$", "", parts)
                         if (all(parts %in% known_groups)) return(parts)
                       }
                       c(NA_character_, NA_character_)
                     }
                     split_parts <- lapply(as.character(tests), split_test_name)
                     sel_comp_dt <- data.table::data.table(
                       Comparison = tests,
                       Group_test = vapply(split_parts, `[`, character(1), 1),
                       Group_ctrl = vapply(split_parts, `[`, character(1), 2)
                     )
                     
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
                     
                     # 7. Update the reactiveVal with raw vectors for instant access.
                     # Only tests matched to samples get an entry.
                     id_mapping <- lapply(res_list, `[[`, "ids_vec")
                     names(id_mapping) <- sel_comp_dt$Comparison
                     id_mapping <- id_mapping[lengths(id_mapping) > 0]
                     
                     test_id_lookup(id_mapping)
                     
                     # 8. Final conversion to data.frame with Row Names, dropping
                     # tests with no matched samples
                     sel_comp <- as.data.frame(sel_comp_dt[N_samples > 0])
                     if (nrow(sel_comp) == 0) {
                       sel_comp <- NULL
                     } else {
                       rownames(sel_comp) <- sel_comp$Comparison
                     }
                   }
                 }
                 if (!is.null(load_timing$start) && identical(load_timing$id, ProjectInfo$ProjectID)) {
                   done <- Sys.time()
                   message(sprintf(
                     "[timing] '%s' (source=%s): %.2f sec from project ID input to DataReactive() ready = %.2f before DataReactive started + %.2f load + %.2f transform",
                     ProjectInfo$ProjectID,
                     if (is.null(ProjectInfo$source)) "local" else ProjectInfo$source,
                     as.numeric(difftime(done, load_timing$start, units = "secs")),
                     as.numeric(difftime(reactive_entry, load_timing$start, units = "secs")),
                     fetch_secs,
                     as.numeric(difftime(done, fetch_start, units = "secs")) - fetch_secs
                   ))
                   load_timing$start <- NULL  # report once per project selection
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
                  "<li>Data Source: ", ProjectInfo$SourceLabel, "</li>",
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
                  "<li>Data Source: ", ProjectInfo$SourceLabel, "</li>",
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
    dplyr::select(-dplyr::any_of(c("Fasta.headers","UniqueID","id")))
  results[,sapply(results,is.numeric)] <- signif(results[,sapply(results,is.numeric)],3)
  DT::datatable(results,  extensions = 'Buttons',
                options = list(
                  dom = 'lBfrtip', buttons = list(list(extend = "csv", text = "Download Page", filename = "Page_Results", exportOptions = list(modifier = list(page = "current")))),
                  pageLength = 15
                ),rownames= T)
})

output$sample <-  DT::renderDT(server=FALSE,{
  req(public_dataset)
  req(DataReactive()$MetaData)
  meta<-DataReactive()$MetaData%>%dplyr::select(-dplyr::any_of(c("Order", "ComparePairs")))
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
                  dom = 'lBfrtip', buttons = list(list(extend = "csv", text = "Download Page", filename = "Page_Data", exportOptions = list(modifier = list(page = "current")))),
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
    dom = 'lBfrtip', buttons = list(list(extend = "csv", text = "Download Page", filename = "Page_ProteinGeneName", exportOptions = list(modifier = list(page = "current")))),
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
