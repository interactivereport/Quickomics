####################################################
working_project=reactiveVal()
observe({
  req(ProjectInfo)
  working_project(ProjectInfo$ProjectID)
})

## =========================
## CORE STATE + LOCK
## =========================
active_group_list <- reactiveVal(list())
all_group_list <- reactiveVal(list())
updating_from_model <- reactiveVal(FALSE)
reset_all <- reactiveVal(FALSE)

state <- reactiveValues(
  meta    = list(),       # named list attr -> ordered groups
  samples = character(),  # ordered sample ids
  tests   = character(),   # ordered comparison ids
  initial_meta = list(),
  initial_samples = character(),
  initial_tests = character()
)

# 2. Initialize the lists when metadata changes
observeEvent(MetaData_long(), {
  # browser()
  md <- MetaData_long()
  # Set initial values
  init_meta <- split(md$group, md$type) %>% lapply(unique) %>% lapply(as.character)
  
  n_samples_total <- length(unique(md$sampleid))
  # compute hide vector ONCE
  hide <- vapply(init_meta, function(groups) {
    if (all(suppressWarnings(!is.na(as.numeric(groups))))) return(TRUE)
    if (length(groups) == 1) return(TRUE)
    if (length(groups) == n_samples_total) return(TRUE)
    FALSE
  }, logical(1))
  # store it
  state$hide_types <- hide
  state$visible_types <- names(state$hide_types)[!state$hide_types]
  
  init_samples <- unique(md$sampleid)
  init_tests <- if(!is.null(all_tests())) all_tests() else character()
  # Update the stable state
  state$meta <- init_meta
  state$samples <- init_samples
  state$tests <- init_tests
  # also store a frozen copy for reset
  state$initial_meta    <- init_meta
  state$initial_samples <- init_samples
  state$initial_tests   <- init_tests
  # Also update the reactiveVals that the UI depends on
  all_group_list(init_meta)
  active_group_list(init_meta)
  all_samples(init_samples)
  sample_order(init_samples)
  all_tests(init_tests)
  test_order(init_tests)
})


## =========================
## HELPERS
## =========================
read_ui_meta <- function(expected_attrs) {
  ui_meta <- lapply(expected_attrs, function(a) input[[paste0("keep_", a)]])
  names(ui_meta) <- expected_attrs
  ui_meta
}

detect_change_source <- function(visible, old_meta, old_samples, old_tests,
                                 ui_meta, ui_samples, ui_tests) {
  # visible <- state$visible_types
  meta_changed <- !identical(lapply(old_meta[visible], as.character),lapply(ui_meta[visible],  as.character))
  # meta_changed    <- !identical(lapply(old_meta, as.character),lapply(ui_meta,  as.character))
  samples_changed <- !identical(as.character(old_samples),as.character(ui_samples))
  tests_changed   <- !identical(as.character(old_tests),as.character(ui_tests))
  # last user action wins: tests > samples > meta
  if (tests_changed)   return("tests")
  if (samples_changed) return("samples")
  if (meta_changed)    return("meta")
  "none"
}


## =========================
## RECONCILIATION FUNCTIONS
## =========================
# A) META → SAMPLES → TESTS
reconcile_from_meta <- function(old_meta, ui_meta, ui_tests, ui_samples,
                                md, all_samples_vec, all_tests_vec) {
  attrs <- names(ui_meta)
  
  # 1) Detect changed attributes (including reorder-only)
  changed_attrs <- attrs[!mapply(
    identical,
    lapply(ui_meta,  as.character),
    lapply(old_meta, as.character)
  )]
  
  # If nothing changed, return as-is
  if (length(changed_attrs) == 0) {
    return(list(meta = ui_meta,
                samples = ui_samples,
                tests = ui_tests))
  }
  
  # 2) Recompute sample order driven by changed attributes
  #    Stable block reorder: UI order of groups → sample blocks
  base_order = unique(c(ui_samples, all_samples_vec))
  
  for (attr in changed_attrs) {
    # UI order for this attribute
    ui_vals <- ui_meta[[attr]]
    
    # For each value, get samples in that group
    blocks <- lapply(ui_vals, function(val) {
      md %>%
        dplyr::filter(type == attr, group == val) %>%
        dplyr::pull(sampleid) %>%
        intersect(base_order)   # preserve internal order
    })
    # Concatenate blocks → new sample order
    new_samples <- unlist(blocks, use.names = FALSE)
  }
  
  # 3) Recompute meta for ALL attributes using reordered md slices
  new_meta <- lapply(attrs, function(attr) {
    # Filter md to this attribute
    md_attr <- md %>% dplyr::filter(type == attr)
    # Reorder md_attr by new_samples (your fix)
    md_attr <- md_attr[match(new_samples, md_attr$sampleid), ]
    # Extract ordered unique groups
    md_attr$group %>%
      unique() %>%
      as.character()
  })
  names(new_meta) <- attrs
  
  # 4) Recompute tests from new samples
  logical_tests <- all_tests_vec[sapply(all_tests_vec, function(t) {
    groups <- unlist(strsplit(as.character(t), "vs"))
    samples_t <- md %>%
      dplyr::filter(group %in% groups) %>%
      dplyr::pull(sampleid) %>%
      as.character()
    all(samples_t %in% new_samples)
  })]
  
  kept_tests  <- ui_tests[ui_tests %in% logical_tests]
  added_tests <- setdiff(logical_tests, kept_tests)
  new_tests   <- c(kept_tests, added_tests)
  
  list(meta = new_meta,
       samples = new_samples,
       tests = new_tests)
}


# B) SAMPLES → META → TESTS
reconcile_from_samples <- function(old_samples, ui_samples, ui_meta, ui_tests,
                                   md, all_samples_vec, all_tests_vec) {
  # 1) Detect if samples actually changed
  samples_changed <- !identical(as.character(old_samples),
                                as.character(ui_samples))
  if (!samples_changed) {
    return(list(meta = ui_meta, samples = ui_samples, tests = ui_tests))
  }
  
  # 2) Samples are driven directly by ui_samples (user action)
  samples <- ui_samples
  # 3) Recompute meta for ALL attributes from samples
  attrs <- names(ui_meta)
  new_meta <- lapply(attrs, function(attr) {
    logical_groups <- md %>%
      dplyr::filter(type == attr, sampleid %in% samples) %>%
      dplyr::pull(group) %>%
      unique() %>%
      as.character()
    # preserve user order
    kept  <- ui_meta[[attr]][ui_meta[[attr]] %in% logical_groups]
    added <- setdiff(logical_groups, kept)
    c(kept, added)
  })
  names(new_meta) <- attrs
  
  # 4) Recompute tests from samples
  logical_tests <- all_tests_vec[sapply(all_tests_vec, function(t) {
    groups <- unlist(strsplit(as.character(t), "vs"))
    samples_t <- md %>%
      dplyr::filter(group %in% groups) %>%
      dplyr::pull(sampleid) %>%
      as.character()
    all(samples_t %in% samples)
  })]
  # preserve UI order
  kept_tests  <- ui_tests[ui_tests %in% logical_tests]
  added_tests <- setdiff(logical_tests, kept_tests)
  new_tests   <- c(kept_tests, added_tests)
  list(meta = new_meta, samples = samples, tests = new_tests)
}


# C) TESTS → SAMPLES → META
reconcile_from_tests <- function(old_tests, ui_tests, ui_meta, ui_samples,
                                 md, all_tests_vec, all_samples_vec) {
  # 1) Detect if tests actually changed
  tests_changed <- !identical(as.character(old_tests),
                              as.character(ui_tests))
  if (!tests_changed) {
    return(list(meta = ui_meta, samples = ui_samples, tests = ui_tests))
  }
  kept_tests    <- ui_tests
  removed_tests <- setdiff(all_tests_vec, kept_tests)
  
  samples_from_test <- function(test_name) {
    groups <- unlist(strsplit(as.character(test_name), "vs"))
    md %>%
      dplyr::filter(group %in% groups) %>%
      dplyr::pull(sampleid) %>%
      as.character()
  }
  
  # 2) Samples in kept and removed tests
  samples_in_kept    <- unique(unlist(lapply(kept_tests, samples_from_test)))
  samples_in_removed <- unique(unlist(lapply(removed_tests, samples_from_test)))
  
  # 3) Samples only in removed tests
  samples_only_removed <- setdiff(samples_in_removed, samples_in_kept)
  
  # 4) Final samples = all samples − samples_only_removed
  final_samples <- all_samples_vec[!(all_samples_vec %in% samples_only_removed)]
  
  # preserve UI order where possible
  final_samples <- unique(c(ui_samples, all_samples_vec))
  final_samples <- final_samples[final_samples %in% (all_samples_vec[!(all_samples_vec %in% samples_only_removed)])]
  
  # 5) Recompute meta from final samples
  attrs <- names(ui_meta)
  
  new_meta <- lapply(attrs, function(attr) {
    logical_groups <- md %>%
      dplyr::filter(type == attr, sampleid %in% final_samples) %>%
      dplyr::pull(group) %>%
      unique() %>%
      as.character()
    
    kept  <- ui_meta[[attr]][ui_meta[[attr]] %in% logical_groups]
    added <- setdiff(logical_groups, kept)
    c(kept, added)
  })
  names(new_meta) <- attrs
  
  list(meta = new_meta, samples = final_samples, tests = kept_tests)
}


# MASTER RECONCILER
reconcile_state <- function(visible, old_meta, old_samples, old_tests,
                            ui_meta, ui_samples, ui_tests,
                            md, all_samples_vec, all_tests_vec) {
  # 1) Detect which source changed
  source <- detect_change_source(visible, old_meta, old_samples, old_tests,ui_meta, ui_samples, ui_tests)
  
  # 2) No change → return old state
  if (source == "none") {
    return(list(meta = old_meta,samples = old_samples,tests   = old_tests))
  }
  
  # 3) META is the driver
  if (source == "meta") {
    out <- reconcile_from_meta(old_meta[visible], ui_meta[visible], ui_tests, ui_samples, md, all_samples_vec, all_tests_vec)
    attrs <- names(old_meta)
    md_filt <- md %>%
      dplyr::filter(sampleid %in% out$samples) %>%
      dplyr::arrange(match(sampleid, out$samples))
    
    calc_meta <- lapply(attrs, function(attr) {
      md_filt %>%
        dplyr::filter(type == attr) %>%
        dplyr::pull(group) %>%
        unique() %>%          # preserves first appearance
        as.character()
    })
    names(calc_meta) <- attrs
    out$meta <- c(out$meta, calc_meta[setdiff(names(calc_meta), visible)])
    return(out)
  }
  
  # 4) SAMPLES is the driver
  if (source == "samples") {
    return(reconcile_from_samples(old_samples, ui_samples, ui_meta, ui_tests, md, all_samples_vec, all_tests_vec))
  }
  
  # 5) TESTS is the driver
  if (source == "tests") {
    return(reconcile_from_tests(old_tests, ui_tests, ui_meta, ui_samples, md, all_tests_vec, all_samples_vec))
  }
  
  # 6) Fallback (should never hit)
  list(meta = old_meta, samples = old_samples, tests = old_tests)
}


## =========================
## MODEL → UI RENDERING
## =========================

observe({
  withProgress(message = "Updating UI...", value = 0, {
    incProgress(0.1)
    
    req(length(state$meta) > 0)
    
    # drive these reactives from state
    active_group_list(state$meta)
    sample_order(state$samples)
    test_order(state$tests)
    
    incProgress(0.2)
    
    # META UI
    output$ui_all_types <- renderUI({
      all_meta <- state$meta
      hide <- state$hide_types
      
      tagList(
        lapply(names(all_meta), function(attr) {
          kept_items    <- state$meta[[attr]]
          removed_items <- setdiff(all_group_list()[[attr]], kept_items)
          # hide row if needed
          row_style <- if (hide[[attr]]) "display:none;" else ""
          tags$div(style = row_style,
                   fluidRow(
                     column(1, tags$b(attr)),
                     column(7, shinyjqui::orderInput(
                       inputId = paste0("keep_", attr),
                       label = "Keep (Drag to Order):",
                       items = kept_items,
                       width = "100%",
                       item_class = "primary",
                       connect = paste0("remove_", attr)
                     )),
                     column(2, shinyjqui::orderInput(
                       inputId = paste0("remove_", attr),
                       label = "Remove:",
                       items = removed_items,
                       width = "100%",
                       placeholder = "Drag here to remove",
                       item_class = "info",
                       connect = paste0("keep_", attr)
                     ))
                   )
          )
        })
      )
    })
    
    incProgress(0.5)
    
    # SAMPLES UI
    samples_all     <- all_samples()
    ordered_samples <- state$samples
    samples_remove  <- samples_all[!(samples_all %in% ordered_samples)]
    
    output$ui_source_s <- renderUI({
      shinyjqui::orderInput('source_s', 'Available Samples:',items = ordered_samples,width = '100%', item_class = 'success', connect = 'dest_s')
    })
    
    output$ui_dest_s <- renderUI({
      shinyjqui::orderInput('dest_s', 'Drag to Remove Samples:',items = samples_remove,width = '100%', placeholder = 'Drag items here...',item_class = 'success', connect = 'source_s')
    })
    
    # TESTS UI
    comps_all     <- all_tests()
    ordered_comps <- state$tests
    comps_remove  <- comps_all[!(comps_all %in% ordered_comps)]
    
    output$ui_source_test <- renderUI({
      shinyjqui::orderInput('source_test', 'Available Comparisons',items = ordered_comps,width = '100%', item_class = 'success', connect = 'dest_test')
    })
    
    output$ui_dest_test <- renderUI({
      shinyjqui::orderInput('dest_test', 'Drag to Remove Comparisons:',items = comps_remove,width = '100%', placeholder = 'Drag items here...',item_class = 'success', connect = 'source_test')
    })
    
    incProgress(1)
  })
})


## =========================
## UI → MODEL RECONCILIATION
## =========================

observe({
  withProgress(message = "Updating the sample meta, sample list and comparison list ...", {
    # browser()
    req(length(all_group_list()) > 0)
    
    expected_attrs <- names(all_group_list())
    req(length(expected_attrs) > 0)
    
    md            <- MetaData_long()
    all_samples_vec <- all_samples()
    all_tests_vec <- all_tests()
    
    # 1) read UI
    ui_meta    <- read_ui_meta(expected_attrs)
    ui_samples <- input$source_s
    ui_tests   <- input$source_test
    
    if (any(sapply(ui_meta, is.null))) return()
    if (is.null(ui_samples) || is.null(ui_tests)) return()
    
    # 2) check if UI has caught up with model
    visible_types <- names(state$hide_types)[!state$hide_types]
    
    ui_matches_model <- (
      identical(
        lapply(ui_meta[visible_types], as.character),
        lapply(state$meta[visible_types], as.character)
      ) &&
        identical(as.character(ui_samples), as.character(state$samples)) &&
        identical(as.character(ui_tests),   as.character(state$tests))
    )
    incProgress(0.3)
    
    # 3) If we are updating the model, WAIT until UI matches
    if (updating_from_model() && !ui_matches_model) {
      return()   # freeze reconciliation until UI catches up
    }
    
    # 4) If UI now matches model, release the lock
    if (updating_from_model() && ui_matches_model) {
      updating_from_model(FALSE)
      return()
    }
    
    if (reset_all() && !ui_matches_model) {
      return()
    }
    reset_all(FALSE)
    
    incProgress(0.8)
    
    # 5) If not updating_from_model, proceed with reconciliation
    new_state <- reconcile_state(
      visible       = state$visible_types,
      old_meta      = state$meta,
      old_samples   = state$samples,
      old_tests     = state$tests,
      ui_meta       = ui_meta,
      ui_samples    = ui_samples,
      ui_tests      = ui_tests,
      md            = md,
      all_samples_vec = all_samples_vec,
      all_tests_vec = all_tests_vec
    )
    
    # 6) identity gate vs current model
    same_meta <- identical(
      lapply(new_state$meta, as.character),
      lapply(state$meta,     as.character)
    )
    same_samples <- identical(
      as.character(new_state$samples),
      as.character(state$samples)
    )
    same_tests <- identical(
      as.character(new_state$tests),
      as.character(state$tests)
    )
    
    if (same_meta && same_samples && same_tests) return()
    
    # 7) commit model update
    updating_from_model(TRUE)
    
    state$meta    <- new_state$meta
    state$samples <- new_state$samples
    state$tests   <- new_state$tests
  })
})

observeEvent(input$reset_all_types, {
  reset_all(TRUE)  
  state$meta    <- state$initial_meta
  state$samples <- state$initial_samples
  state$tests   <- state$initial_tests
})

shared_header_content <- reactive({
  req(state$meta, state$samples, state$tests)
  total_groups  <- sum(lengths(state$initial_meta))
  kept_groups   <- sum(lengths(state$meta))
  total_samples <- length(state$initial_samples)
  kept_samples  <- length(state$samples)
  total_tests   <- length(state$initial_tests)
  kept_tests    <- length(state$tests)
  summary_text <- paste0(
    "Selected ", kept_groups,  " / ", total_groups,  " Groups; ",
    kept_samples, " / ", total_samples, " Samples; ",
    kept_tests,   " / ", total_tests,   " Comparisons. ",
    "(Update Selection at: Top Menu → Groups and Samples.)"
  )
  tagList(
    tags$p(summary_text),
    tags$hr()
  )
})

output$selectGroupSample <- renderText({
  total_groups  <- sum(lengths(state$initial_meta))
  kept_groups   <- sum(lengths(state$meta))
  total_samples <- length(state$initial_samples)
  kept_samples  <- length(state$samples)
  total_tests   <- length(state$initial_tests)
  kept_tests    <- length(state$tests)
  paste0(
    "Selected ", kept_groups,  " / ", total_groups,  " Groups; ",
    kept_samples, " / ", total_samples, " Samples; ",
    kept_tests,   " / ", total_tests,   " Comparisons."
  )
})


output$summaryDetail <- renderPrint({
  #### Meta detail ####
  meta_detail <- vapply(names(state$meta), function(attr) {
    initial <- state$initial_meta[[attr]]
    current <- state$meta[[attr]]
    removed <- setdiff(initial, current)
    added   <- setdiff(current, initial)
    
    paste0(
      attr, ": ",
      length(current), "/", length(initial),
      if (length(removed) > 0)
        paste0(" | removed: ", paste(removed, collapse=", ")),
      if (length(added) > 0)
        paste0(" | added: ", paste(added, collapse=", "))
    )
  }, character(1))
  
  #### Comparison detail ####
  removed_tests <- setdiff(state$initial_tests, state$tests)
  added_tests   <- setdiff(state$tests, state$initial_tests)
  
  test_detail <- paste0(
    "Comparisons: ",
    length(state$tests), "/", length(state$initial_tests),
    if (length(removed_tests) > 0)
      paste0(" | removed: ", paste(removed_tests, collapse=", ")),
    if (length(added_tests) > 0)
      paste0(" | added: ", paste(added_tests, collapse=", "))
  )
  
  # Use cat to print cleanly to the verbatim box
  cat("--- Meta Detail ---\n")
  cat(meta_detail, sep = "\n")
  cat("\n--- Comparison Detail ---\n")
  cat(test_detail)
})

filter_data_long <- function(data_long, active_group_list, sample_order) {
  # 1. Filter by metadata attributes
  filtered <- data_long
  for (attr in names(active_group_list)) {
    filtered <- filtered[filtered[[attr]] %in% active_group_list[[attr]], ]
  }
  # 2. Filter by sample order
  filtered <- filtered[filtered$sampleid %in% sample_order, ]
  # 3. Preserve sample order
  filtered$sampleid <- factor(filtered$sampleid, levels = sample_order)
  filtered <- filtered[order(filtered$sampleid), ]
  filtered[] <- lapply(filtered[], function(col) {
    if (is.factor(col) || is.character(col)) {
      col <- as.character(col)    
      factor(col, levels = unique(col)) 
    } else {
      col                                    
    }
  })
  filtered
}

DataQCReactive <- reactive({
  DataIn = DataReactive()
  if (is.null(DataIn$groups)) {
    DataIn = DataReactive()
    results_long = DataIn$results_long
    ProteinGeneName = DataIn$ProteinGeneName
    return(list('tmp_data_wide'=NULL,
                'tmp_data_long'=NULL,  
                'tmp_results_long' = results_long, 
                'tmp_group' = NULL, 
                'tmp_sampleid'=NULL, 
                "MetaData"=NULL, 
                'ProteinGeneName' = ProteinGeneName)
    )
  } else {
    req(length(active_group_list()) > 0 )
    DataIn = DataReactive()
    results_long = DataIn$results_long
    ProteinGeneName = DataIn$ProteinGeneName
    MetaData = DataIn$MetaData
    data_long = DataIn$data_long
    data_wide = DataIn$data_wide
    selected_groups <- active_group_list()
    all_groups <- all_group_list()
    
    input_samples = sample_order() # input$QC_samples
    
    tmp_data_long = filter_data_long(data_long, selected_groups, input_samples)
    
    MetaData <- MetaData %>%
      dplyr::filter(sampleid %in% input_samples) %>%
      dplyr::mutate(sampleid = factor(sampleid, levels = input_samples)) %>%
      dplyr::arrange(sampleid)
    
    MetaData[] <- lapply(MetaData[], function(col) {
      if (is.factor(col) || is.character(col)) {
        col <- as.character(col)    
        factor(col, levels = unique(col)) 
      } else {
        col                                    
      }
    })
    
    tmp_sampleid = MetaData$sampleid
    tmp_data_wide = data_wide[, as.character(tmp_sampleid), drop = FALSE] %>% as.matrix()
    input_tests <- test_order()
    tmp_results_long <- results_long %>%
      dplyr::filter(test %in% input_tests) %>%
      dplyr::mutate(test = factor(test, levels = input_tests)) %>%
      dplyr::arrange(test)
    ProteinGeneName_filtered <- ProteinGeneName[ProteinGeneName$UniqueID %in% rownames(data_wide),]
    
    return(list('tmp_data_wide'=tmp_data_wide,
                'tmp_data_long'=tmp_data_long,  
                'tmp_results_long' = tmp_results_long, 
                'tmp_group' = selected_groups, 
                'tmp_sampleid'=tmp_sampleid, 
                "MetaData"=MetaData, 
                'ProteinGeneName' = ProteinGeneName_filtered)
    )
  }
})