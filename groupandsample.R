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
  md <- data.table::as.data.table(MetaData_long())
  n_samples_total <- data.table::uniqueN(md$sampleid)
  
  # Hide logic: 
  # - Only 1 group
  # - Groups are as many as samples (unique for every sample)
  # - All groups are numeric (e.g., age, weight)
  meta_summary <- md[, .(
    groups = list(as.character(unique(group))),
    is_numeric = all(!is.na(suppressWarnings(as.numeric(unique(group)))))
  ), by = type]
  
  meta_summary[, hide := (
    is_numeric | 
      lengths(groups) == 1 | 
      lengths(groups) == n_samples_total
  )]
  
  init_meta <- setNames(meta_summary$groups, meta_summary$type)
  hide_vec  <- setNames(meta_summary$hide, meta_summary$type)
  
  state$hide_types    <- hide_vec
  state$visible_types <- names(hide_vec)[!hide_vec]
  
  state$meta    <- init_meta
  state$samples <- unique(md$sampleid)
  state$tests   <- if (!is.null(all_tests())) all_tests() else character()
  
  state$initial_meta    <- init_meta
  state$initial_samples <- state$samples
  state$initial_tests   <- state$tests
  
  # Also update the reactiveVals that the UI depends on
  all_group_list(init_meta)
  active_group_list(init_meta)
  all_samples(state$samples)
  sample_order(state$samples)
  all_tests(state$tests)
  test_order(state$tests)
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
reconcile_from_meta <- function(old_meta, ui_meta, ui_tests, ui_samples, md, all_samples_vec, all_tests_vec) {
  attrs <- names(ui_meta)
  
  # 1) Detect changed attributes
  changed_attrs <- attrs[!mapply(identical, lapply(ui_meta, as.character), lapply(old_meta, as.character))]
  
  if (length(changed_attrs) == 0) {
    return(list(meta = ui_meta, samples = ui_samples, tests = ui_tests))
  }
  
  # Convert md to data.table for high-speed indexing
  data.table::setDT(md) 
  base_order <- unique(c(ui_samples, all_samples_vec))
  new_samples <- ui_samples # Default fallback
  
  # 2) Recompute sample order driven by changed attributes
  for (attr in changed_attrs) {
    ui_vals <- ui_meta[[attr]]
    
    # Fast binary search subsetting
    blocks <- lapply(ui_vals, function(val) {
      relevant_samples <- md[type == attr & group == val, sampleid]
      base_order[base_order %in% relevant_samples]
    })
    
    new_samples <- unlist(blocks, use.names = FALSE)
  }
  
  # 3) Recompute meta for ALL attributes using reordered md slices
  new_meta <- lapply(attrs, function(attr) {
    md[type == attr][match(new_samples, sampleid), unique(as.character(group))]
  })
  names(new_meta) <- attrs

  # 4) Recompute tests from samples
  logical_tests <- all_tests_vec[vapply(all_tests_vec, function(t) {
    # Get pre-calculated IDs for this specific test
    required_samples <- test_id_lookup()[[t]]
    
    # If the test exists in our lookup, check if all its samples are in the current set
    if (!is.null(required_samples)) {
      return(all(required_samples %in% new_samples))
    } else {
      return(FALSE)
    }
  }, FUN.VALUE = logical(1))]
  
  # preserve UI order
  kept_tests  <- ui_tests[ui_tests %in% logical_tests]
  added_tests <- setdiff(logical_tests, kept_tests)
  new_tests   <- c(kept_tests, added_tests)
  
  list(meta = new_meta, samples = new_samples, tests = new_tests)
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
  md_dt <- data.table::as.data.table(md)
  active_md <- md_dt[sampleid %in% ui_samples]
  
  # 3) Recompute meta from samples
  attrs <- names(ui_meta)
  new_meta <- lapply(attrs, function(attr) {
    logical_groups <- unique(as.character(active_md[type == attr, group]))
    current_ui_attr <- ui_meta[[attr]]
    kept <- current_ui_attr[current_ui_attr %in% logical_groups]
    added <- data.table::setdiff(logical_groups, kept)
    return(c(kept, added))
  })
  names(new_meta) <- attrs  
  
  # 4) Recompute tests from samples
  logical_tests <- all_tests_vec[vapply(all_tests_vec, function(t) {
    # Get pre-calculated IDs for this specific test
    required_samples <- test_id_lookup()[[t]]
    
    # If the test exists in our lookup, check if all its samples are in the current set
    if (!is.null(required_samples)) {
      return(all(required_samples %in% samples))
    } else {
      return(FALSE)
    }
  }, FUN.VALUE = logical(1))]
  
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
    
    # Use data.table for a single-pass extraction
    md_dt <- data.table::as.data.table(md)
    data.table::setkey(md_dt, sampleid)
    md_filt <- md_dt[.(out$samples), nomatch = NULL]
    
    attrs <- names(old_meta)
    hidden_attrs <- setdiff(attrs, visible)
    
    if (length(hidden_attrs) > 0) {
      calc_dt <- md_filt[type %in% hidden_attrs]
      calc_meta <- split(calc_dt$group, calc_dt$type)
      calc_meta <- lapply(calc_meta, function(x) unique(as.character(x)))
      out$meta <- c(out$meta, calc_meta)
    }
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
  total_samples <- length(state$initial_samples)
  kept_samples  <- length(state$samples)
  total_tests   <- length(state$initial_tests)
  kept_tests    <- length(state$tests)
  summary_text <- paste0(
    "Selected ", kept_samples, " / ", total_samples, " Samples; ",
    kept_tests,   " / ", total_tests,   " Comparisons. ",
    "(Update Selection at: Top Menu → Groups and Samples.)"
  )
  tagList(
    tags$p(summary_text),
    tags$hr()
  )
})

output$selectGroupSample <- renderText({
  total_samples <- length(state$initial_samples)
  kept_samples  <- length(state$samples)
  total_tests   <- length(state$initial_tests)
  kept_tests    <- length(state$tests)
  paste0(
    "Selected ", kept_samples, " / ", total_samples, " Samples; ",
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
   if (!data.table::is.data.table(data_long)) {
    data.table::setDT(data_long)
  }
  
  data.table::setkey(data_long, sampleid)
  filtered <- data_long[.(sample_order), nomatch = NULL]
  setDT(filtered)
  cols_to_fix <- names(filtered)[vapply(filtered, function(x) is.character(x) || is.factor(x), logical(1))]
  for (col in cols_to_fix) {
    # Extract the column values once to avoid repeated indexing
    vals <- as.character(filtered[[col]])
    data.table::set(filtered, j = col, value = factor(vals, levels = unique(vals)))
  }  
  return(as.data.frame(filtered))
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
    
    system.time({
      tmp_data_long = filter_data_long(data_long, selected_groups, input_samples)
    })
    
    # 1. Ensure MetaData is a data.table
    md_dt <- data.table::as.data.table(MetaData)
    
    # 2. Filter and Arrange in one step using keyed subsetting
    # This filters for the samples and forces the order of input_samples
    data.table::setkey(md_dt, sampleid)
    md_dt <- md_dt[.(input_samples), nomatch = NULL]
    
    # 3. Optimize the column-wise factor conversion
    # Identify character or factor columns once
    cols_to_fix <- names(md_dt)[vapply(md_dt, function(x) is.character(x) || is.factor(x), logical(1))]
    
    # Update columns in-place using the 'set' function (zero-copy)
    for (col in cols_to_fix) {
      # Convert to character first, then to factor with unique levels
      data.table::set(md_dt, j = col, value = factor(as.character(md_dt[[col]]), 
                                                        levels = unique(as.character(md_dt[[col]]))))
    }
    
    tmp_sampleid = md_dt$sampleid
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
                "MetaData"= as.data.frame(md_dt), 
                'ProteinGeneName' = ProteinGeneName_filtered)
    )
  }
})