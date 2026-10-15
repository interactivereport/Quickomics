 ###########################################################################################################
##This software belongs to Biogen Inc. All right reserved.
##
##@file: TimeSeriesmodule.R
##@Developer : Wei Li (wei.li@Biogen.com)
##@Date : 11/06/2025
##@version 1.0
###########################################################################################################
## correlation analysis
##########################################################################################################
#pkgs:  

library(DESeq2)
library(ggplot2)
library(DEGreport)

get_count_mtx <- function(data_long, linear_base, linear_small_value) {
  if (is.na(linear_base)) linear_base = 2
  if (is.na(linear_small_value)) linear_small_value = 1
  data_long_linear <- data_long %>% mutate(expr=pmax(linear_base^expr-linear_small_value, 0))
  
  data_wide_linear <- data_long_linear %>%
    dplyr::select(UniqueID, sampleid, expr) %>%
    tidyr::pivot_wider(
      names_from = sampleid,
      values_from = expr
    )
  data_wide_linear <- as.data.frame(data_wide_linear)
  rownames(data_wide_linear) <- data_wide_linear$UniqueID
  data_wide_linear$UniqueID <- NULL
  data_wide_linear[is.na(data_wide_linear)] <- 0
  
  data_wide_linear[] <- lapply(data_wide_linear, function(col) {
    if (is.numeric(col)) {
      as.integer(round(col))
    } else {
      col
    }
  })
  return(data_wide_linear)
}

parse_time <- function(x) {
  # Extract numeric part
  value <- as.numeric(sub("([0-9.]+).*", "\\1", x))
  
  # Extract unit part
  unit <- tolower(sub("[0-9.]+\\s*(.*)", "\\1", x))
  unit <- trimws(unit)
  
  # Unit conversion map (to hours)
  unit_map <- c(
    # seconds
    "s" = 1/3600, "sec" = 1/3600, "second" = 1/3600, "seconds" = 1/3600,
    # minutes
    "min" = 1/60, "mins" = 1/60, "minute" = 1/60, "minutes" = 1/60,
    # hours
    "h" = 1, "hr" = 1, "hrs" = 1, "hour" = 1, "hours" = 1,
    # days
    "d" = 24, "day" = 24, "days" = 24,
    # weeks
    "wk" = 24*7, "week" = 24*7, "weeks" = 24*7,
    # months
    "m" = 24*30, "mo" = 24*30, "mon" = 24*30, "month" = 24*30, "months" = 24*30,
    # years
    "yr" = 24*365, "year" = 24*365, "years" = 24*365
  )
  
  # Match units
  multiplier <- unit_map[unit]
  if (any(is.na(multiplier))) {
    stop("Unknown time units: ", paste(unique(unit[is.na(multiplier)]), collapse = ", "))
  }
  
  # Convert to hours
  hours <- value * multiplier
  
  # Canonical unit mapping
  canonical_unit <- c(
    "s"="s","sec"="s","second"="s","seconds"="s",
    "m"="min","min"="min","mins"="min","minute"="min","minutes"="min",
    "h"="h","hr"="h","hrs"="h","hour"="h","hours"="h",
    "d"="d","day"="d","days"="d",
    "wk"="wk","week"="wk","weeks"="wk",
    "mo"="mo","mon"="mo","month"="mo","months"="mo",
    "yr"="yr","year"="yr","years"="yr"
  )
  
  # Determine majority unit
  majority_unit <- names(sort(table(canonical_unit[unit]), decreasing = TRUE))[1]
  
  # Conversion factors (hours → target unit)
  unit_to_hours <- c(
    "s" = 1/3600,
    "min" = 1/60,
    "h" = 1,
    "d" = 24,
    "wk" = 24*7,
    "mo" = 24*30,
    "yr" = 24*365
  )
  
  # Convert hours → majority unit
  converted <- hours / unit_to_hours[majority_unit]
  
  # ---- NEW: format to max 3 decimals, no trailing zeros ----
  converted_fmt <- sub("\\.?0+$", "", format(round(converted, 3), nsmall = 3))
  
  # Build canonical labels
  canonical <- paste0(converted_fmt, majority_unit)
  
  # Order levels by numeric time
  ord <- order(hours)
  
  # Return UNORDERED factor with ordered levels
  factor(canonical, levels = unique(canonical[ord]), ordered = FALSE)
}

# Sample table and LRT model formulas shared by the DESeq2 and limma tests.
# test = "interaction": genes whose time trend differs between levels of the
#   last Condition variable (full ~ Condition + Time + lastCondition:Time,
#   reduced ~ Condition + Time).
# test = "time": genes that change over time, adjusting for Condition
#   (full ~ Condition + Time, reduced ~ Condition; with no Condition,
#   ~ Time vs ~ 1).
ts_prepare_design <- function(MetaData, samples, time_variable, Condition, test = "interaction") {
  # -----------------------------
  # 1. Prepare metadata
  # -----------------------------
  sample_meta <- MetaData %>%
    filter(sampleid %in% samples) %>%
    arrange(match(sampleid, samples))
  # Convert all factor columns to character (DEGreport safety)
  sample_meta[] <- lapply(sample_meta, function(x) {
    if (is.factor(x)) as.character(x) else x
  })
  # Create ordered time factor
  if (!is.numeric(sample_meta[[time_variable]])) {
    sample_meta$Unified_Time <- parse_time(sample_meta[[time_variable]])
  } else {
    sample_meta$Unified_Time <- factor(sample_meta[[time_variable]], levels = sort(unique(sample_meta[[time_variable]])))
  }
  # -----------------------------
  # 1b. Factor all Condition columns
  # -----------------------------
  Condition <- setdiff(Condition, time_variable)
  for (col in Condition) {
    sample_meta[[col]] <- factor(sample_meta[[col]])
  }
  if (nlevels(droplevels(sample_meta$Unified_Time)) < 2) {
    return(list(status = "error", msg = "The time variable needs at least 2 time points in the selected samples."))
  }
  if (test == "interaction" && length(Condition) == 0) {
    return(list(status = "error", msg = "Time x Condition test needs at least one Condition Variable. Use the Time only test to find genes that change over time without a condition."))
  }
  # Constant "condition" so degPatterns and the plots/heatmaps downstream
  # (which all expect one condition column) still work with no condition.
  if (length(Condition) == 0) sample_meta[[ts_no_condition_col]] <- "All"

  # Set rownames
  rownames(sample_meta) <- sample_meta$sampleid

  # -----------------------------
  # 2. Build model formulas
  # -----------------------------
  if (test == "time") {
    # Main effects (if any) stay in both models; Time is what gets tested.
    main_effects <- if (length(Condition) > 0) paste(Condition, collapse = " + ") else "1"
    full_formula    <- as.formula(paste("~", if (length(Condition) > 0) paste(main_effects, "+ Unified_Time") else "Unified_Time"))
    reduced_formula <- as.formula(paste("~", main_effects))
  } else {
    # Last covariate is the one to test interaction
    interaction_cov <- tail(Condition, 1)
    # Main effects: all covariates
    main_effects <- paste(Condition, collapse = " + ")

    # Full model: all covariates + time + lastCov:time
    full_formula <- as.formula(
      paste("~", main_effects, "+ Unified_Time +",
            paste0(interaction_cov, ":Unified_Time"))
    )

    # Reduced model: drop interaction
    reduced_formula <- as.formula(
      paste("~", main_effects, "+ Unified_Time")
    )
  }
  list(status = "success", sample_meta = sample_meta, Condition = Condition,
       full_formula = full_formula, reduced_formula = reduced_formula)
}

# Detailed message shown when the full model cannot be fitted. detail adds
# method-specific information (e.g. which limma coefficients are not estimable).
ts_design_error_msg <- function(method, d, detail = NULL) {
  method_name <- if (identical(method, "limma")) "limma" else "DESeq2"
  causes <- c(
    "- Some combinations of your variables (Condition x Time) have zero samples.",
    "- One variable (e.g., 'Batch') perfectly overlaps with another (e.g., 'Condition') or with Time.",
    "- A column in your metadata contains only one unique value in the selected samples."
  )
  if (identical(method, "limma")) causes <- c(causes,
    "- There are not more samples than model coefficients, so limma has no residual degrees of freedom to estimate the variance.")
  paste(c(
    sprintf("Design Matrix Error (%s): The variables in your formula are confounded, so the model cannot be fitted.", method_name),
    sprintf("Full model: %s   Reduced model: %s", paste(deparse(d$full_formula), collapse = ""), paste(deparse(d$reduced_formula), collapse = "")),
    detail,
    "Possible causes:",
    causes,
    "Please simplify your design (remove a Condition Variable, or use the Time only test) or check your sample metadata for redundant columns."
  ), collapse = "\n")
}

# DESeq2 LRT on counts (RNA-seq). Returns the common time-course result:
# results (one row per gene, with pvalue/padj), expr (VST matrix used by
# clustering and Similar Pattern Search), sample_meta, plus dds.
run_timecourse_DEG <- function(MetaData, count_mtx,time_variable,Condition, test = "interaction") {
  library(DESeq2)
  library(dplyr)
  library(tibble)

  d <- ts_prepare_design(MetaData, colnames(count_mtx), time_variable, Condition, test)
  if (d$status != "success") return(d)
  sample_meta <- d$sample_meta
  count_mtx <- count_mtx[, sample_meta$sampleid]

  # -----------------------------
  # 3. Build DESeq2 object
  # -----------------------------
  dds <- tryCatch({
    DESeqDataSetFromMatrix(
      countData = count_mtx,
      colData   = sample_meta,
      design    = d$full_formula
    )
  }, error = function(e) e)

  # Errors here come from the design (not full rank, a variable with one
  # value, ...); show them with the detailed design message.
  if (inherits(dds, "error")) {
    detail <- if (!grepl("full rank", conditionMessage(dds)))
      paste("Model error:", gsub("\\s+", " ", conditionMessage(dds)))
    return(list(status = "error", msg = ts_design_error_msg("deseq2", d, detail)))
  }
    
  # Filter low-count genes
  dds <- dds[rowSums(counts(dds)) > 1, ]
  # LRT test
  dds <- DESeq(dds, test = "LRT", reduced = d$reduced_formula)

  res <- as.data.frame(results(dds))
  return(list(status = "success", method = "deseq2", dds = dds, test = test,
              results = res, mean_col = "baseMean",
              expr = assay(vst(dds, blind = TRUE)), expr_label = "VST",
              sample_meta = as.data.frame(colData(dds)),
              full_formula = deparse(d$full_formula), reduced_formula = deparse(d$reduced_formula)))
}

# limma moderated F-test of the same full vs reduced models, for continuous
# log-scale data (proteomics, microarray, metabolomics). expr_mtx is genes x
# samples on the log scale, NA = missing. Genes are kept only with at least
# min_per_group non-missing values in every condition x time group, so every
# group mean exists for clustering and Similar Pattern Search; remaining NAs
# are left to limma (it fits each gene on the samples it has).
run_timecourse_limma <- function(MetaData, expr_mtx, time_variable, Condition, test = "interaction", min_per_group = 2) {
  library(limma)

  d <- ts_prepare_design(MetaData, colnames(expr_mtx), time_variable, Condition, test)
  if (d$status != "success") return(d)
  sample_meta <- d$sample_meta
  expr_mtx <- as.matrix(expr_mtx[, sample_meta$sampleid, drop = FALSE])

  cond_var <- if (length(d$Condition) > 0) tail(d$Condition, 1) else ts_no_condition_col
  grp <- paste(sample_meta[[cond_var]], sample_meta$Unified_Time, sep = "_")
  grp_idx <- split(seq_len(ncol(expr_mtx)), grp)
  n_ok <- sapply(grp_idx, function(ix) rowSums(!is.na(expr_mtx[, ix, drop = FALSE])))
  if (is.null(dim(n_ok))) n_ok <- matrix(n_ok, nrow = nrow(expr_mtx))
  # A group with fewer samples than min_per_group only needs all its samples.
  need <- pmin(min_per_group, lengths(grp_idx))
  keep <- apply(sweep(n_ok, 2, need, ">="), 1, all)
  n_total <- nrow(expr_mtx)
  expr_mtx <- expr_mtx[keep, , drop = FALSE]
  if (nrow(expr_mtx) < 2) {
    return(list(status = "error", msg = sprintf("Only %d genes/proteins have at least %d values in every condition x time group.", nrow(expr_mtx), min_per_group)))
  }

  design_full <- tryCatch(model.matrix(d$full_formula, data = sample_meta), error = function(e) e)
  if (inherits(design_full, "error")) {
    return(list(status = "error", msg = ts_design_error_msg("limma", d, paste("Model error:", gsub("\\s+", " ", conditionMessage(design_full))))))
  }
  design_red  <- model.matrix(d$reduced_formula, data = sample_meta)
  not_est <- limma::nonEstimable(design_full)
  if (qr(design_full)$rank < ncol(design_full) || length(not_est) > 0) {
    return(list(status = "error", msg = ts_design_error_msg("limma", d,
      if (length(not_est) > 0) paste("Coefficients that cannot be estimated:", paste(not_est, collapse = ", ")))))
  }
  if (nrow(design_full) <= ncol(design_full)) {
    return(list(status = "error", msg = ts_design_error_msg("limma", d,
      sprintf("The model has %d coefficients but only %d samples.", ncol(design_full), nrow(design_full)))))
  }
  test_coefs <- setdiff(colnames(design_full), colnames(design_red))

  fit <- lmFit(expr_mtx, design_full)
  fit <- eBayes(fit, trend = TRUE)
  tt <- topTable(fit, coef = test_coefs, number = Inf, sort.by = "none")
  # With a single tested coefficient (e.g. 2 time points) topTable reports a
  # moderated t instead of F; F = t^2 and the p-values are the same test.
  F_stat <- if (!is.null(tt$F)) tt$F else tt$t^2
  res <- data.frame(AveExpr = tt$AveExpr, F = F_stat, pvalue = tt$P.Value, padj = tt$adj.P.Val,
                    row.names = rownames(tt))

  # For clustering/similarity only: fill each NA with its own condition x
  # time group mean, which leaves every group mean unchanged.
  expr_filled <- expr_mtx
  for (g in unique(grp)) {
    ix <- which(grp == g)
    sub <- expr_filled[, ix, drop = FALSE]
    m <- rowMeans(sub, na.rm = TRUE)
    na_pos <- which(is.na(sub), arr.ind = TRUE)
    if (nrow(na_pos) > 0) sub[na_pos] <- m[na_pos[, 1]]
    expr_filled[, ix] <- sub
  }

  return(list(status = "success", method = "limma", test = test,
              results = res, mean_col = "AveExpr",
              expr = expr_filled, expr_label = "log expression",
              sample_meta = sample_meta, n_filtered = n_total - nrow(expr_mtx), min_per_group = min_per_group,
              full_formula = deparse(d$full_formula), reduced_formula = deparse(d$reduced_formula)))
}

run_DEG_cluster <- function(expr_mat, sample_meta, sig_gene, condition_var) {
  cluster_rlog <- expr_mat[sig_gene, ]
  
  clusters <- degPatterns(
    cluster_rlog,
    metadata = sample_meta,
    time     = "Unified_Time",
    col      = condition_var,
    plot     = FALSE
  )
  
  return(list(
    rlog = expr_mat,
    cluster_rlog = cluster_rlog,
    sample_meta = sample_meta,
    clusters = clusters
  ))
}

ts_no_condition_col <- "AllSamples"

# Positive whole number from a numericInput, falling back to a default
# while the box is empty or invalid.
ts_pos_int <- function(x, default) {
  x <- suppressWarnings(as.integer(x))
  if (length(x) != 1 || is.na(x) || x < 1) default else x
}

# Mean VST expression per (condition, time) group, the same collapsing
# degPatterns does before clustering. Returns genes x groups, with groups
# ordered by condition then time, plus the group table.
ts_group_profiles <- function(vst_mat, sample_meta, cond_var) {
  sample_meta <- as.data.frame(sample_meta)
  cond <- as.character(sample_meta[[cond_var]])
  cond_levels <- if (is.factor(sample_meta[[cond_var]])) levels(droplevels(sample_meta[[cond_var]])) else unique(cond)
  time <- droplevels(as.factor(sample_meta$Unified_Time))
  groups <- expand.grid(time = levels(time), condition = cond_levels, stringsAsFactors = FALSE)
  groups$label <- paste(groups$condition, groups$time, sep = "_")
  key <- paste(cond, as.character(time), sep = "_")
  groups <- groups[groups$label %in% key, , drop = FALSE]
  prof <- sapply(groups$label, function(g) rowMeans(vst_mat[, key == g, drop = FALSE]))
  if (is.null(dim(prof))) prof <- matrix(prof, nrow = nrow(vst_mat), dimnames = list(rownames(vst_mat), groups$label))
  groups$time <- factor(groups$time, levels = levels(time))
  list(profiles = prof, groups = groups)
}

# Rank genes by how closely their group-mean profile follows the anchor
# gene's. Kendall is what degPatterns clusters on; with few time points
# Kendall has many ties, so Pearson on the same profile breaks them, then
# change_distance (RMS difference of the mean-centred profiles, i.e. how
# close the size of each change is to the anchor's). Both r values are
# compared at 2 decimals. With only 2 groups every correlation is +1 or -1,
# so change_distance does the ranking.
find_similar_genes <- function(profiles, anchor, method = "kendall", direction = "similar", candidates = NULL) {
  stopifnot(anchor %in% rownames(profiles))
  a <- profiles[anchor, ]
  pool <- if (is.null(candidates)) rownames(profiles) else intersect(candidates, rownames(profiles))
  pool <- setdiff(pool, anchor)
  m <- profiles[pool, , drop = FALSE]
  m <- m[apply(m, 1, function(x) sd(x) > 0), , drop = FALSE]
  r_main <- if (method == "pearson") cor(t(m), a) else cor(t(m), a, method = method)
  r_pear <- cor(t(m), a)
  sgn <- if (direction == "opposite") -1 else 1
  a_c <- sgn * (a - mean(a))
  m_c <- m - rowMeans(m)
  change_dist <- sqrt(rowMeans(sweep(m_c, 2, a_c)^2))
  df <- data.frame(UniqueID = rownames(m), r = as.vector(r_main), pearson_r = as.vector(r_pear),
                   change_distance = as.vector(change_dist), stringsAsFactors = FALSE)
  # r values that agree to 2 decimals count as tied, so among genes with the
  # same shape the one whose changes are closest in size comes first.
  # Comparing r at full precision let 7th-decimal differences decide, which
  # put nearly flat genes ahead of genes changing as much as the anchor.
  df <- df[order(-sgn * round(df$r, 2), -sgn * round(df$pearson_r, 2), df$change_distance), ]
  df$rank <- seq_len(nrow(df))
  rownames(df) <- NULL
  df
}

# ---- UI Module ----
TimeSeries_ui <- function(id) {
  ns <- shiny::NS(id)
  fluidRow(
    column(3,
           wellPanel(
             column(width=12,uiOutput(ns("selectGroupSampleTimeSeries"))),
             radioButtons(ns("ts_genelable"),label="Select Gene Label",inline = TRUE, choices=c("Gene.Name","UniqueID"), selected="Gene.Name"),
             conditionalPanel(ns = ns, "input.TimeSeries_tabset =='Time-Series DE Analysis'",
                              selectInput(ns("sel_time_var"), label="Select Time Variable", choices=NULL),
                              selectizeInput(ns("sel_condition_var"), label="Select Condition Variable (optional for Time only)", choices=NULL,multiple=TRUE),
                              radioButtons(ns("ts_method"), label="Data Type / Method",
                                           choices = c("RNA-seq counts (DESeq2 LRT)" = "deseq2",
                                                       "Continuous log values, e.g. proteomics (limma F-test)" = "limma"),
                                           selected = "deseq2"),
                              conditionalPanel(ns = ns, "input.ts_method == 'limma'",
                                               numericInput(ns("ts_min_per_group"), label = "Min. Non-missing Values per Condition x Time Group",
                                                            value = 2, min = 1, step = 1)),
                              radioButtons(ns("ts_test"), label="Test",
                                           choices = c("Time x Condition (trend differs between conditions)" = "interaction",
                                                       "Time only (changes over time)" = "time"),
                                           selected = "interaction"),
                              actionButton(ns("compute_DE"), "Run Time-Series DE Analysis", style="color: #0961E3; background-color: #F6E98C ; border-color: #2e6da4"),
                              tags$hr()
             ),
             conditionalPanel(ns = ns, "input.TimeSeries_tabset =='Time-Series DE Analysis' || input.TimeSeries_tabset =='DEG Clustering Analysis' || input.TimeSeries_tabset =='Similar Pattern Search'",
                              sliderInput(ns("padjCutoff"), label="Set the Adjusted p-value Cutoff:", min = 0, max = 1, step = 0.01, value = 0.05),
                              sliderInput(ns("top_n"), label="Select Number of Top DEGs for Clustering", min = 300, max = 1500, step = 100, value = 1000),
                              span(textOutput(ns("filtered_DEG")), style = "color:blue; font-size:15px; font-family:arial; font-style:italic")
             ),
             conditionalPanel(ns = ns, "input.TimeSeries_tabset =='Similar Pattern Search'",
                              tags$hr(),
                              selectizeInput(ns("sim_gene"), label="Gene to Match", choices=NULL, multiple=FALSE,
                                             options = list(placeholder = "Type to search")),
                              radioButtons(ns("sim_pool"), label="Search Among", inline = TRUE,
                                           choices = c("All tested genes" = "all", "Time-series DEGs above" = "deg"), selected = "all"),
                              radioButtons(ns("sim_direction"), label="Find", inline = TRUE,
                                           choices = c("Similar pattern" = "similar", "Opposite pattern" = "opposite"), selected = "similar"),
                              radioButtons(ns("sim_method"), label="Similarity", inline = TRUE,
                                           choices = c("Kendall" = "kendall", "Pearson" = "pearson", "Spearman" = "spearman"), selected = "kendall"),
                              sliderInput(ns("sim_top_n"), label="Number of Similar Genes to List", min = 10, max = 500, step = 10, value = 50),
                              sliderInput(ns("sim_plot_n"), label="Number of Top Similar Genes to Plot", min = 1, max = 100, step = 1, value = 10),
                              radioButtons(ns("sim_plot_scale"), label="Plot Scale", inline = TRUE,
                                           choices = c("Centered log2 (keeps change size)" = "centered", "Z-score" = "zscore"), selected = "centered"),
                              actionButton(ns("compute_similar"), "Find Similar Genes", style="color: #0961E3; background-color: #F6E98C ; border-color: #2e6da4")
             ),
             conditionalPanel(ns = ns, "input.TimeSeries_tabset =='DEG Clustering Analysis'",
                              actionButton(ns("compute_cluster"), "Run Clustering", style="color: #0961E3; background-color: #F6E98C ; border-color: #2e6da4"),
                              conditionalPanel(ns = ns, "input.ts_cluster_tabset =='Time Series Cluster Plot'",
                                               tags$hr(),
                                               tags$p(
                                                 "Customize Plot:",
                                                 style = "font-size: 20px; font-weight: bold;"
                                               ),
                                               sliderInput(ns("plot_width"), label="Plot Width", min = 600, max = 2000, step = 100, value = 1200),
                                               sliderInput(ns("plot_height"), label="Plot Height", min = 600, max = 2000, step = 100, value = 800),
                                               sliderInput(ns("strip_text"), "Facet Title Size", min = 12, max = 25, value = 18),
                                               column(width=6,sliderInput(ns("axis_title"), "Axis Title Size", min = 12, max = 25, value = 20)),
                                               column(width=6,sliderInput(ns("axis_text"), "Axis Text Size", min = 12, max = 25, value = 18)),
                                               column(width=6,sliderInput(ns("legend_title"), "Legend Title Size", min = 12, max = 25, value = 20)),
                                               column(width=6,sliderInput(ns("legend_text"), "Legend Text Size", min = 12, max = 25, value = 18))
                                               ),
                              conditionalPanel(ns = ns, "input.ts_cluster_tabset =='Time‑series Sample Heatmap'",
                                               selectizeInput(ns("ts_hm_annot_sample"), label="Annotate Samples", choices=NULL, multiple = TRUE)
                                               ),
                              conditionalPanel(ns = ns, "input.ts_cluster_tabset =='Time‑series Group Heatmap'",
                                               selectizeInput(ns("ts_hm_annot_group"), label="Annotate Groups", choices=NULL, multiple = TRUE)
                                               ),
                              conditionalPanel(ns = ns, "input.ts_cluster_tabset =='Time‑series Sample Heatmap' || input.ts_cluster_tabset =='Time‑series Group Heatmap'",
                                               radioButtons(ns("ts_hm_gene_cluster"), label="Cluster genes", inline = TRUE, choices = c("Yes","No"), selected = "No"),
                                               column(width=5,sliderInput(ns("ts_hm_x_fontsize"), "Column Font Size:", min = 0, max = 20, step = 1, value = 15)),
                                               column(width=5,sliderInput(ns("ts_hm_y_fontsize"), "Row Font Size:", min = 0, max = 20, step = 1, value = 10)),
                                               sliderInput(ns("ts_hm_N_genes"), "Max Number of Genes to Label:", min = 0, max = 500, step = 10, value = 100),
                                               h5("After changing parameters, please click Plot/Refresh button in the plot panel to generate heatmap."),
                                               radioButtons(ns("ts_hm_more_options"), label="Show More Options", inline = TRUE, choices = c("Yes","No"), selected = "No"),
                                               conditionalPanel(ns = ns, "input.ts_hm_more_options=='Yes'",
                                                                radioButtons(ns("ts_hm_annot_color"),label="Color Setting for Annotations", choices=c("Auto-Set by Rand. Seed","Select Palette"), selected="Auto-Set by Rand. Seed"),
                                                                conditionalPanel(ns=ns, "input.ts_hm_annot_color=='Auto-Set by Rand. Seed'",
                                                                                 numericInput(ns("ts_hm_seed"),label="Random Seed for Color Palettes for Categories", min=1, max= 5000, value=123, step=1)),
                                                                conditionalPanel(ns=ns, "input.ts_hm_annot_color=='Select Palette'",
                                                                                 selectizeInput(ns("ts_hm_cat_pal"), label="Color Palettes (One per Category Annotation)",
                                                                                                choices=c("Dark2", "Accent",  "Set1", "Set2", "Set3", "npg", "nejm", "lancet", "jama", "d3", "uchicago"), multiple = TRUE),
                                                                                 selectizeInput(ns("ts_hm_num_pal"), label="Set Max Colors for Numerical Annotations",
                                                                                                choices=c("Dark2", "Accent",  "Set1", "Set2", "Set3", "npg", "nejm", "lancet", "jama", "d3", "uchicago"), selected="Set1", multiple = FALSE
                                                                                 )
                                                                ),
                                                                tags$hr(),
                                                                sliderInput(ns("ts_hm_height"), "Heatmap Height:", min = 200, max = 3000, step = 50, value = 800),
                                                                fluidRow(
                                                                  column(width=4,colourInput(ns("ts_hm_lowColor"), "Low", "blue")),
                                                                  column(width=4,colourInput(ns("ts_hm_midColor"), "Mid", "white")),
                                                                  column(width=4,colourInput(ns("ts_hm_highColor"), "High", "red"))
                                                                ),
                                                                fluidRow(
                                                                  textInput(ns("ts_hm_row_title"), "Row Title", width = "100%"),
                                                                  sliderInput(ns("ts_hm_row_title_font_size"), "Row Title Font Size:", min = 0, max = 30, step = 1, value = 16),
                                                                  textInput(ns("ts_hm_column_title"), "Column Title", width = "100%"),
                                                                  sliderInput(ns("ts_hm_column_title_font_size"), "Column Title Font Size:", min = 0, max = 30, step = 1, value = 16)
                                                                ),
                                               )
                                               )
             )
           )
    ),
    column(9,
           tabsetPanel(id = ns("TimeSeries_tabset"),
                       tabPanel(title="Time-Series DE Analysis",
                                br(),
                                tags$details(
                                  tags$summary(
                                    tags$strong("Click to view Time-Series DE Methodology (DESeq2 LRT / limma) & Example", 
                                                style = "color: #28a745; cursor: pointer;")
                                  ),
                                  tags$div(
                                    style = "padding: 15px; background-color: #f4faf6; border-left: 5px solid #28a745; margin-top: 10px;",
                                    tags$p("Both methods compare a full model to a reduced model. A gene is significant when the term(s) dropped from the reduced model explain a significant amount of its variation over time. ",
                                           tags$strong("The p-values describe the gene's whole trend, not a single fold change.")),
                                    tags$p(tags$strong("1. Data Type / Method (option on the left):")),
                                    tags$ul(
                                      tags$li(tags$strong("RNA-seq counts (DESeq2 LRT):"), " the project's log values (e.g. log2 TPM) are converted back to the linear scale (e.g. TPM), rounded to whole numbers and treated as counts, since raw counts are not available in the project; DESeq2's likelihood ratio test (LRT) then compares the two models. For RNA-seq data; because the input is TPM rather than real counts, p-values are approximate. Columns: baseMean (mean of the normalized values), stat (LRT chi-square), pvalue, padj."),
                                      tags$li(tags$strong("Continuous log values (limma F-test):"), " the same two models are compared with limma's moderated F-test (eBayes with intensity trend) on the project's log values. For proteomics, microarray, metabolomics or other normalized log-scale data; it can also be used on RNA-seq log values (limma-trend). Columns: AveExpr (mean log expression), F, pvalue, padj. ",
                                              "Missing values are kept as missing; genes/proteins with fewer than the set number of values in any condition x time group (or missing in any sample of a smaller group) are left out, so every group mean exists for clustering and Similar Pattern Search.")
                                    ),
                                    tags$p("The default method is chosen from the expression unit (TPM/FPKM/CPM/counts: DESeq2, otherwise limma) and can be changed."),
                                    tags$p(tags$strong("2. Test (option on the left), available with both methods:")),
                                    tags$ul(
                                      tags$li(tags$strong("Time x Condition:"), " finds genes whose time trend is different between the levels of the last Condition Variable. Needs at least one Condition Variable."),
                                      tags$li(tags$strong("Time only:"), " finds genes that change over time, adjusting for the Condition Variables. Condition Variables are optional.")
                                    ),
                                    tags$p(tags$strong("Example:"), " Time Variable ", tags$strong("Age"), ", Condition Variables ", tags$strong("Tissue, Treatment"), "."),
                                    tags$ul(
                                      tags$li("Time x Condition: full model ", tags$strong("Tissue + Treatment + Age + Treatment:Age"), ", reduced model ", tags$strong("Tissue + Treatment + Age"), ". Significant genes have different Age trends between Treatment groups."),
                                      tags$li("Time only: full model ", tags$strong("Tissue + Treatment + Age"), ", reduced model ", tags$strong("Tissue + Treatment"), ". Significant genes change with Age. With no Condition Variable: ", tags$strong("Age"), " vs ", tags$strong("1"), ".")
                                    ),
                                    tags$p("Time is treated as a categorical factor (one level per time point). The DEGs from either method and test can be clustered in DEG Clustering Analysis, and any tested gene can be used in Similar Pattern Search.")
                                  )
                                ),
                                br(),
                                actionButton(ns("save_ts_deg"), "Save to output"),
                                br(),br(),
                                tags$style(HTML(paste0("#", ns("ts_deg"), ".shiny-output-error-validation { white-space: pre-line; }"))),
                                DT::DTOutput(ns("ts_deg"))),
                       tabPanel(title="DEG Clustering Analysis",
                                tabsetPanel(id=ns("ts_cluster_tabset"),
                                            tabPanel(title = "Time Series Cluster Plot",
                                                     br(),
                                                     tags$details(
                                                       # This is the clickable header that is always visible
                                                       tags$summary(
                                                         tags$strong("Click to view Clustering Methodology (degPatterns & DIANA)", 
                                                                     style = "color: #28a745; cursor: pointer;")
                                                       ),
                                                       # Everything inside this div is hidden until expanded
                                                       tags$div(
                                                         style = "padding: 15px; background-color: #f4faf6; border-left: 5px solid #28a745; margin-top: 10px;",
                                                         tags$p(tags$strong("The degPatterns() R function in DEGreport library is used to do the DEG clustering.")),
                                                         tags$p(
                                                           "It can work with one or more groups with 2 or more several time points. Before calculating the genes similarity among samples, all samples inside the same time point (time parameter) and group (col parameter) are collapsed together, and the mean value is the representation of the group for the gene abundance. Then, all pair-wise gene expression is calculated using cor.test R function using kendall as the statistical method. A distance matrix is created from those values. After that, cluster::diana() is used for the clustering of gene-gene distance matrix and cut the tree using the divisive coefficient of the clustering, giving as well by diana."
                                                         ),
                                                         tags$p(tags$strong("How cluster::diana() Calculates Clusters")),
                                                         tags$p(
                                                           'The cluster::diana() function, which stands for DIvisive ANAlysis, is a "top-down" hierarchical clustering algorithm. Unlike common "bottom-up" (agglomerative) methods like hclust, it starts with all your data in one single cluster and systematically splits them. The algorithm follows a divisive process to build a hierarchy:'
                                                         ), 
                                                         tags$p(
                                                           tags$strong("Start at the Top:"), 
                                                           " Initially, all observations (genes, in your case) are in one large cluster."
                                                         ),
                                                         tags$p(
                                                           tags$strong('Find the "Splinter" Group:'), 
                                                           " In each step, the algorithm identifies the cluster with the largest diameter (the maximum distance between any two elements)."
                                                         ),
                                                         tags$p(
                                                           tags$strong("The Seed:"), 
                                                           'Within that cluster, it finds the most "disparate" element—the one with the highest average distance to all other elements in the cluster. This element starts a new "splinter group."'
                                                         ),
                                                         tags$p(
                                                           tags$strong("Reassignment:"), 
                                                           " The algorithm then looks at all other elements in the original cluster. If an element is closer to the splinter group than to the remainder of the original group, it is moved to the splinter group."
                                                         ),
                                                         tags$p(
                                                           tags$strong("Iterate:"), 
                                                           " This process repeats until every observation is its own individual cluster, creating a dendrogram (tree)."
                                                         ),
                                                         tags$p(
                                                           tags$strong("Distance Metric:"), 
                                                           " By default, it uses Euclidean distance, though the degPatterns function specifically calculates a distance matrix using Kendall correlation before passing it to diana()."
                                                         )
                                                       )
                                                     ),
                                                     br(),
                                                     actionButton(ns("ts_cluster_plot"), "Save to output"),
                                                     br(),br(),
                                                     plotOutput(ns("ts_plot"), width = '1200px', height = "1200px")),
                                            tabPanel(title="DEG Cluster Result",
                                                     br(),
                                                     DT::DTOutput(ns("ts_cluster_table")),
                                                     tags$script(HTML(sprintf("
                                                     $(document).on('click', '[id^=plotExp_]', function() {
                                                     var genes = $(this).attr('data-genes');
                                                     Shiny.setInputValue('%s', genes, {priority: 'event'});
                                                     });", ns("plotExp_trigger"))))
                                            ),
                                            tabPanel(title="Time‑series Sample Heatmap",
                                                     br(),
                                                     actionButton(ns("ts_tpm_gct"), "Save sample normalized log2TPM GCT data file to output"),
                                                     br(),
                                                     actionButton(ns("ts_plot_sample_heatmap"), "Plot/Refresh", style="color: #0961E3; background-color: #F6E98C ; border-color: #2e6da4"),
                                                     actionButton(ns("save_ts_sample_heatmap"), "Save to output"),
                                                     plotOutput(ns("ts_sample_heatmap"))
                                                     ),
                                            tabPanel(title="Time‑series Group Heatmap",
                                                     br(),
                                                     actionButton(ns("ts_zscore_gct"), "Save group z-score GCT data file to output"),
                                                     br(),
                                                     actionButton(ns("ts_plot_group_heatmap"), "Plot/Refresh", style="color: #0961E3; background-color: #F6E98C ; border-color: #2e6da4"),
                                                     actionButton(ns("save_ts_group_heatmap"), "Save to output"),
                                                     plotOutput(ns("ts_group_heatmap"))
                                                     )
                                )
                       ),
                       tabPanel(title="Similar Pattern Search",
                                br(),
                                tags$details(
                                  tags$summary(
                                    tags$strong("Click to view Similar Pattern Search Methodology",
                                                style = "color: #28a745; cursor: pointer;")
                                  ),
                                  tags$div(
                                    style = "padding: 15px; background-color: #f4faf6; border-left: 5px solid #28a745; margin-top: 10px;",
                                    tags$p("Finds genes whose expression changes over time the same way (or the opposite way) as the seed gene (the gene chosen in Gene to Match). Run the Time-Series DE Analysis first; this uses its samples, time variable and last Condition Variable."),
                                    tags$p("Each gene's expression (VST or log values) is averaged per time point, or per condition and time point when a Condition Variable is selected (the same profile DEG Clustering uses), then every gene's profile is correlated with the seed gene's profile. ",
                                           tags$strong("Kendall"), " (default) only looks at the up/down ordering, the same similarity DEG Clustering uses, so expression level and fold-change size do not matter. ",
                                           tags$strong("Pearson"), " also rewards similar relative step sizes. With few time points Kendall has many ties; ties are broken by Pearson r."),
                                    tags$p("Without a Condition Variable (Time only test), the profile is the time-point means T1..Tn. With a condition selected, the profile covers every condition (e.g. Ctrl T1..Tn then Treated T1..Tn), so a match must follow the seed gene in each condition; this also applies to a Time only test run with a condition selected."),
                                    tags$p(tags$strong("change_distance"), " measures how close the size of the gene's changes is to the seed gene's (log2 scale, lower = closer). Both profiles are centred on their own mean, which removes expression level; change_distance is the root mean square of the differences between the two centred profiles:"),
                                    tags$ul(
                                      tags$li("Example: seed gene means 4.22, 3.62, 3.80 (centred +0.34, -0.26, -0.08); matched gene means 3.75, 3.55, 3.61 (centred +0.12, -0.09, -0.03)."),
                                      tags$li("Differences: -0.23, +0.17, +0.06; change_distance = sqrt((0.051 + 0.029 + 0.003) / 3) = 0.17."),
                                      tags$li("0 = identical changes (only the expression level differs). A gene that does not change at all scores the seed gene's own RMS change, so values close to the largest in the table usually mean the gene barely changes.")
                                    ),
                                    tags$p(tags$strong("Ranking:"), " genes are sorted by the chosen correlation, then Pearson r, both compared at 2 decimals (0.996 and 0.999 count as tied), and then by change_distance. With few time points many genes reach r = 1.00, so change_distance usually decides among them."),
                                    tags$p(tags$strong("How to pick the top genes from the table:")),
                                    tags$ol(
                                      tags$li("Shape: keep genes with r close to 1 (e.g. 0.9 or higher; close to -1 for an Opposite search)."),
                                      tags$li("Size of change: among those, prefer a small change_distance. change_distance is how far the gene's centred profile is from the seed gene's, so a larger value means it is further from the seed gene."),
                                      tags$li("Evidence: check DE_padj. A gene with a small DE_padj really changes over time; a large DE_padj means its profile may be noise even if it matches well."),
                                      tags$li("Abundance: check baseMean (DESeq2) or AveExpr (limma). Very low values are noisier."),
                                      tags$li("Search Among Time-series DEGs only to apply the DE filter automatically, then use the table's column sorting and search box, Copy seed + selected genes, or Create Expression Plot.")
                                    ),
                                    tags$p(tags$strong("Search Among:"), " all genes tested in the DE run, or only the Time-series DEGs passing the cutoffs above. Each gene's DE padj is shown for reference."),
                                    tags$p(tags$strong("Expression used:"), " VST for the DESeq2 method, the project's own log values for the limma method."),
                                    tags$p(tags$strong("VST"), " (variance stabilizing transformation, DESeq2::vst) puts count-like values on a log2-like scale: the app's log values (e.g. log2 TPM) are converted back to the linear scale (e.g. TPM), rounded and treated as counts (raw counts are not available), normalized with DESeq2 size factors, and transformed with vst(blind = TRUE). For moderate and high values VST is close to log2 of the normalized values; for low values it shrinks values together so noisy low-expression genes do not swing widely. Differences of 1 are therefore about a 2-fold change, and changes of low-expression genes are slightly compressed. The group means in the plot and in change_distance are VST means.")
                                  )
                                ),
                                br(),
                                actionButton(ns("save_sim_plot"), "Save plot to output"),
                                actionButton(ns("save_sim_table"), "Save table to output"),
                                br(), br(),
                                span(textOutput(ns("sim_msg")), style = "color:blue; font-size:15px; font-family:arial; font-style:italic"),
                                # The plot area appears only after the first search, so the
                                # tab has no large empty box before that.
                                conditionalPanel(ns = ns, "input.compute_similar > 0",
                                  plotly::plotlyOutput(ns("sim_plot"), width = "1000px", height = "500px")
                                ),
                                br(), br(),
                                uiOutput(ns("sim_actions")),
                                br(),
                                DT::DTOutput(ns("sim_table"))
                       )
           )         
    )
)}

# ---- Server Module ----
TimeSeries_server <- function(id, parent_session) {
  moduleServer(id,
               function(input, output, session) {
                 # See the matching comment in correlation.R -- actionButton
                 # click counts get bookmarked and restoring a non-zero count
                 # silently re-triggers that button's own observeEvent.
                 session$setBookmarkExclude(c("compute_DE", "compute_cluster", "save_ts_deg",
                                               "ts_cluster_plot", "ts_tpm_gct", "ts_plot_sample_heatmap",
                                               "ts_zscore_gct", "ts_plot_group_heatmap",
                                               "compute_similar", "save_sim_plot", "save_sim_table",
                                               "save_ts_sample_heatmap", "save_ts_group_heatmap"))

                 # ts_hm_annot_group/ts_hm_annot_sample only get populated once
                 # the user re-runs clustering (compute_cluster, deliberately
                 # excluded from bookmark restore above), which can happen well
                 # after any fixed restore window -- captured here and consumed
                 # by the observeEvent(DataClusterReactive(), ...) below whenever
                 # that eventually fires, same pattern as WGCNA's base_<attr>.
                 restored_ts_hm_annot <- reactiveVal(NULL)

                 # sel_time_var/sel_condition_var get their choices populated
                 # (and selected forced to NULL) by observeEvent(DataQCReactive(),
                 # ...) below, which also fires again at restore time. Reapply
                 # unconditionally for a few seconds after restore, same pattern
                 # as app.R.
                 session$onRestored(function(state) {
                   restored <- state$input
                   if (!is.null(restored$ts_hm_annot_group) || !is.null(restored$ts_hm_annot_sample)) {
                     restored_ts_hm_annot(list(group = restored$ts_hm_annot_group, sample = restored$ts_hm_annot_sample))
                   }

                   pending <- list(
                     sel_time_var = list(value = restored$sel_time_var, apply = function(v) {
                       req(DataQCReactive())
                       MetaData <- DataQCReactive()$MetaData
                       MetaData_clean <- MetaData[, sapply(MetaData, function(x) length(unique(x)) > 1)]
                       attrs <- sort(setdiff(colnames(MetaData_clean), c("sampleid", "Order", "ComparePairs")))
                       updateSelectInput(session, "sel_time_var", choices = attrs, selected = v)
                     }),
                     ts_method = list(value = restored$ts_method, apply = function(v) {
                       updateRadioButtons(session, "ts_method", selected = v)
                     }),
                     sel_condition_var = list(value = restored$sel_condition_var, apply = function(v) {
                       req(DataQCReactive())
                       MetaData <- DataQCReactive()$MetaData
                       MetaData_clean <- MetaData[, sapply(MetaData, function(x) length(unique(x)) > 1)]
                       attrs <- sort(setdiff(colnames(MetaData_clean), c("sampleid", "Order", "ComparePairs")))
                       updateSelectizeInput(session, "sel_condition_var", choices = attrs, selected = v)
                     })
                   )
                   pending <- Filter(function(p) !is.null(p$value), pending)
                   if (length(pending) > 0) {
                     attempts_left <- 10  # ~3s at 300ms
                     restore_observer <- NULL
                     restore_observer <- observe({
                       invalidateLater(300, session)
                       isolate({
                         attempts_left <<- attempts_left - 1
                         for (nm in names(pending)) {
                           p <- pending[[nm]]
                           tryCatch(p$apply(p$value), error = function(e) NULL)
                         }
                         if (attempts_left <= 0) restore_observer$destroy()
                       })
                     })
                   }
                 })

                 linear_small_value <- reactiveVal()
                 linear_base <- reactiveVal()
                 time_col <- reactiveVal()
                 cond_col <- reactiveVal()
                 df_normalized_log2TPM <- reactiveVal()
                 df_cluster_normalized_wide <- reactiveVal()
                 cluster_genes_ordered <- reactiveVal()
                 df_cluster_sample_meta <- reactiveVal()
                 df_cluster_group_meta <- reactiveVal()
                 
                 output$selectGroupSampleTimeSeries <- renderUI(shared_header_content())

                 observe( {
                   expU = exp_unit()
                   linear_base(as.numeric(str_replace(str_extract(expU,"log\\d+"),"log","")))
                   linear_small_value(as.numeric(str_replace(str_split_fixed(expU, "\\+", 2)[2], "\\)", "")))
                   unit=str_replace_all(str_extract(expU, "\\(.+\\+"), "(\\(|\\+)", "")
                 })
                 
                 # Default the method from the expression unit: count-like
                 # units (TPM/FPKM/CPM/counts) -> DESeq2, anything else
                 # (e.g. log2 intensity) -> limma.
                 observeEvent(exp_unit(), {
                   expU <- exp_unit()
                   req(expU)
                   method <- if (grepl("TPM|FPKM|RPKM|CPM|count", expU, ignore.case = TRUE)) "deseq2" else "limma"
                   updateRadioButtons(session, "ts_method", selected = method)
                 })

                 observeEvent(DataQCReactive(), {
                   req(DataQCReactive())
                   DataIn = DataQCReactive()
                   MetaData=DataIn$MetaData
                   MetaData_clean <- MetaData[, sapply(MetaData, function(x) length(unique(x)) > 1)]
                   attributes=sort(setdiff(colnames(MetaData_clean), c("sampleid", "Order", "ComparePairs")))
                   updateSelectInput(session, "sel_time_var", choices=attributes, selected=NULL)  
                   updateSelectizeInput(session,'sel_condition_var',choices=attributes, selected=NULL)
                 })
                 
                 DEGReactive <- eventReactive(input$compute_DE, {
                   withProgress(message = 'Processing...', value = 0, {
                     log_info("Time Course Analysis: computing DEGs (test=", input$ts_test, ", time=", input$sel_time_var,
                              ", condition=", paste(input$sel_condition_var, collapse=","), ").")
                     req(DataQCReactive())
                     req(input$sel_time_var)
                     DataIn = DataQCReactive()
                     data_long <- DataIn$tmp_data_long
                     MetaData <- DataIn$MetaData
                     ProteinGeneName <- DataIn$ProteinGeneName

                     if (identical(input$ts_method, "limma")) {
                       # limma works on the app's log values directly; missing
                       # values stay NA (not 0 as in the DESeq2 count matrix).
                       expr_mtx <- data_long %>%
                         dplyr::select(UniqueID, sampleid, expr) %>%
                         tidyr::pivot_wider(names_from = sampleid, values_from = expr) %>%
                         tibble::column_to_rownames("UniqueID") %>%
                         as.matrix()
                       min_n <- suppressWarnings(as.integer(input$ts_min_per_group))
                       if (is.na(min_n) || min_n < 1) min_n <- 1
                       res <- run_timecourse_limma(MetaData,
                                                   expr_mtx,
                                                   time_variable = input$sel_time_var,
                                                   Condition = input$sel_condition_var,
                                                   test = input$ts_test,
                                                   min_per_group = min_n)
                     } else {
                       count_mtx <- get_count_mtx(data_long, linear_base(), linear_small_value())

                       res <- run_timecourse_DEG(MetaData,
                                          count_mtx,
                                          time_variable = input$sel_time_var,
                                          Condition = input$sel_condition_var,
                                          test = input$ts_test
                                          )
                     }
                     # Remember the variables this result was built with, so
                     # clustering/similar-gene search stay consistent with it
                     # even if the inputs change afterwards.
                     conds <- setdiff(input$sel_condition_var, input$sel_time_var)
                     res$time_var <- input$sel_time_var
                     res$cond_var <- if (length(conds) > 0) tail(conds, 1) else ts_no_condition_col
                     if (identical(res$status, "success")) {
                       log_info("Time Course Analysis: DEG computation succeeded.")
                     } else {
                       log_warn("Time Course Analysis: DEG computation failed: ", res$msg)
                     }
                     res
                   })
                 })
                 
                 filtered_DE <- reactive({
                   req(DEGReactive())
                   res <- DEGReactive()
                   
                   if (res$status == "success") {
                     df<- res$results %>%
                       data.frame() %>%
                       rownames_to_column("gene") %>%
                       as_tibble() %>%
                       filter(padj < input$padjCutoff) %>%
                       arrange(padj) %>%
                       head(n = input$top_n)
                     
                     if (input$ts_genelable == 'Gene.Name') {
                       req(DataQCReactive())
                       ProteinGeneName <- DataQCReactive()$ProteinGeneName
                       df <- df %>% 
                         dplyr::left_join(ProteinGeneName %>% dplyr::select(UniqueID, Gene.Name), by = c("gene" = "UniqueID")) %>%
                         dplyr:: rename(UniqueID = gene) %>%
                         dplyr:: filter(!is.na(Gene.Name), Gene.Name != "") %>%
                         dplyr:: group_by(Gene.Name) %>%
                         dplyr:: slice_max(.data[[res$mean_col]], n = 1, with_ties = FALSE) %>%
                         dplyr:: ungroup() %>%
                         dplyr:: rename(gene = Gene.Name) %>%
                         dplyr:: relocate(gene, .before = 1)
                     }
                     
                     nrow = nrow(df)
                     test_name <- paste0(ifelse(identical(res$test, "time"), "Time only", "Time x Condition"),
                                         ", ", ifelse(identical(res$method, "limma"), "limma", "DESeq2"))
                     if (identical(res$method, "limma") && res$n_filtered > 0)
                       test_name <- sprintf("%s; %d genes/proteins with < %d values in a group left out", test_name, res$n_filtered, res$min_per_group)
                     max_padj = if (nrow > 0) max(df$padj) else NA
                     list(
                       df = df,
                       msg = if (nrow == 0) sprintf("No genes pass padj < %s (%s).", input$padjCutoff, test_name) else
                         sprintf("Top %s DEGs (%s) are selected with maximum padj %.3f.", nrow, test_name, max_padj)
                     )
                   } else {
                     list(
                       df = NULL,
                       msg = res$msg
                     )
                   }
                 })
                     
                 # Errors: first line here, the full multi-line message in the DE table area.
                 # Read the message first: inside sub() (an S4 generic here) a pending req()
                 # would surface as an error instead of a blank box.
                 output$filtered_DEG <- renderText({
                   msg <- filtered_DE()$msg
                   sub("\n.*", "", msg)
                 })
                     
                 observeEvent(input$save_ts_deg, {
                   log_info("Time Course Analysis: DEG data table saved to output.")
                   df <- filtered_DE()$df
                   if (input$ts_genelable == 'Gene.Name') {
                     df <- df %>%
                       dplyr::select(-UniqueID)
                   }
                   # Own sheet name: "DEG_data" is used by the DEGs tab and
                   # would overwrite (or be overwritten by) that table.
                   saved_table$TimeSeries_DE <- as.data.frame(df)
                   showNotification("Time-Series DE table saved to output.", type = "message")
                 })
                 
                 output$ts_deg<- DT::renderDT({
                   req(filtered_DE())
                   df <- filtered_DE()$df
                   
                   # Use validate to show the (multi-line) error message in the UI
                   validate(need(!is.null(df), filtered_DE()$msg))
                   
                   if (input$ts_genelable == 'Gene.Name') {
                     df <- df %>%
                       dplyr::select(-UniqueID)
                   }
                   num_cols <- names(df)[sapply(df, is.numeric)]
                   DT::datatable(df, extensions = 'Buttons', options = list(
                     dom = "Blfrtip", buttons = c("csv", "excel", "print"), pageLength = 20), rownames= FALSE) %>% 
                     formatSignif(columns= num_cols, digits=3)
                 })
                 
                 DataClusterReactive <- eventReactive(input$compute_cluster, {
                   withProgress(message = 'Processing...', value = 0, {
                     log_info("Time Course Analysis: computing DEG clusters.")
                     req(DEGReactive())
                     req(filtered_DE())
                     
                     res <- DEGReactive()
                     filtered_dds <- filtered_DE()$df
                     sig_gene <- filtered_dds$gene
                     if (input$ts_genelable == 'Gene.Name') {
                       sig_gene <- filtered_dds$UniqueID
                     }
                     time_col(res$time_var)
                     cond_col(res$cond_var)
                     res_cluster <- run_DEG_cluster(res$expr, res$sample_meta, sig_gene, cond_col())
                     
                     df_cluster_normalized <- res_cluster$clusters$normalized
                     df_counts <- res_cluster$clusters$counts
                     df_counts_unique <- df_counts[unique(df_cluster_normalized$genes), ]

                     df_wide <- df_cluster_normalized %>%
                       dplyr::select(genes, merge, value, cluster) %>%
                       pivot_wider(names_from = merge, values_from = value) %>%
                       tibble::column_to_rownames('genes')
                     
                     df_cluster <- res_cluster$clusters$df
                     genes_ordered <- df_cluster[order(df_cluster$cluster), "genes"]
                     
                     # No condition selected -> constant AllSamples column, which has no tmp_group entry.
                     cond_levels <- DataQCReactive()$tmp_group[[cond_col()]]
                     if (is.null(cond_levels)) cond_levels <- unique(as.character(df_cluster_normalized[[cond_col()]]))
                     group_meta <- df_cluster_normalized %>% 
                       dplyr::select(all_of(c('merge', time_col(), "Unified_Time", cond_col()))) %>% 
                       distinct() %>% 
                       dplyr::mutate(!!cond_col() := factor(!!sym(cond_col()), levels = cond_levels)) %>%
                       dplyr::arrange(!!sym(cond_col()), Unified_Time) %>%
                       dplyr::mutate(!!cond_col() := as.character(!!sym(cond_col())))
                     
                     df_wide <- df_wide[genes_ordered, c('cluster', group_meta$merge)]
                     df_cluster_normalized_wide(df_wide)
                     
                     df_counts_unique <- df_counts_unique[genes_ordered, sample_order()]
                     df_normalized_log2TPM(df_counts_unique)
                     
                     cluster_genes_ordered(genes_ordered)
                     df_cluster_sample_meta(res_cluster$sample_meta)
                     df_cluster_group_meta(group_meta)
                     res_cluster
                   })
                 })
                 
                 cluster_table <- eventReactive(input$compute_cluster, {
                   withProgress(message = "Processing...", value = 0, {
                     req(DataClusterReactive())
                     req(DataQCReactive())
                     ProteinGeneName <- DataQCReactive()$ProteinGeneName
                     DataCluster <-DataClusterReactive()
                     clusters <- DataCluster$clusters
                     df = clusters$df

                     if (input$ts_genelable == 'Gene.Name') {
                       df = df %>% 
                         dplyr::left_join(ProteinGeneName %>% dplyr::select(UniqueID, Gene.Name), by = c("genes" = "UniqueID")) %>%
                         dplyr:: select(-genes) %>%
                         dplyr:: filter(!is.na(Gene.Name), Gene.Name != "") %>%
                         dplyr:: rename(genes = Gene.Name) %>%
                         dplyr:: relocate(genes, .before = 1)
                     }
                     df = df %>%
                       dplyr:: group_by(cluster) %>%
                       dplyr:: summarise(
                         n_gene = n(),
                         genes = paste(genes, collapse = ","), .groups = "drop"
                       ) %>%   
                       dplyr::mutate(
                         copy = purrr::map_chr(seq_len(n()), ~ paste0(
                           as.character(
                             rclipButton(
                               paste0("clipbtn_", .x),
                               label = "Copy all genes",
                               clipText = genes[.x],
                               icon = icon("copy", lib = "glyphicon"),
                               class = "btn-primary btn-sm"
                             )
                           ),
                           " ",
                           sprintf(
                             '<button id="plotExp_%s" class="btn btn-info btn-sm" data-genes="%s">Create Expression Plot</button>',
                             .x,
                             genes[.x]
                           )
                         ))
                       )  %>%
                       dplyr::select(cluster, n_gene, copy, genes)
                   })
                 })
                 
                 output$ts_cluster_table<- DT::renderDT({
                   req(cluster_table())
                   df <- cluster_table()
                   DT::datatable(df, 
                                 rownames = FALSE,
                                 escape = FALSE,
                                 selection = "none",
                                 colnames=c("Cluster", "Number of genes", "Action","Genes in cluster"),
                                 extensions = 'Buttons', 
                                 options = list(dom = "Blfrtip", 
                                                buttons = c("csv", "excel", "print")
                                                )
                                 ) 
                 })
                 
                 observeEvent(input$plotExp_trigger, {
                   gene_list <- input$plotExp_trigger
                   updateNavbarPage(parent_session, inputId = "menu", selected = "Exp_Plot")
                   updateTabsetPanel(parent_session, inputId = "expression_tabset", selected = "Searched Expression Data")
                   parent_session$onFlushed(function() {
                     updateSelectInput(parent_session, inputId = "exp_subset", selected = "Upload Genes")
                     updateTextAreaInput(parent_session, inputId = "exp_list", value = gene_list)
                   }, once = TRUE)
                 })
                 
                 cluster_plot <- eventReactive(input$compute_cluster, {
                   withProgress(message = "Processing...", value = 0, {
                     req(DataClusterReactive())
                     DataCluster <-DataClusterReactive()
                     clusters <- DataCluster$clusters
                     
                     interaction_cov <- cond_col()
                     
                     clusters$plot$data[[interaction_cov]] <- factor(
                       clusters$plot$data[[interaction_cov]], 
                       levels = unique(clusters$plot$data[[interaction_cov]])
                     )
                     clusters$plot
                   })
                 })
                 
                 output$ts_plot<- renderPlot({
                   cluster_plot() +
                     theme(
                       text = element_text(size = input$base_text),        # all text
                       axis.text = element_text(size = input$axis_text),   # tick labels
                       axis.title = element_text(size = input$axis_title),  # axis titles
                       strip.text = element_text(size = input$strip_text),  # facet titles
                       legend.text = element_text(size = input$legend_text),
                       legend.title = element_text(size = input$legend_title)
                     )
                 },
                 width  = function() input$plot_width,
                 height = function() input$plot_height
                 )
                 
                 observeEvent(input$ts_cluster_plot, {
                   log_info("Time Series Cluster Plot saved to output.")
                   saved.num <- length(saved_plots$ts_cluster_plot) + 1
                   p <- cluster_plot() +
                     theme(
                       text = element_text(size = input$base_text),        # all text
                       axis.text = element_text(size = input$axis_text),   # tick labels
                       axis.title = element_text(size = input$axis_title),  # axis titles
                       strip.text = element_text(size = input$strip_text),  # facet titles
                       legend.text = element_text(size = input$legend_text),
                       legend.title = element_text(size = input$legend_title)
                     )
                   
                   saved_plots$ts_cluster_plot[[saved.num]] <- p
                 })
                 
                 ############## Similar Pattern Search
                 # Expression of every gene tested in the current DE run (VST
                 # for DESeq2, the app's log values for limma).
                 ts_vst <- reactive({
                   req(DEGReactive())
                   res <- DEGReactive()
                   req(identical(res$status, "success"))
                   res$expr
                 })

                 ts_gene_choices <- reactive({
                   req(DEGReactive())
                   res <- DEGReactive()
                   req(identical(res$status, "success"))
                   ids <- rownames(res$expr)
                   labels <- ids
                   if (input$ts_genelable == 'Gene.Name') {
                     ProteinGeneName <- DataQCReactive()$ProteinGeneName
                     gn <- ProteinGeneName$Gene.Name[match(ids, ProteinGeneName$UniqueID)]
                     labels <- ifelse(is.na(gn) | gn == "", ids, gn)
                     # A gene name shared by several IDs needs the ID to tell them apart.
                     dup <- labels %in% labels[duplicated(labels)]
                     labels[dup] <- paste0(labels[dup], " (", ids[dup], ")")
                   }
                   o <- order(labels)
                   setNames(ids[o], labels[o])
                 })

                 observeEvent(ts_gene_choices(), {
                   choices <- ts_gene_choices()
                   # Keep the current gene if it is still in this run.
                   sel <- isolate(input$sim_gene)
                   if (is.null(sel) || !(sel %in% choices)) sel <- character(0)
                   updateSelectizeInput(session, "sim_gene", choices = choices, selected = sel, server = TRUE)
                 })

                 SimilarReactive <- eventReactive(input$compute_similar, {
                   validate(need(!is.null(input$compute_DE) && input$compute_DE > 0, "Please run the Time-Series DE Analysis first."))
                   req(DEGReactive())
                   res <- DEGReactive()
                   validate(need(identical(res$status, "success"), "The Time-Series DE Analysis did not finish; please fix it first."))
                   validate(need(!is.null(input$sim_gene) && input$sim_gene != "", "Please select a gene to match."))
                   withProgress(message = "Finding similar genes...", value = 0, {
                     anchor <- input$sim_gene
                     vst_mat <- ts_vst()
                     validate(need(anchor %in% rownames(vst_mat), "The seed gene was not tested in the current DE run."))
                     sample_meta <- res$sample_meta
                     prof <- ts_group_profiles(vst_mat, sample_meta, res$cond_var)
                     validate(need(sd(prof$profiles[anchor, ]) > 0, "The seed gene has the same mean expression at every time point."))
                     candidates <- NULL
                     if (input$sim_pool == "deg") {
                       df_deg <- filtered_DE()$df
                       candidates <- if ("UniqueID" %in% names(df_deg)) df_deg$UniqueID else df_deg$gene
                       validate(need(length(candidates) > 0, "No Time-series DEGs pass the current cutoffs."))
                     }
                     sim <- find_similar_genes(prof$profiles, anchor, method = input$sim_method,
                                               direction = input$sim_direction, candidates = candidates)
                     lrt <- res$results
                     ProteinGeneName <- DataQCReactive()$ProteinGeneName
                     sim <- sim %>%
                       dplyr::mutate(Gene.Name = ProteinGeneName$Gene.Name[match(UniqueID, ProteinGeneName$UniqueID)],
                                     DE_padj = lrt[UniqueID, "padj"],
                                     mean_expr = lrt[UniqueID, res$mean_col]) %>%
                       dplyr::select(rank, UniqueID, Gene.Name, r, pearson_r, change_distance, DE_padj, mean_expr) %>%
                       head(ts_pos_int(input$sim_top_n, 50))
                     # With Pearson as the method, r already is the Pearson r;
                     # drop the tiebreak copy so the table has one pearson_r.
                     if (identical(input$sim_method, "pearson")) sim$pearson_r <- NULL
                     names(sim)[names(sim) == "r"] <- paste0(input$sim_method, "_r")
                     names(sim)[names(sim) == "mean_expr"] <- res$mean_col
                     log_info("Time Course Analysis: similar pattern search for ", anchor, " (", input$sim_method, ", ",
                              input$sim_direction, ", pool=", input$sim_pool, ") returned ", nrow(sim), " genes.")
                     list(table = sim, profiles = prof$profiles, groups = prof$groups, anchor = anchor,
                          cond_var = res$cond_var, time_var = res$time_var, method = input$sim_method,
                          direction = input$sim_direction, expr_label = res$expr_label)
                   })
                 })

                 sim_label <- function(ids) {
                   if (input$ts_genelable == 'Gene.Name') {
                     ProteinGeneName <- DataQCReactive()$ProteinGeneName
                     gn <- ProteinGeneName$Gene.Name[match(ids, ProteinGeneName$UniqueID)]
                     ifelse(is.na(gn) | gn == "", ids, gn)
                   } else ids
                 }

                 output$sim_msg <- renderText({
                   sim <- SimilarReactive()
                   grp_name <- if (identical(sim$cond_var, ts_no_condition_col)) "time-point" else "condition x time"
                   msg <- sprintf("%s genes with the most %s pattern to %s (%s on %s means of %s).",
                                  nrow(sim$table), ifelse(sim$direction == "opposite", "opposite", "similar"),
                                  sim_label(sim$anchor), tools::toTitleCase(sim$method), grp_name, sim$expr_label)
                   if (ncol(sim$profiles) < 3) msg <- paste0(msg, " Only 2 ", grp_name, " groups: every correlation is +1 or -1, so genes are ranked by direction and then by how close their change size is (change_distance).")
                   msg
                 })

                 # Long table of the plotted profiles, shared by the static
                 # ggplot (saved to output) and the interactive plotly view.
                 sim_plot_data <- reactive({
                   sim <- SimilarReactive()
                   n_plot <- min(ts_pos_int(input$sim_plot_n, 10), nrow(sim$table))
                   ids <- c(sim$anchor, head(sim$table$UniqueID, n_plot))
                   pm <- sim$profiles[ids, , drop = FALSE]
                   z <- if (input$sim_plot_scale == "zscore") t(scale(t(pm))) else pm - rowMeans(pm)
                   groups <- sim$groups
                   df <- data.frame(UniqueID = rep(ids, times = ncol(z)),
                                    label = groups$label[rep(seq_len(ncol(z)), each = length(ids))],
                                    z = as.vector(z), stringsAsFactors = FALSE)
                   df <- dplyr::left_join(df, groups, by = "label")
                   df$condition <- factor(df$condition, levels = unique(groups$condition))
                   lab <- sim_label(ids)
                   # Two IDs with the same gene name get the ID appended so
                   # legend entries and hover labels stay distinguishable.
                   dup <- lab %in% lab[duplicated(lab)] & lab != ids
                   lab[dup] <- paste0(lab[dup], " (", ids[dup], ")")
                   df$gene_label <- lab[match(df$UniqueID, ids)]
                   val_name <- if (input$sim_plot_scale == "zscore") "z-score" else "centred value"
                   tab <- sim$table
                   r_col <- paste0(sim$method, "_r")
                   i <- match(df$UniqueID, tab$UniqueID)
                   df$hover <- ifelse(df$UniqueID == sim$anchor,
                                      sprintf("<b>%s</b> (seed gene)<br>%s<br>%s: %.3g", df$gene_label, df$label, val_name, df$z),
                                      sprintf("<b>%s</b><br>rank %s, %s = %.3f, change_distance = %.3g<br>%s<br>%s: %.3g",
                                              df$gene_label, tab$rank[i], r_col, tab[[r_col]][i], tab$change_distance[i],
                                              df$label, val_name, df$z))
                   list(df = df, n_plot = n_plot, anchor = sim$anchor, time_var = sim$time_var,
                        ylab = if (input$sim_plot_scale == "zscore") paste("Z-score of mean", sim$expr_label) else paste("Mean", sim$expr_label, "(centred)"),
                        title = paste0(sim_label(sim$anchor), " and top ", n_plot, " ",
                                       ifelse(sim$direction == "opposite", "opposite", "similar"), " genes"))
                 })

                 sim_plot <- reactive({
                   sim <- SimilarReactive()
                   n_plot <- min(ts_pos_int(input$sim_plot_n, 10), nrow(sim$table))
                   ids <- c(sim$anchor, head(sim$table$UniqueID, n_plot))
                   # Z-scores of the best matches are often identical to the
                   # seed gene's (that is what a near-perfect correlation
                   # means), so they would sit exactly under its line. The
                   # centred scale keeps each gene's change size apart.
                   pm <- sim$profiles[ids, , drop = FALSE]
                   z <- if (input$sim_plot_scale == "zscore") t(scale(t(pm))) else pm - rowMeans(pm)
                   expr_label <- sim$expr_label
                   groups <- sim$groups
                   df <- data.frame(UniqueID = rep(ids, times = ncol(z)),
                                    label = groups$label[rep(seq_len(ncol(z)), each = length(ids))],
                                    z = as.vector(z), stringsAsFactors = FALSE)
                   df <- dplyr::left_join(df, groups, by = "label")
                   df$condition <- factor(df$condition, levels = unique(groups$condition))
                   df$Gene <- factor(ifelse(df$UniqueID == sim$anchor, "Seed gene", "Matched genes"),
                                     levels = c("Matched genes", "Seed gene"))
                   p <- ggplot(df, aes(x = time, y = z, group = UniqueID)) +
                     geom_line(data = df[df$Gene == "Matched genes", ], aes(color = Gene), alpha = 0.7, linewidth = 0.8) +
                     geom_point(data = df[df$Gene == "Matched genes", ], aes(color = Gene), alpha = 0.7, size = 1.5) +
                     # Dashed so matched lines lying right underneath stay visible.
                     geom_line(data = df[df$Gene == "Seed gene", ], aes(color = Gene), linewidth = 1.4, linetype = "dashed") +
                     scale_color_manual(values = c("Matched genes" = "grey40", "Seed gene" = "red")) +
                     labs(x = sim$time_var, color = NULL,
                          y = if (input$sim_plot_scale == "zscore") paste("Z-score of mean", expr_label) else paste("Mean", expr_label, "(centred)"),
                          title = paste0(sim_label(sim$anchor), " and top ", n_plot, " ",
                                         ifelse(sim$direction == "opposite", "opposite", "similar"), " genes")) +
                     theme_bw(base_size = 16) + theme(legend.position = "bottom")
                   if (length(unique(groups$condition)) > 1) p <- p + facet_wrap(~ condition, nrow = 1)
                   p
                 })

                 # Interactive view: one trace per gene so hovering names the
                 # gene; hovering a line highlights it and dims the others.
                 output$sim_plot <- plotly::renderPlotly({
                   pd <- sim_plot_data()
                   df <- pd$df
                   conds <- levels(df$condition)
                   genes <- c(setdiff(unique(df$UniqueID), pd$anchor), pd$anchor)
                   panels <- lapply(conds, function(cd) {
                     d <- df[df$condition == cd, ]
                     p <- plotly::plot_ly()
                     for (g in genes) {
                       dg <- d[d$UniqueID == g, ]
                       dg <- dg[order(dg$time), ]
                       if (nrow(dg) == 0) next
                       is_anchor <- g == pd$anchor
                       p <- plotly::add_trace(p, x = as.character(dg$time), y = dg$z, type = "scatter",
                                              mode = if (is_anchor) "lines" else "lines+markers",
                                              name = if (is_anchor) paste(dg$gene_label[1], "(seed)") else dg$gene_label[1],
                                              legendgroup = g, showlegend = cd == conds[1],
                                              text = dg$hover, hoverinfo = "text",
                                              line = list(color = if (is_anchor) "red" else "rgba(102,102,102,0.8)",
                                                          width = if (is_anchor) 3 else 1.5,
                                                          dash = if (is_anchor) "dash" else "solid"),
                                              marker = list(color = "rgba(102,102,102,0.8)", size = 5))
                     }
                     plotly::layout(p, xaxis = list(title = pd$time_var, type = "category",
                                                    categoryorder = "array", categoryarray = levels(df$time)),
                                    annotations = if (length(conds) > 1) list(list(text = cd, x = 0.5, y = 1.06, xref = "paper",
                                                                                   yref = "paper", showarrow = FALSE, xanchor = "center")))
                   })
                   fig <- if (length(panels) > 1) plotly::subplot(panels, nrows = 1, shareY = TRUE, titleX = TRUE) else panels[[1]]
                   fig <- plotly::layout(fig, title = list(text = pd$title, x = 0.02),
                                         yaxis = list(title = pd$ylab),
                                         legend = list(title = list(text = "Gene")),
                                         hovermode = "closest")
                   htmlwidgets::onRender(fig, "
                     function(el) {
                       var base = el.data.map(function(t) { return t.line.width; });
                       el.on('plotly_hover', function(e) {
                         var grp = el.data[e.points[0].curveNumber].legendgroup;
                         var op = el.data.map(function(t) { return t.legendgroup === grp ? 1 : 0.15; });
                         var w = el.data.map(function(t, k) { return t.legendgroup === grp ? base[k] + 2 : base[k]; });
                         Plotly.restyle(el, {'opacity': op, 'line.width': w});
                       });
                       el.on('plotly_unhover', function() {
                         Plotly.restyle(el, {'opacity': 1, 'line.width': base});
                       });
                     }")
                 })

                 observeEvent(input$save_sim_table, {
                   if (is.null(input$compute_similar) || input$compute_similar == 0) {
                     showNotification("Click Find Similar Genes first.", type = "warning")
                     return()
                   }
                   sim <- SimilarReactive()
                   log_info("Time Course Analysis: similar pattern table saved to output.")
                   saved_table$TimeSeries_SimilarPattern <- data.frame(seed_gene = sim_label(sim$anchor), sim$table,
                                                                       check.names = FALSE)
                   showNotification("Similar Pattern table saved to output.", type = "message")
                 })

                 observeEvent(input$save_sim_plot, {
                   log_info("Time Course Analysis: similar pattern plot saved to output.")
                   saved.num <- length(saved_plots$ts_similar_plot) + 1
                   saved_plots$ts_similar_plot[[saved.num]] <- sim_plot()
                 })

                 output$sim_actions <- renderUI({
                   sim <- SimilarReactive()
                   genes <- paste(unique(sim_label(c(sim$anchor, sim$table$UniqueID))), collapse = ",")
                   tagList(
                     rclipButton("clipbtn_sim", label = "Copy seed + selected genes", clipText = genes,
                                 icon = icon("copy", lib = "glyphicon"), class = "btn-primary btn-sm"),
                     # Uses the same plotExp_ click handler as the DEG Cluster Result table.
                     HTML(sprintf('<button id="plotExp_sim" class="btn btn-info btn-sm" data-genes="%s">Create Expression Plot</button>',
                                  paste(unique(sim_label(c(sim$anchor, head(sim$table$UniqueID, ts_pos_int(input$sim_plot_n, 10))))), collapse = ",")))
                   )
                 })

                 output$sim_table <- DT::renderDT({
                   sim <- SimilarReactive()
                   df <- sim$table
                   num_cols <- names(df)[sapply(df, is.numeric) & names(df) != "rank"]
                   DT::datatable(df, rownames = FALSE, extensions = 'Buttons',
                                 options = list(dom = "Blfrtip", buttons = c("csv", "excel", "print"), pageLength = 20)) %>%
                     formatSignif(columns = num_cols, digits = 3)
                 })

                 ##############                            
                 ts_TPM_gct <- reactive({
                   DataIn = DataQCReactive()
                   ProteinGeneName <- DataIn$ProteinGeneName

                   DataCluster <-DataClusterReactive()
                   clusters <- DataCluster$clusters
                   sample_meta <- df_cluster_sample_meta() # DataCluster$sample_meta
                   df_counts <- df_normalized_log2TPM()
                   
                   df_cluster <- clusters$df
                   genes_ordered <- df_cluster[order(df_cluster$cluster), "genes"]
                   
                   row_meta <- as.data.frame(df_counts) %>%
                     tibble::rownames_to_column('UniqueID') %>%
                     dplyr::select(UniqueID) %>%
                     left_join(df_cluster, , by = c("UniqueID" = "genes")) %>%
                     left_join(ProteinGeneName, by = 'UniqueID')
                   
                   gct_data <- create_gct_object(df_counts, row_meta, sample_meta)
                   gct_data
                 })
                 
                 observeEvent(input$ts_tpm_gct, {
                   log_info("Time Series cluster sample normalized log2TPM GCT saved to output.")
                   saved_gcts$ts_tpm_gct <- ts_TPM_gct()
                 })
                 
                 ts_zscore_gct <- reactive({
                   DataIn = DataQCReactive()
                   ProteinGeneName <- DataIn$ProteinGeneName
                   
                   DataCluster <-DataClusterReactive()
                   clusters <- DataCluster$clusters

                   df_cluster_wide <- df_cluster_normalized_wide()
                   
                   cluster_genes <- cluster_genes_ordered()
                   
                   row_meta <- df_cluster_wide %>%
                     tibble::rownames_to_column('UniqueID') %>%
                     dplyr::select(UniqueID, cluster) %>%
                     left_join(ProteinGeneName, by = 'UniqueID')
                   
                   df_cluster_wide$cluster <-NULL
                   
                   group_meta <- df_cluster_group_meta()
                   col_meta = group_meta[, c(time_col(), "Unified_Time", cond_col(), 'merge')] 
                   
                   gct_data <- create_gct_object(df_cluster_wide, row_meta, col_meta)
                   gct_data
                 })
                 
                 observeEvent(input$ts_zscore_gct, {
                   log_info("Time Series cluster group z-score GCT saved to output.")
                   saved_gcts$ts_zscore_gct <- ts_zscore_gct()
                 })
                 ############## 
                 observeEvent(DataClusterReactive(), {
                   req(DataClusterReactive())
                   time_col <- time_col()
                   cond_col <- cond_col()

                   # Prefer a bookmarked selection over the plain defaults
                   # below, if one exists -- see restored_ts_hm_annot's
                   # definition above.
                   snap <- restored_ts_hm_annot()

                   group_default <- c(time_col, "Unified_Time", cond_col)
                   group_selected <- if (!is.null(snap) && length(snap$group) > 0) {
                     intersect(snap$group, group_default)
                   } else {
                     character(0)
                   }
                   if (length(group_selected) == 0) group_selected <- cond_col
                   updateSelectizeInput(session,'ts_hm_annot_group',choices=group_default, selected=group_selected)

                   req(DataQCReactive())
                   MetaData <- df_cluster_sample_meta() # DataQCReactive()$MetaData
                   attributes=sort(setdiff(colnames(MetaData), c("sampleid", "Order", "ComparePairs", "replaceable","sizeFactor") ))
                   sample_selected <- if (!is.null(snap) && length(snap$sample) > 0) {
                     intersect(snap$sample, attributes)
                   } else {
                     character(0)
                   }
                   if (length(sample_selected) == 0) sample_selected <- c(time_col, "Unified_Time", cond_col)
                   updateSelectizeInput(session,'ts_hm_annot_sample',choices=attributes, selected = sample_selected)

                   if (!is.null(snap)) restored_ts_hm_annot(NULL)
                 })
                 
#################################################
                 heatmap_sample_plot <- eventReactive(input$ts_plot_sample_heatmap, {
                   withProgress(message = "Processing...", value = 0, {
                     log_info("Time Series: rendering sample heatmap.")
                     req(DataClusterReactive())
                     req(input$ts_hm_annot_sample)
                     DataCluster <-DataClusterReactive()
                     clusters <- DataCluster$clusters
                     
                     gene_clusters <- clusters$df %>% 
                       arrange(cluster)
                     
                     plot_matrix_wide <- df_normalized_log2TPM() %>%
                       as.matrix()
                     
                     plot_matrix_scaled <- t(scale(t(plot_matrix_wide)))

                     if (input$ts_genelable == 'Gene.Name') {
                       ProteinGeneName <- DataQCReactive()$ProteinGeneName 
                       rownames(ProteinGeneName) = ProteinGeneName$UniqueID
                       rownames(plot_matrix_scaled) <- ProteinGeneName[rownames(plot_matrix_scaled),'Gene.Name']
                     }
                     
                     sample_meta <- df_cluster_sample_meta()
                     
                     col_anno=NULL #column annotation
                     if (!is.null(input$ts_hm_annot_sample)) {
                       sel_col=match(input$ts_hm_annot_sample, names(sample_meta))
                       df_annot=sample_meta[, sel_col, drop=FALSE]
                       # After defining df_annot
                       # Convert character columns to factors based on appearance order
                       df_annot[] <- lapply(df_annot, function(x) {
                         if (is.character(x) || is.factor(x)) {
                           return(factor(as.character(x), levels = unique(x)))
                         } else {
                           return(x)
                         }
                       })
                       col_anno=HeatmapAnnotation(df = df_annot)
                       
                       is_num <- sapply(df_annot, is.numeric)
                       num_cols <- names(df_annot)[is_num]
                       cat_cols <- names(df_annot)[!is_num]
                       if (input$ts_hm_annot_color=="Auto-Set by Rand. Seed") { #color with palette selected using random seed (default)
                         discrete_palettes=c("Dark2", "Accent",  "Set2", "Set3", "npg", "nejm", "lancet", "jama", "d3", "uchicago") #removed Set1 for numeric
                         set.seed(input$ts_hm_seed)
                         pal_cat_assigned <- sample(discrete_palettes, length(cat_cols), replace = (length(cat_cols) > length(discrete_palettes)))
                         num_palette <- "Set1"
                       } else if (input$ts_hm_annot_color=="Select Palette") { #color with user selected palettes
                         if (length(input$ts_hm_cat_pal) == 0) log_warn("Time Series heatmap: no color palette selected for category annotations.")
                         validate(need(length(input$ts_hm_cat_pal)>0,message = "Please select color palettes for category annotations"))
                         pal_cat_assigned <- rep(input$ts_hm_cat_pal,  length.out=length(cat_cols) )
                         num_palette <- input$ts_hm_num_pal
                       }
                       color_num_assigned <- get_palette(num_palette , length(num_cols))
                       color_list <- imap(df_annot, function(val, col_name) {
                         idx <- match(col_name, if (is.numeric(val)) num_cols else cat_cols)
                         if (is.numeric(val))
                           hm_m_color(df_annot, col_name, high_col = color_num_assigned[idx])
                         else
                           hm_c_color(df_annot, col_name, pal_cat_assigned[idx])
                       })
                       col_anno=HeatmapAnnotation(df = df_annot, col=color_list, show_annotation_name = TRUE)
                     }
                     
                     # 2. Prepare Row Annotations (Multiple)
                     anno_row <- data.frame(Cluster = as.factor(gene_clusters$cluster))
                     rownames(anno_row) <- gene_clusters$genes
                     
                     cluster_names <- unique(anno_row$Cluster)
                     n_clusters <- length(cluster_names)
                     colors <- get_palette("Set3", n_clusters)
                     names(colors) <- sort(cluster_names)
                     
                     # 4. Create the final color_list
                     color_list <- list(Cluster = colors)    
                     
                     row_anno <- rowAnnotation(Cluster = anno_row$Cluster, col = color_list)
                     
                     # 5. Define the Main Color Palette
                     data_range=quantile(unlist(plot_matrix_scaled), probs=c(0.01, 0.5, 0.99), na.rm=T)
                     col_fun=colorRamp2(data_range, c(input$ts_hm_lowColor,input$ts_hm_midColor, input$ts_hm_highColor) )
                     
                     # 6. Set lay out options
                     cexCol = as.numeric(as.character(input$ts_hm_x_fontsize))
                     cexRow = as.numeric(as.character(input$ts_hm_y_fontsize))
                     
                     labCol = TRUE
                     labRow = TRUE
                     # cat("pheatmap ", dim(data.in), date(), "\n") #debug
                     if (cexRow  == 0 | nrow(plot_matrix_scaled) > input$ts_hm_N_genes) {
                       labRow = FALSE
                       cexRow = 5
                     }
                     if (cexCol == 0) {
                       labCol = FALSE
                       cexCol  = 5
                     }
                     
                     time_col <- time_col()
                     condition_var <- cond_col()
                     
                     # 4. Create the Complex Heatmap
                     Heatmap(
                       plot_matrix_scaled, 
                       name = "normalized expression",
                       row_title=input$ts_hm_row_title,
                       row_title_gp = gpar(fontsize = input$ts_hm_row_title_font_size), 
                       column_title=input$ts_hm_column_title,
                       column_title_gp = gpar(fontsize = input$ts_hm_column_title_font_size),
                       row_names_gp = gpar(fontsize = cexRow),
                       column_names_gp = gpar(fontsize = cexCol),
                       
                       # --- NESTED CLUSTERING ---
                       cluster_rows = (input$ts_hm_gene_cluster == "Yes"),
                       show_row_dend = (input$ts_hm_gene_cluster == "Yes"),           # Shows the dendrogram for each sub-cluster
                       cluster_columns = FALSE, 
                       show_column_dend = FALSE,
                       column_labels = sample_meta$sampleid,
                       
                       # Annotations
                       top_annotation = col_anno,
                       left_annotation = row_anno,
                       
                       # Visuals
                       col = col_fun,
                       show_row_names = labRow,
                       show_column_names = labCol,
                       row_names_side="right",
                       
                       # Split by Cluster (Optional but recommended for clarity)
                       row_split = anno_row$Cluster,      # Existing row split
                       column_split = df_annot[[condition_var]],  # NEW: Splits columns by WT, Het, KO
                       column_gap = unit(4, "mm"),        # Gap between genotype groups
                       row_gap = unit(4, "mm"),           # Gap between clusters
                       border = TRUE
                     )
                   })
                 })
                 
                 output$ts_sample_heatmap <- renderPlot({
                   heatmap_sample_plot()
                 }, height = function() {
                   input$ts_hm_height 
                 })                 

                 # Saves the heatmap last drawn with Plot/Refresh.
                 observeEvent(input$save_ts_sample_heatmap, {
                   if (is.null(input$ts_plot_sample_heatmap) || input$ts_plot_sample_heatmap == 0) {
                     showNotification("Click Plot/Refresh to draw the heatmap first.", type = "warning")
                     return()
                   }
                   log_info("Time Series: sample heatmap saved to output.")
                   saved.num <- length(saved_plots$ts_sample_heatmap) + 1
                   saved_plots$ts_sample_heatmap[[saved.num]] <- heatmap_sample_plot()
                 })

#################################################
                 
                 heatmap_group_plot <- eventReactive(input$ts_plot_group_heatmap, {
                   withProgress(message = "Processing...", value = 0, {
                     log_info("Time Series: rendering group heatmap.")
                     req(DataClusterReactive())
                     req(input$ts_hm_annot_group)
                     DataCluster <-DataClusterReactive()
                     clusters <- DataCluster$clusters
                     
                     
                     gene_clusters <- clusters$df %>% 
                       arrange(cluster)
                     
                     plot_matrix_wide <- df_cluster_normalized_wide() %>%
                       dplyr::select(-cluster) %>%
                       as.matrix()
                     
                     if (input$ts_genelable == 'Gene.Name') {
                       ProteinGeneName <- DataQCReactive()$ProteinGeneName 
                       rownames(ProteinGeneName) = ProteinGeneName$UniqueID
                       rownames(plot_matrix_wide) <- ProteinGeneName[rownames(plot_matrix_wide),'Gene.Name']
                     }
                     
                     group_meta <- df_cluster_group_meta()
                     
                    col_anno=NULL #column annotation
                     if (!is.null(input$ts_hm_annot_group)) {
                       sel_col=match(input$ts_hm_annot_group, names(group_meta))
                       df_annot=group_meta[, sel_col, drop=FALSE]
                       # After defining df_annot
                       # Convert character columns to factors based on appearance order
                       df_annot[] <- lapply(df_annot, function(x) {
                         if (is.character(x)) {
                           return(factor(x, levels = unique(x)))
                         } else {
                           return(x)
                         }
                       })
                       col_anno=HeatmapAnnotation(df = df_annot)
                       
                       is_num <- sapply(df_annot, is.numeric)
                       num_cols <- names(df_annot)[is_num]
                       cat_cols <- names(df_annot)[!is_num]
                       if (input$ts_hm_annot_color=="Auto-Set by Rand. Seed") { #color with palette selected using random seed (default)
                         discrete_palettes=c("Dark2", "Accent",  "Set2", "Set3", "npg", "nejm", "lancet", "jama", "d3", "uchicago") #removed Set1 for numeric
                         set.seed(input$ts_hm_seed)
                         pal_cat_assigned <- sample(discrete_palettes, length(cat_cols), replace = (length(cat_cols) > length(discrete_palettes)))
                         num_palette <- "Set1"
                       } else if (input$ts_hm_annot_color=="Select Palette") { #color with user selected palettes
                         if (length(input$ts_hm_cat_pal) == 0) log_warn("Time Series heatmap: no color palette selected for category annotations.")
                         validate(need(length(input$ts_hm_cat_pal)>0,message = "Please select color palettes for category annotations"))
                         pal_cat_assigned <- rep(input$ts_hm_cat_pal,  length.out=length(cat_cols) )
                         num_palette <- input$ts_hm_num_pal
                       }
                       color_num_assigned <- get_palette(num_palette , length(num_cols))
                       color_list <- imap(df_annot, function(val, col_name) {
                         idx <- match(col_name, if (is.numeric(val)) num_cols else cat_cols)
                         if (is.numeric(val))
                           hm_m_color(df_annot, col_name, high_col = color_num_assigned[idx])
                         else
                           hm_c_color(df_annot, col_name, pal_cat_assigned[idx])
                       })
                       col_anno=HeatmapAnnotation(df = df_annot, col=color_list, show_annotation_name = TRUE)
                     }
                     
                     # 2. Prepare Row Annotations (Multiple)
                     anno_row <- data.frame(Cluster = as.factor(gene_clusters$cluster))
                     rownames(anno_row) <- gene_clusters$genes
                     
                     cluster_names <- unique(anno_row$Cluster)
                     n_clusters <- length(cluster_names)
                     colors <- get_palette("Set3", n_clusters)
                     names(colors) <- sort(cluster_names)
                     
                     # 4. Create the final color_list
                     color_list <- list(Cluster = colors)    
                     
                     row_anno <- rowAnnotation(Cluster = anno_row$Cluster, col = color_list)
                     
                     # 5. Define the Main Color Palette
                     data_range=quantile(unlist(plot_matrix_wide), probs=c(0.01, 0.5, 0.99), na.rm=T)
                     col_fun=colorRamp2(data_range, c(input$ts_hm_lowColor,input$ts_hm_midColor, input$ts_hm_highColor) )
                     
                     # 6. Set lay out options
                     cexCol = as.numeric(as.character(input$ts_hm_x_fontsize))
                     cexRow = as.numeric(as.character(input$ts_hm_y_fontsize))
                     
                     labCol = TRUE
                     labRow = TRUE
                     # cat("pheatmap ", dim(data.in), date(), "\n") #debug
                     if (cexRow  == 0 | nrow(plot_matrix_wide) > input$ts_hm_N_genes) {
                       labRow = FALSE
                       cexRow = 5
                     }
                     if (cexCol == 0) {
                       labCol = FALSE
                       cexCol  = 5
                     }
                     
                     time_col <- time_col()
                     condition_var <- cond_col()
                     
                     # 4. Create the Complex Heatmap
                     Heatmap(
                       plot_matrix_wide, 
                       name = "Z-score",
                       row_title=input$ts_hm_row_title,
                       row_title_gp = gpar(fontsize = input$ts_hm_row_title_font_size), 
                       column_title=input$ts_hm_column_title,
                       column_title_gp = gpar(fontsize = input$ts_hm_column_title_font_size),
                       row_names_gp = gpar(fontsize = cexRow),
                       column_names_gp = gpar(fontsize = cexCol),
                       
                       # --- NESTED CLUSTERING ---
                       cluster_rows = (input$ts_hm_gene_cluster == "Yes"),
                       show_row_dend = (input$ts_hm_gene_cluster == "Yes"),           # Shows the dendrogram for each sub-cluster
                       cluster_columns = FALSE, 
                       show_column_dend = FALSE,
                       column_labels = group_meta[[time_col]],
                       
                       # Annotations
                       top_annotation = col_anno,
                       left_annotation = row_anno,
                       
                       # Visuals
                       col = col_fun,
                       show_row_names = labRow,
                       show_column_names = labCol,
                       row_names_side="right",
                       
                       # Split by Cluster (Optional but recommended for clarity)
                       row_split = anno_row$Cluster,      # Existing row split
                       column_split = df_annot[[condition_var]],  # NEW: Splits columns by WT, Het, KO
                       column_gap = unit(4, "mm"),        # Gap between genotype groups
                       row_gap = unit(4, "mm"),           # Gap between clusters
                       border = TRUE
                     )
                   })
                 })
                 
                 output$ts_group_heatmap <- renderPlot({
                   heatmap_group_plot()
                 }, height = function() {
                   input$ts_hm_height 
                 })                 

                 # Saves the heatmap last drawn with Plot/Refresh.
                 observeEvent(input$save_ts_group_heatmap, {
                   if (is.null(input$ts_plot_group_heatmap) || input$ts_plot_group_heatmap == 0) {
                     showNotification("Click Plot/Refresh to draw the heatmap first.", type = "warning")
                     return()
                   }
                   log_info("Time Series: group heatmap saved to output.")
                   saved.num <- length(saved_plots$ts_group_heatmap) + 1
                   saved_plots$ts_group_heatmap[[saved.num]] <- heatmap_group_plot()
                 })
               })
}



