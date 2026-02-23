#!/usr/bin/env Rscript

# ================================
# Standalone WGCNA Analysis Script
# ================================

config_check <- function(config) {
  errors <- c()
  # -----------------------------
  # Helper: check required field
  # -----------------------------
  check_required <- function(value, name) {
    if (is.null(value) || value == "" || is.na(value)) {
      errors <<- c(errors, paste0("Missing required parameter: ", name))
      return(FALSE)
    }
    TRUE
  }
  crit <- config$critical
  # project_name
  if (check_required(crit$project_name, "critical$project_name")) {
    if (!is.character(crit$project_name))
      errors <- c(errors, "critical$project_name must be a string.")
  }
  # expression_datafile
  if (check_required(crit$expression_datafile, "critical$expression_datafile")) {
    if (!file.exists(crit$expression_datafile)) {
      errors <- c(errors, paste0(
        "File not found: critical$expression_datafile = ",
        crit$expression_datafile
      ))
    }
  }
  # outout_dir
  if (check_required(crit$outout_dir, "critical$outout_dir")) {
    dir.create(crit$outout_dir, showWarnings = FALSE)
  }
  basic <- config$basic
  if (!is.null(basic$top_variable_genes)) {
    if (!is.numeric(basic$top_variable_genes)) {
      errors <- c(errors, "basic$top_variable_genes must be numeric.")
    }
  }
  # ============================================================
  # FINAL RESULT
  # ============================================================
  if (length(errors) == 0) {
    return("Config OK")
  } else {
    return(errors)
  }
}

load_expression_data <- function(data_file) {
  ext <- tolower(tools::file_ext(data_file))
  if (ext == "rdata") {
    env <- new.env()
    load(data_file, envir = env)
   if (!"data_wide" %in% ls(env)) {
      stop("RData file does not contain an object named 'data_wide'.")
    }
    return(env$data_wide)
  } else if (ext %in% c("csv", "txt")) {
    df <- data.table::fread(data_file)
    df <- as.data.frame(df)
    rownames(df) <- df[[1]]
    df[[1]] <- NULL 
    return(df)
  }
  stop("Unsupported file extension: ", ext,
       ". Expected .RData, .csv, or .txt")
}

clean_expression_data <- function(df) {
  df <- as.data.frame(df)
  numeric_cols <- sapply(df, is.numeric)
  df <- df[, numeric_cols, drop = FALSE]
  is_inf <- is.infinite(as.matrix(df))
  if (any(is_inf)) {
    message("Converting ", sum(is_inf), " infinite values to NA")
    df[is_inf] <- NA
  }
  keep <- complete.cases(df)
  df <- df[keep, , drop = FALSE]
  return(df)
}

prepare_soft_threshold_table <- function(sft, picked_power) {
  # Extract the fit indices table
  df <- sft$fitIndices
  # Add derived columns
  df$Power  <- df[, 1]
  df$SFT_R2 <- -sign(df[, 3]) * df[, 2]
  df$MeanK  <- df[, 5]
  df$Picked <- ifelse(df$Power == picked_power, "picked", "")
  # Round numeric columns except Power
  num_cols <- sapply(df, is.numeric)
  num_cols["Power"] <- FALSE
  df[num_cols] <- lapply(df[num_cols], function(x) round(x, 2))
  return(df)
}


save_qc_plots_pdf <- function(df, output_file) {
  library(ggplot2)
  library(patchwork)
  
  # Base themes with extra bottom margin to avoid cropping
  base_theme <- theme_minimal() +
    theme(
      plot.margin = margin(t = 20, r = 25, b = 20, l = 25)
    )
  
  # Plot 1: Scale Independence
  p1 <- ggplot(df, aes(x = Power, y = SFT_R2)) +
    geom_point(color = "red") +
    geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
    labs(
      title = "Scale Independence",
      x = "Soft Threshold (power)",
      y = "Signed R^2"
    ) +
    base_theme
  
  # Plot 2: Mean Connectivity
  p2 <- ggplot(df, aes(x = Power, y = MeanK)) +
    geom_point(color = "red") +
    labs(
      title = "Mean Connectivity",
      x = "Soft Threshold (power)",
      y = "Mean Connectivity"
    ) +
    base_theme
  
  # Practical decision rules text (plain, left-aligned)
  rules_text <- paste(
    "Practical decision rules:",
    "",
    "- If R² rises and then plateaus, pick the lowest beta at the plateau.",
    "- If R² rises slowly but connectivity is still good, pick beta = 6-10 (signed).",
    "- If R² is flat and connectivity collapses at high beta, pick beta = 6-8.",
    "- If R² keeps rising up to the max tested power, expand the tested range (1-30).",
    "- Ensure mean connectivity is not too low. It should > ~ 20 (for typical datasets).",
    "- If estimated power is NA, use beta = 6-12 for signed networks.",
    sep = "\n"
  )
  
  combined <- (p1 + p2) +
    plot_annotation(
      caption = rules_text,
      theme = theme(
        plot.caption = element_text(
          hjust = 0,          # left align
          vjust = 1,
          size  = 12,
          margin = margin(t = 20, r = 25, b = 20, l = 25)
        ),
        plot.margin = margin(t = 20, r = 25, b = 20, l = 25)
      )
    )
  
  # Larger device to avoid cropping
  pdf(output_file, width = 12, height = 8)
  print(combined)
  dev.off()
}


generate_soft_threshold_pdf <- function(sft,
                                        powers,
                                        picked_power,
                                        output_file = "soft_threshold_report.pdf") {
  
  # ---- Dependencies ----
  stopifnot(
    requireNamespace("ggplot2"),
    requireNamespace("patchwork"),
    requireNamespace("rmarkdown"),
    requireNamespace("knitr"),
    requireNamespace("dplyr")
  )
  
  # ---- Prepare data ----
  table_df <- prepare_soft_threshold_table(sft, picked_power)

  # ---- Plots ----
  library(ggplot2)
  library(patchwork)
  
  p1 <- ggplot(table_df, aes(x = Power, y = SFT_R2, label = powers)) +
    geom_point(color = "red") +
    geom_text(vjust = -0.5, color = "red", size = 3) +
    geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
    labs(
      x = "Soft Threshold (power)",
      y = "Scale Free Topology Model Fit, signed R^2",
      title = "Scale Independence"
    ) +
    theme_minimal(base_size = 11)
  
  p2 <- ggplot(table_df, aes(x = Power, y = MeanK, label = powers)) +
    geom_point(color = "red") +
    geom_text(vjust = -0.5, color = "red", size = 3) +
    labs(
      x = "Soft Threshold (power)",
      y = "Mean Connectivity",
      title = "Mean Connectivity"
    ) +
    theme_minimal(base_size = 11)
  
  combined_plot <- p1 + p2
  
  # ---- Text blocks ----
  explanation_text_top <- "
In WGCNA, soft-thresholding power (beta) controls how strongly you emphasize large correlations
and suppress small ones when building the adjacency matrix.

The table below shows the result for all beta candidates tested.
**The highlighted row corresponds to the selected beta.**
"
  
  explanation_text_bottom <- "
**Practical decision rules:**

- If R² rises and then plateaus, pick the lowest beta at the plateau.  
- If R² rises slowly but connectivity is still good, pick beta = 6–10 (signed).  
- If R² is flat and connectivity collapses at high beta, pick beta = 6–8.  
- If R² keeps rising up to the max tested power, expand the tested range (1–30).
- Ensure mean connectivity is not too low. It should > ~ 20 (for typical datasets).
- If estimated power is NA, use beta = 6–12 for signed networks.  
"
  
  # ---- R Markdown template ----
  rmd_template <- '
---
title: "Soft Thresholding Power Report"
output:
  pdf_document:
    latex_engine: xelatex
geometry: margin=1in
fontsize: 11pt
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE, message = FALSE, warning = FALSE)
library(dplyr)
library(knitr)
```

# Soft-thresholding power (beta) explanation

`r explanation_text_top`

```{r table}
table_df %>%
  kable(align = "r")
```
# 
# \\newpage

# Diagnostic plots

```{r plots, fig.width=9, fig.height=4}
combined_plot
```

# Practical decision rules

`r explanation_text_bottom`
'
  
  # ---- Write temp Rmd and render ----
  tmp_rmd <- tempfile(fileext = ".Rmd")
  
  # Environment for Rmd execution
  env <- new.env(parent = globalenv())
  env$table_df <- table_df
  env$combined_plot <- combined_plot
  env$explanation_text_top <- explanation_text_top
  env$explanation_text_bottom <- explanation_text_bottom
  
  writeLines(rmd_template, con = tmp_rmd)
  
  rmarkdown::render(
    input = tmp_rmd,
    output_file = normalizePath(output_file, mustWork = FALSE),
    envir = env,
    quiet = TRUE
  )
  
  message("PDF generated: ", normalizePath(output_file, mustWork = FALSE))
}


args <- commandArgs(trailingOnly = TRUE)

# -----------------------------------------------------
# Help information and copy config.yml template file.
# -----------------------------------------------------
if (length(args) < 1) {
  cat("
  Usage: Rscript run_wgcna.R <config.yml>
  A template config.yml file has been copied to the current path.
  \n")
  args <- commandArgs(trailingOnly = FALSE)
  script_path <- dirname(sub("--file=", "", args[grep("--file=", args)]))
  file.copy(file.path(script_path, 'run_wgcna_config_template.yml'), "./config.yml")
  quit(save = "no")
}

# ---- config file check ---- 
suppressMessages({
  library(WGCNA)
  library(dplyr)
  library(tibble)
  library(yaml)
  library(data.table)
})
args <- commandArgs(trailingOnly = TRUE)

config <- read_yaml(args[1])

check_result <- config_check(config)

if (identical(check_result, "Config OK")) {
  message("Configuration validated successfully.")
} else {
  message("Configuration errors:")
  print(check_result)
}

# ---- read in data ---- 
project_name <- config$critical$project_name
data_file    <- config$critical$expression_datafile
out_dir      <- config$critical$outout_dir

data_wide <- load_expression_data(data_file)
data_wide <- clean_expression_data(data_wide)

diff <- apply(data_wide, 1, sd, na.rm = TRUE)/(rowMeans(data_wide) + median(rowMeans(data_wide)))
data_wide=data_wide[order(diff, decreasing=TRUE), ]                            

number_top_genes <- as.numeric(config$basic$top_variable_genes)

if (0 < number_top_genes && number_top_genes < 1) {
  data_wide <- data_wide[1:floor(nrow(data_wide)*number_top_genes), ]
  dataExpr <- as.data.frame(t(data_wide))
} else if (1 <= number_top_genes && number_top_genes < nrow(data_wide)) {
  dataExpr <- data_wide %>%
    dplyr::select(where(is.numeric)) %>%
    { as.data.frame(t(.)) } %>%
    { .[, 1L:number_top_genes]}
} else if (number_top_genes >= nrow(data_wide)) {
  dataExpr <- data_wide %>%
    dplyr::select(where(is.numeric)) %>%
    { as.data.frame(t(.)) }
}

if (as.numeric(config$critical$WGCNAThreads) == 0) {
  enableWGCNAThreads()
} else {
  enableWGCNAThreads(as.numeric(config$critical$WGCNAThreads))
}

# ---- Pick soft-thresholding power ---- 
t2 <- Sys.time()

if (!is.null(config$basic$power)) {
  powers <- c(1L:10L, seq(12L, 20L, by = 2L))
} else {
  power_limit <- as.integer(config$advanced$power_limit) 
  if (power_limit <= 10L) {
    powers <- c(1L:power_limit)
  } else {
    powers <- c(1L:10L, seq(12L, power_limit, by = 2L))
  }
}

cor <- stats::cor
sft <- WGCNA::pickSoftThreshold(
  dataExpr,
  dataIsExpr = TRUE,
  powerVector = powers,
  corFnc = cor,
  corOptions = list(use = "p"),
  networkType = config$advanced$networkType
)

if (is.null(config$basic$power)) {
  if (!is.na(sft$powerEstimate)) {
    message("**** Pick power from sft$powerEstimate ****")
    picked_power <- sft$powerEstimate
  } else {
    message("**** Pick power of default value ****")
    picked_power <- 6L
  }
  t3 <- Sys.time()
  cat(paste0("Computing softpower: ", round(difftime(t3, t2, units = "mins"), 2), " min\n"))
} else {
  picked_power <- config$basic$power
  message("**** Use pre-assigned power ****")
  t3 <- Sys.time()
}

cor <- WGCNA::cor 
netwk <- blockwiseModules(
  dataExpr,
  power = picked_power,
  networkType = config$advanced$networkType,
  deepSplit = config$advanced$deepSplit,
  pamRespectsDendro = config$advanced$pamRespectsDendro,
  minModuleSize = config$advanced$minModuleSize,
  maxBlockSize = ncol(dataExpr),
  reassignThreshold = 0,
  mergeCutHeight = config$advanced$mergeCutHeight,
  saveTOMs = config$basic$saveTOMs_files,
  saveTOMFileBase = config$basic$saveTOMFileBase,
  numericLabels = config$advanced$numericLabels,
  loadTOM = config$basic$loadTOM,
  verbose = 3L
)

t4 <- Sys.time()
cat(paste0("Run WGCNA: ", round(difftime(t4, t3, units = "mins"), 2), " min\n"))
cor <- stats::cor 
# -----------------------------
# Step 5: Save results
# -----------------------------
QC_type <- config$critical$QC_report_type
message("QC report type: ", QC_type)

# Always save RData files (unchanged)
# File 1: load_<project_name>.RData (contains dataExpr and picked_power)
save(dataExpr, picked_power, sft, file = file.path(out_dir, paste0("load_", project_name, ".RData")))
# File 2: wgcna_<project_name>.RData (contains netwk)
save(netwk, file = file.path(out_dir, paste0("wgcna_", project_name, ".RData")))

# QC branching
if (QC_type == "None") {
  message("Skipping QC output (QC_report_type = None).")
} else if (QC_type == "Files") {
  message("Generating QC files (CSV + combined PDF plots).")
  # 1. Save soft-threshold table
  qc_table_file <- file.path(out_dir, paste0(project_name, "_soft_threshold_table.csv"))
  # ---- Prepare data ----
  table_df <- prepare_soft_threshold_table(sft, picked_power)
  write.csv(table_df, qc_table_file, row.names = FALSE)
  
  # 2. Save combined QC plots
  qc_plot_file <- file.path(out_dir, paste0(project_name, "_soft_threshold_QC_plots.pdf"))
  # df <- sft$fitIndices
  save_qc_plots_pdf(table_df, qc_plot_file)
} else if (QC_type == "PDF") {
  message("Rendering full QC PDF report via Rmd.")
  generate_soft_threshold_pdf(
    sft = sft,
    powers = powers,
    picked_power = picked_power,
    output_file = file.path(out_dir, paste0(project_name, "_soft_threshold_QC_report.pdf"))
  )
} 
cat("✅ WGCNA analysis complete. Results saved as load_", project_name, ".RData and wgcna_", project_name, ".RData\n", sep = "")
