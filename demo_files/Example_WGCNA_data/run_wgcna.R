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
  df <- sft$fitIndices
  df$Power  <- df[, 1]
  df$SFT_R2 <- -sign(df[, 3]) * df[, 2]
  df$MeanK  <- df[, 5]
  
  table_df <- df
  table_df$Picked <- ifelse(table_df$Power == picked_power, "picked", "")
  
  # Format numeric columns to 2 decimals (except Power)
  num_cols <- sapply(table_df, is.numeric)
  num_cols["Power"] <- FALSE
  table_df[num_cols] <- lapply(table_df[num_cols], function(x) round(x, 2))
  
  # ---- Plots ----
  library(ggplot2)
  library(patchwork)
  
  p1 <- ggplot(df, aes(x = Power, y = SFT_R2, label = powers)) +
    geom_point(color = "red") +
    geom_text(vjust = -0.5, color = "red", size = 3) +
    geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
    labs(
      x = "Soft Threshold (power)",
      y = "Scale Free Topology Model Fit, signed R^2",
      title = "Scale Independence"
    ) +
    theme_minimal(base_size = 11)
  
  p2 <- ggplot(df, aes(x = Power, y = MeanK, label = powers)) +
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

enableWGCNAThreads()

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

generate_soft_threshold_pdf(
  sft = sft,
  powers = powers,
  picked_power = picked_power,
  output_file = file.path(out_dir, "soft_threshold_QC_report.pdf")
)

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
# File 1: load_<project_name>.RData (contains dataExpr and picked_power)
save(dataExpr, picked_power, sft, file = file.path(out_dir, paste0("load_", project_name, ".RData")))

# File 2: wgcna_<project_name>.RData (contains netwk)
save(netwk, file = file.path(out_dir, paste0("wgcna_", project_name, ".RData")))

cat("✅ WGCNA analysis complete. Results saved as load_", project_name, ".RData and wgcna_", project_name, ".RData\n", sep = "")
