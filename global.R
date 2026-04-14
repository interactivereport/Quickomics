#source("renv.restore.R")

options(stringsAsFactors=F)
options(ggrepel.max.overlaps = Inf)

suppressPackageStartupMessages({
	library(shiny)
	library(shinythemes)
	library(shinyjs)
	library(plotly)
	#library(heatmaply)
	library(reshape2)
	#library(dplyr)
	#library(tidyr)
	#library(tibble)
	#library(stringr)
	#library(ggplot2)
	library(tidyverse)
	library(gplots)
	library(ggpubr)
	library(gridExtra)
	library(ggrepel)
	#library(dendextend)
	library(DT)
	#library(data.table)
	library(RColorBrewer)
	library(pheatmap)
  options(rgl.useNULL = TRUE)
	library(rgl)
	library(car)
	library(colourpicker)
	library(VennDiagram)
	#library(pathview)
	library(factoextra)
	#library(cluster)
	#library(NbClust)
	library(Mfuzz)
	library(openxlsx)
	#library(rmarkdown)
	library(visNetwork)
	library(cowplot)
	library(circlize)
	library(ComplexHeatmap)
	library(svglite)
	library(pathview)
	library(shinyjqui)
	library(ggrastr)
	library(biomaRt)
	library(networkD3)
	library(Hmisc)
	library(stringi)
	library(plyr)
	library(ggExtra)
  library(png)
  library(psych)
  library(broom)
  library(fgsea)
  #library(DEGreport)
  library(rclipboard)
  library(WGCNA)
  library(DEGreport)
  library(gprofiler2)
  #library(cmapR)
})



#load("db/hgnc.RData")
#load("db/kegg.pathways.RData")
#load("db/gmtlist.RData")
#load("db/mouse_rat_genes_map2_human.RData")
homologs=readRDS("db/Homologs.rds") #cross species gene symbol mapping file from Ensembl Gene 110
source("PC_Covariates.R")
source("write_gct.R")
homolog_mapping<-function(genelist, species1, species2, homologs) {
  if (species2=="human") {
    genelist2=toupper(genelist)
  } else {
    genelist2=str_to_title(genelist)
  }
  if (tolower(species1) %in% c("human", "mouse", "rat")) {
    species1=tolower(species1)
    lookup<-homologs%>%filter(Species1==species1, Species2==species2)
    df<-data.frame(genelist, genelist2)%>%left_join(lookup%>%transmute(genelist=symbol1, mapped_symbol=symbol2))%>%
      mutate(mapped_symbol=ifelse(is.na(mapped_symbol), genelist2, mapped_symbol))
    genelist2<-df$mapped_symbol
  }
  return(genelist2)
}

GeneFilter <- function(results_long, test_sel, p_sel, direction, pvalcut, FCcut, sel_label) {
  results_long_filtered <- results_long %>%
    dplyr::mutate_if(is.factor, as.character) %>%
    dplyr::filter(if (!is.na(test_sel)) {test == test_sel} else TRUE) %>%
    dplyr::filter(if (p_sel=="Padj") {Adj.P.Value < pvalcut} else {P.Value < pvalcut}) %>%
    dplyr::filter(if (direction=="Up") {logFC >= FCcut} else if (direction=="Down") {logFC <= -FCcut} else {abs(logFC) >= FCcut}) %>%
    dplyr::arrange(desc(abs(logFC)))  %>%
    dplyr::mutate(labelid = !!sym(sel_label)) %>%
    dplyr::filter(!is.na(labelid) & labelid != "") %>%
    as.data.frame()
  return(results_long_filtered)
}


ProcessUploadGeneList <- function(gene_list) {
  gene_list <- unlist(strsplit(gene_list, "[,\n]"))
  gene_list <- gsub(" ", "", gene_list, fixed = TRUE)
  gene_list <- unique(gene_list[gene_list != ""])
  return(gene_list)
}

GetProteinGeneNames <- function(species) {
  load("db/ProteinGeneName.RData")
  return(ProteinGeneName[[species]])
}

LoadedData <- reactiveValues()
GetGeneSetNames <- function() {
  if(!("gmtlist" %in% names(LoadedData))){
    load("db/gmtlist.RData")
    LoadedData[["gmtlist"]]  <-  gmtlist
  } else {
    gmtlist <- LoadedData[["gmtlist"]]
  }
  
  if(!("kegg.pathways" %in% names(LoadedData))){
    load("db/kegg.pathways.RData")
    LoadedData[["kegg.pathways"]]  <-  kegg.pathways
  } else{
    kegg.pathways <- LoadedData[["kegg.pathways"]]
  }
  
  genesetnames <- c()
  for (setname in names(gmtlist)) {
    genesetnames <- c(genesetnames, names(gmtlist[[setname]]))
  }
  genesetnames <- c(names(kegg.pathways$human$kg.sets), genesetnames)
  return(genesetnames)
}

GetGenesFromGeneSet <- function(sel_geneset) {
  if(!("hgnc" %in% names(LoadedData))){
    load("db/hgnc.RData")
    LoadedData[["hgnc"]]  <-  hgnc
  } else {
    hgnc <- LoadedData[["hgnc"]]
  }
  kegg.pathways <- LoadedData[["kegg.pathways"]]
  gmtlist <- LoadedData[["gmtlist"]]
  
  if (sel_geneset %in% names(kegg.pathways$human$kg.sets)) {
    geneset_genes <- kegg.pathways$human$kg.sets[[sel_geneset]]
  }	else {
    for (setname in names(gmtlist)) {
      if (sel_geneset %in% names(gmtlist[[setname]])) {
        geneset_genes <- gmtlist[[setname]][[sel_geneset]]
        break
      }
    }
  }
  geneset_genenames <- hgnc %>%
    dplyr::filter(entrez_id %in% geneset_genes) %>%
    dplyr::pull(symbol)
  return(geneset_genenames)
}

#' Create a GCT object with strict ID matching
#' @param exp_data A data frame or matrix of expression values
#' @param row_meta A data frame containing gene/protein annotations
#' @param col_meta A data frame containing sample annotations
#' @return A GCT object
create_gct_object <- function(exp_data, row_meta, col_meta) {
 # library(cmapR)
  # 1. Prepare the Expression Matrix
  my_matrix <- as.matrix(exp_data)
  my_matrix[is.na(my_matrix)] <- 0
  
  # 2. Internal function to find, rename, AND reorder metadata
  align_metadata <- function(meta_df, target_names, type = "row") {
    meta_df <- as.data.frame(meta_df)
    match_idx <- NULL
    
    # Find the column that contains all the IDs (regardless of order)
    for (i in seq_len(ncol(meta_df))) {
      column_values <- as.character(meta_df[[i]])
      # Check if all matrix IDs exist in this metadata column
      if (all(target_names %in% column_values)) {
        match_idx <- i
        break
      }
    }
    if (is.null(match_idx)) {
      stop(paste("Mapping Failed: No column in", type, 
                 "metadata contains all IDs found in the expression matrix."))
    }
    id_col_name <- colnames(meta_df)[match_idx]
    message(paste0("Match found: Using '", id_col_name, "' to align ", type, " metadata."))
    
    # Rename the matching column to 'id' (removing existing 'id' column if it's the wrong one)
    if ("id" %in% colnames(meta_df) && id_col_name != "id") {
      meta_df$id <- NULL
      match_idx <- which(colnames(meta_df) == id_col_name)
    }
    colnames(meta_df)[match_idx] <- "id"
    
    # CRITICAL STEP: Reorder metadata rows to match the matrix order
    reorder_idx <- match(target_names, meta_df$id)
    meta_df <- meta_df[reorder_idx, , drop = FALSE]
    meta_df <- meta_df[, c("id", setdiff(colnames(meta_df), "id")), drop = FALSE] 
    
    return(meta_df)
  }  
  
  # 3. Process Row Metadata
  # Keep only columns where the number of unique values is greater than 1
  row_meta <- row_meta[sapply(row_meta, function(x) length(unique(x)) > 1)]
  row_metadata <- align_metadata(row_meta, rownames(my_matrix), "row")
  
  # 4. Process Column Metadata
  # First, clean the specific columns you identified in your environment
  cols_to_remove <- intersect(c("Order", "ComparePairs"), colnames(col_meta))
  col_meta[] <- lapply(col_meta, as.character)
  col_metadata <- col_meta
  if (length(cols_to_remove) > 0) {
    col_metadata <- col_metadata[, !colnames(col_metadata) %in% cols_to_remove]
  }
  col_metadata <- align_metadata(col_metadata, colnames(my_matrix), "column")
  
  # 5. Assemble the GCT object
  my_gct <- new("GCT", 
                mat = my_matrix, 
                rdesc = row_metadata, 
                cdesc = col_metadata)
  return(my_gct)
}

#functions to assign heatmap colors
hm_m_color<-function(df, var, low_col="white", high_col=color, min=0, max=0.999) { #numeric annotations
  data<-df[var]%>%unlist%>%unname
  q<-quantile(data, c(min, max) )
  col_fun = colorRamp2(c(q[1], q[2]), c(low_col, high_col))
  return(col_fun)
}
hm_c_color<-function(df, var, colPal, sort=T) { #category annotation
  cat<-df[var]%>%unlist%>%unname%>%unique
  if (sort) {cat<-sort(cat) }
  Nc=length(cat)
  colorSet=get_palette(colPal, Nc)
  names(colorSet)=cat
  return(colorSet)
}


options(shiny.maxRequestSize = 40*1024^2)  #upload files up to 30 Mb


mycss <- "select ~ .selectize-control .selectize-input {
        max-height: 200px;
        overflow-y: auto;
        }
    .shiny-notification{
    position: fixed;
    top: 33%;
    left: 33%;
    right: 33%;
    }"

saved_projects = read.csv("data/saved_projects.csv")
projects = saved_projects$ProjectID
names(projects) = saved_projects$ShortNames

pub_projects=NULL
if (file.exists("unlisted/projects_public.csv")) { #load selected projects from unlisted folder for Venn across projects
  Pinfo=read.csv("unlisted/projects_public.csv")
  pub_projects=Pinfo$ProjectID
  names(pub_projects)=Pinfo$ShortName
}

html_geneset0 =  '
<script>
var GENESET_DEFAULT_TABLE = "";
var GENESET_DEFAULT_SPECIES = "human";
var GENESET_DEFAULT_MIN = 10;
var GENESET_DEFAULT_MAX = 500;
</script>

<div id="div_geneset3" class="div_geneset">
<div class="my-3 dropdown btn-group">
<input style="display: none;" id="Geneset_Name1" name="Geneset_Name1" placeholder="Start typing to enter or select a geneset" type="text" class="form-control form-control-sm geneset_name" />
<a href="Javascript: void(0);" class="btn btn-primary btn_browse_geneset"><i class="fas fa-search"></i> Select Geneset</a>
<div id="Dropdown1" class="geneset_dropdown dropdown-menu"></div>
</div>

<input class="geneset_id" type="hidden" id="Geneset_ID1" name="Geneset_ID1" />
<label class="control-label" for="geneset_list">List of genes to label (UniqueID, Gene.Name or Protein.ID)</label>
<textarea id="geneset_list" name="geneset_list" class="form-control shiny-bound-input geneset_genes my-3" rows="6" cols="5"></textarea>
</div>
'

html_geneset_hm0 =str_replace_all(html_geneset0, "geneset_list", "geneset_list_hm")
html_geneset_exp0 =str_replace_all(html_geneset0, "geneset_list", "geneset_list_exp")



footer_text = '
<script>
var GENESET_ACTION_URL = "//bxngs.com/bxomics/api/genesets3.php";
var MY_SECRET_ID = /PHPSESSID=([^;]+)/i.test(document.cookie) ? RegExp.$1 : false;
</script>
<link href="//bxngs.com/bxomics/api/genesets3.css" rel="stylesheet">
<script src="//bxngs.com/bxomics/api/genesets3.js"></script>
<hr>
<div align="center" style="font-size:11px">QuickOmics ver3.0 Developed by:
Benbo Gao, Xinmin Zhang and Baohong Zhang<br><a href="https://github.com/interactivereport/Quickomics/">More information at GitHub</a> | <a href="https://interactivereport.github.io/Quickomics/tutorial/docs/introduction.html">Tutorial</a>
</div>
'

config=NULL
server_dir=NULL
test_dir=NULL
gmt_file_info=NULL
system_info=NULL
public_dataset=TRUE
if (file.exists("config.csv")) { #load optional configuration file
  config=read_csv("config.csv")
  N=match("server_dir", config$category)
  if (!is.na(N)) {server_dir=config$value[N]}
  N=match("test_dir", config$category)
  if (!is.na(N)) {test_dir=config$value[N]}
  N=match("geneset_api", config$category) #if using internal api for geneset
  if (!is.na(N)) {footer_text=str_replace_all(footer_text, "//bxngs.com/bxomics/api", config$value[N])}
  N=match("gmt_file_info", config$category)
  if (!is.na(N)) {gmt_file_info=config$value[N]}
  N=match("system_info", config$category)
  if (!is.na(N)) {system_info=config$value[N]}
  N=match("restricted_dataset", config$category)
  if (!is.na(N)) {public_dataset=ifelse(config$value[N] == "yes", FALSE, TRUE)}
  #browser() #debug
}
