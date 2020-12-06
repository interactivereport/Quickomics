###########################################################################################################
## Proteomics Visualization R Shiny App
##
##This software belongs to Biogen Inc. All right reserved.
##
##@file: global.R
##@Developer : Benbo Gao (benbo.gao@Biogen.com)
##@Date : 5/16/2018
##@version 1.0
###########################################################################################################
options(stringsAsFactors=F)

suppressPackageStartupMessages({
	library(shiny)
	library(shinythemes)
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
	#library(networkD3)
	#library(Hmisc)
})


load("db/hgnc.RData")
load("db/kegg.pathways.RData")
load("db/gmtlist.RData")
load("db/mouse_rat_genes_map2_human.RData")


#source("scatter3d.R",local = TRUE)
#source("pathviewfun.R",local = TRUE)
#source("mfuzzfun.R",local = TRUE)

	
ORAEnrichment <- function(deGenes,universe, gsets, logFC, Dir="Both"){
  deGenes = deGenes[which(deGenes %in% universe)]
  tmp = rep(NA, length(gsets))
  ora.stats = data.frame(p.value=tmp, p.adj = tmp, DeGeneNum=tmp,UpGene= tmp, DownGene=tmp, SetNum = tmp)
  totalDE = length(deGenes)
  n = length(universe) - totalDE
  
  for (j in 1:length(gsets)){
    gset = gsets[[j]]
    DEinS = intersect(gset, deGenes)
    logFCinS = logFC[DEinS]
    totalDEinS = length(intersect(gset, deGenes))
    totalSinUniverse = length(intersect(gset, universe))
    
    N_q=totalDEinS- 0.5
    if (Dir=="Up") {N_q=length(logFCinS[logFCinS > 0])-0.5 
		} else if (Dir=="Down") {N_q=length(logFCinS[logFCinS < 0])-0.5}
		
		ora.stats[j, "p.value"] = phyper(q = N_q, m=totalDE, n = n, k = totalSinUniverse, lower.tail = FALSE)
		ora.stats[j, "DeGeneNum"] = totalDEinS
		ora.stats[j, "SetNum"] = length(gset)
		ora.stats[j, "UpGene"] = length(logFCinS[logFCinS > 0])
		ora.stats[j, "DownGene"] = length(logFCinS[logFCinS < 0])

	}
	ora.stats[, "p.adj"] = p.adjust(ora.stats[, "p.value"], method = "BH")
	row.names(ora.stats) = names(gsets)

	ora.stats = ora.stats[order(ora.stats[,"p.value"]),]
	ora.stats = cbind(Rank=seq(1, nrow(ora.stats)),ora.stats)
	return(ora.stats)
}

options(shiny.maxRequestSize = 40*1024^2)  #upload files up to 30 Mb


mycss <- "select ~ .selectize-control .selectize-input {
        max-height: 100px;
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
<link rel="stylesheet" type="text/css" href="//bxngs.com/bxomics/api/datatables/datatables.min.css"/>
<script type="text/javascript" src="//bxngs.com/bxomics/api/datatables/datatables.min.js"></script>
<script type="text/javascript" src="url-input-binding.js"></script>

<script>
var GENESET_ACTION_URL = "//bxngs.com/bxomics/api/genesets3.php";
var MY_SECRET_ID = /PHPSESSID=([^;]+)/i.test(document.cookie) ? RegExp.$1 : false;
</script>
<link href="//bxngs.com/bxomics/api/genesets3.css" rel="stylesheet">
<script src="//bxngs.com/bxomics/api/genesets3.js"></script>
<hr>
<div align="center" style="font-size:11px">QuickOmics ver1.0 Developed by:
Benbo Gao, Xinmin Zhang and Baohong Zhang<br><a href="https://github.com/interactivereport/Quickomics/">More information at GitHub</a>
</div>
'
