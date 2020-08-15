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


suppressPackageStartupMessages({
	library(shiny)
	library(shinythemes)
	library(plotly)
	#library(heatmaply)
	library(reshape2)
	library(dplyr)
	library(tidyr)
	library(tibble)
	library(ggplot2)
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
	#library(car)
	library(colourpicker)
	library(VennDiagram)
	#library(pathview)
	library(factoextra)
	#library(cluster)
	#library(NbClust)
	#library(Mfuzz)
	library(openxlsx)
	#library(rmarkdown)
	library(visNetwork)
	#library(networkD3)
	#library(Hmisc)
})


load("db/hgnc.RData")
load("db/kegg.pathways.RData")
load("db/gmtlist.RData")



source("scatter3d.R",local = TRUE)
source("pathviewfun.R",local = TRUE)
source("mfuzzfun.R",local = TRUE)

	
ORAEnrichment <- function(deGenes,universe, gsets, logFC){
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
		ora.stats[j, "p.value"] = phyper(q = totalDEinS- 0.5, m=totalDE, n = n, k = totalSinUniverse, lower.tail = FALSE)
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
