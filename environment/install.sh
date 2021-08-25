#!/bin/bash

~/conda/bin/mamba create -y -n QUICKO -c conda-forge -c bioconda -c anaconda -c r r-shiny r-shinythemes r-shinyjs plotly  r-plotly r-reshape2 r-tidyverse r-gplots r-ggpubr r-gridextra r-ggrepel  r-rcolorbrewer  r-pheatmap r-rgl r-car r-colourpicker r-venndiagram  r-factoextra r-openxlsx r-visnetwork r-cowplot r-circlize bioconductor-complexheatmap bioconductor-interactivecomplexheatmap r-svglite r-shinyjqui r-hmisc r-ggrastr r-ggextra r-networkd3 r-vctrs r-ragg r-textshaping r-biocmanager bioconductor-mfuzz

#R --vanilla -e "source('../.Rprofile'); BiocManager::install('pathview')"
#Manually installl pathview
# BiocManager::install("pathview")
#manually install psych
#manually install BiocManager::install("biomaRt")
