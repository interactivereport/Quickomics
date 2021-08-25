source('.Rprofile')
###########################################################################################################
## Proteomics Visualization R Shiny App
##
##This software belongs to Biogen Inc. All right reserved.
##
##@file: server.R
##@Developer : Benbo Gao (benbo.gao@Biogen.com)
##@Date : 5/16/2018
##@version 1.0
###########################################################################################################
options(shiny.host = '0.0.0.0')
options(shiny.port = 8023)

function(input, output, session) {
source("inputdata.R",local = TRUE)
source("process_uploaded_files.R",local = TRUE)
source("qcplot.R",local=TRUE)
source("volcano.R",local = TRUE)
source("heatmap.R",local = TRUE)
source("barboxplot.R",local = TRUE)
source("venn.R",local = TRUE)
source("geneset.R",local = TRUE)
source("pattern.R",local = TRUE) 
source("vennprojects.R",local = TRUE)
source("network.R",local = TRUE)
source("help.R",local = TRUE)
source("output.R",local = TRUE)
source("scurve.R",local = TRUE)  
}

