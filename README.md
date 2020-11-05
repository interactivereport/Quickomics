# Quickomics
Smart query tool to explore omics data in an intuitive, interactive and informative manner.

# Requirments
```
shiny >= v1.4.0.2 
```
# Installation
1) Instal the following R packages:
```
cran_packages=c("shiny", "shinythemes", "plotly", "reshape2", "tidyverse", "gplots", "ggpubr", "gridExtra", "ggrepel",
"RcolorBrewer", "pheatmap", "rgl", "car", "colourpicker", "VennDiagram", "factoextra",  "openxlsx", "visNetwork", "cowplot", "circlize", 
"ComplexHeatmap", "svglite", "shinyjgui", "Hmisc")
#Note: Hmisc is not required to run the Shiny App, but is needed to prepare network data from expression matrix.
install.packages(cran_packages, repos="http://cran.r-project.org/")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("Mfuzz"))
```
2) Clone Quickomics GitHub repository
e.g.
```
git clone https://github.com/YOUR-USERNAME/YOUR-REPOSITORY
```
see more at https://docs.github.com/en/free-pro-team@latest/github/creating-cloning-and-archiving-repositories/cloning-a-repository

3) Launch the R Shiny App
Check the follwing web links on various options to launch the app. 
https://shiny.rstudio.com/articles/running.html
https://shiny.rstudio.com/deploy/

# Prepare your own data
Two files are needed for each data set, one contains the main data, the other contains network information. For the pre-loaded datasets, the main data files are located in the data/ folder, while the network files are located at networkdata/ folder. 
To prepare the main data files, you need to get the following
1) 
