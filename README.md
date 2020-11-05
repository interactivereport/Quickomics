# Quickomics
Smart query tool to explore omics data in an intuitive, interactive and informative manner.

# Requirments
```
shiny >= v1.4.0.2 
```
# Installation
### 1) Instal the following R packages:
```R
cran_packages=c("shiny", "shinythemes", "plotly", "reshape2", "tidyverse", "gplots", "ggpubr", "gridExtra", "ggrepel",
"RcolorBrewer", "pheatmap", "rgl", "car", "colourpicker", "VennDiagram", "factoextra",  "openxlsx", "visNetwork",
"cowplot", "circlize", "ComplexHeatmap", "svglite", "shinyjgui", "Hmisc")
#Note: Hmisc is not required to run the Shiny App, but is needed to prepare network data from expression matrix.
install.packages(cran_packages, repos="http://cran.r-project.org/")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("Mfuzz"))
```
### 2) Clone Quickomics GitHub repository
e.g.
```
git clone https://github.com/YOUR-USERNAME/YOUR-REPOSITORY
```
see more at https://docs.github.com/en/free-pro-team@latest/github/creating-cloning-and-archiving-repositories/cloning-a-repository

### 3) Launch the R Shiny App
Check the follwing web links on various options to launch the app. 
* https://shiny.rstudio.com/articles/running.html
* https://shiny.rstudio.com/deploy/

# Prepare your own data
Two files are needed for each data set, one contains the main data, the other contains network information. For the pre-loaded datasets, the main data files are located in the Quickomics/data/ folder, while the network files are located at Quickomics/networkdata/ folder. 
To prepare the main data files, you need to get the following R data frame objects. 
1. MetaData. Must have sampleid, group, Order, ComparePaires columns. Additional metadata columns about the samples can be added. sampleid values should match with expression data. The group is how you divide samples into biological meaning groups. Order is how you order the group. CompareParies are the comparisons performed. 
2. ProteinGeneName. Must have UniqueID and Gene.Name column. UniqueID values should match with the gene IDs in the data files. Gene.Name should be official gene symbols. Additional columns about the proteins or genes can be added.
3. data_wide. This is the expression matrix; rows are genes, columns are samples. Samples must match sampleid from MetaData, gene IDs must match with UniqueID from ProteinGeneName. 
4. data_long. data_wide converted into long format, with the following four columns: UniqueID, sampleid, expr, group. 

