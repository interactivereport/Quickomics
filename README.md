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
"cowplot", "circlize", "ComplexHeatmap", "svglite", "shinyjgui", "Hmisc", "ggrastr")
#Note: Hmisc is not required to run the Shiny App, but is needed to prepare network data from expression matrix.
install.packages(cran_packages, repos="http://cran.r-project.org/")  #choose repos based on your location if needed

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
Two files are needed for each data set, one contains the main data, the other contains network information. For the pre-loaded datasets, the main data files are located in the Quickomics/data/ folder, while the network files are located at Quickomics/networkdata/ folder.  One can check the Rdata files (e.g. Mouse_microglia_RNA-Seq.RData) in the data folder to see some examples. 
To prepare the main data files, you need to get the following R data frame objects. 
1. MetaData. Must have sampleid, group, Order, ComparePaires columns. Additional metadata columns about the samples can be added. sampleid values should match with expression data. The group is how you divide samples into biological meaning groups. Order is how you order the group. CompareParies are the comparisons performed. 
2. ProteinGeneName. Must have UniqueID and Gene.Name column. UniqueID values should match with the gene IDs in the data files. Gene.Name should be official gene symbols. Additional columns about the proteins or genes can be added.
3. data_wide. This is the expression matrix; rows are genes, columns are samples. Samples must match sampleid values from MetaData, gene IDs must match with UniqueID from ProteinGeneName. For RNA-Seq, we recommend using log2(TPM+1) as the expression values in data_wide. 
4. data_long. data_wide converted into long format, with the following four columns: UniqueID, sampleid, expr, group. The group values must match thoese listed in MetaData.
5. results_long. The comparison results in long format. The columns are: UniqueID, test, Adj.P.Value, P.Value and logFC. UniqueID must match with UniqueID values from ProteinGeneName; test column has the comparison names, must match values ComparePairs column of MetaData; the other values are typically computed from packages like DESeq2 or limmma, but the spelling must be changed to Adj.P.Value, P.Value and logFC.
6. data_results. This is a summary table, it starts with UniqueID and Gene.Name columns, then the intensity (max or mean expression value from data_wide for each gene), followed by mean and SD expresion values for each group, then followed by comparison data (test name added as prefix).

The network data file can be computed from data_wide using Hmisc package. See the example R code mentioned below. 

## R codes to preare Rdata file from RNA-Seq results
We have provided the example input files (TPM and count matrix files, sample grouping file, comprisons list file) and the R code to generate the main data and network data files. These files are located at demo_files/Example_RNA_Seq_data/
* rsem_TPM.txt The TPM matrix. One can also use RPKM matrix if needed. 
* rsem_expected_count.txt  The gene count matrix. We used RSEM counts in this case, but gene count results from other packages (like kallisto, salmon, featureCount) can be used as well. 
* grpID.txt This file lists the group information for each sample
* comparison.txt This list lists the comparisons to perform (group 1 vs group 2 in each row).
* RNA_Seq_raw2quickomics.R The R code to read the data, run DEG using DESeq2, and create main data file. Then the code will use data_wide to generate network data. 