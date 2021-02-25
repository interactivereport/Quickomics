# Quickomics: Smart query tool to explore omics data in an intuitive, interactive and informative manner

Demo site: http://quickomics.bxgenomics.com

Tutorial: https://bit.ly/3rXIyhL

Interactive Figure 1: https://interactivereport.github.io/Quickomics/Figure1.html
![https://interactivereport.github.io/Quickomics/Figure1.html](https://interactivereport.github.io/Quickomics/Figure1.png?raw=true "Quickomics")

Fig. 1. Selected Quickomics functions applied to a dataset of microglial RNA-seq gene expression from three mouse genotypes over time. A) PCA based on full dataset highlights primary sample separation by mouse age at which the cells were isolated. (B) Volcano plot visualizes differentially expressed genes, most of which show reduced expression in 2mo KO compared to 2mo_WT microglia. For spacing purpose, absolute log2FC (Fold Change) and negative log10 adjusted p-value are capped at 1.5 and 15, respectively. (C) Correlation analysis between two comparisons shows that aging and Cx3cr1-KO have a similar effect on gene expression. (D) Pattern clustering identifies subsets of genes with similar expression over the samples. The clustering is mostly driven by age, with the KO genotype having a similar, but smaller effect. (E) Heatmap of all samples allows the identification of gene clusters with expression regulated by age and/or genotype. Key genes and the pathways they belong to are highlighted on the right. (F) After pathway enrichment analysis, KEGG pathways (Kanehisa and Goto, 2000) of interest can be displayed in a cellular context. The color bars with each stripe representing one comparison show log2 fold changes in various comparisons, allowing project-wide insights for patterns of expression. (G) Correlation network shows potential links between genes of interest.


# Installation
### 1) Install the following R packages:
```R
cran_packages=c("shiny", "shinythemes", "shinyjs", "plotly", "reshape2", "tidyverse", "gplots", "ggpubr", "gridExtra", "ggrepel",
"RcolorBrewer", "pheatmap", "rgl", "car", "colourpicker", "VennDiagram", "factoextra",  "openxlsx", "visNetwork",
"cowplot", "circlize", "ComplexHeatmap", "svglite", "shinyjgui", "Hmisc", "ggrastr", "ggExtra")
#Note: Hmisc is not required to run the Shiny app but is needed to prepare network data from expression matrix.
install.packages(cran_packages, repos="http://cran.r-project.org/")  #choose repos based on your location if needed

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("Mfuzz"))

Requirements: shiny >= v1.4.0.2 
```
### 2) Clone Quickomics GitHub repository
```
git clone https://github.com/interactivereport/Quickomics.git
```
see more at https://docs.github.com/en/free-pro-team@latest/github/creating-cloning-and-archiving-repositories/cloning-a-repository

### 3) Launch the R Shiny App
Check the following web links on various options to launch the app. 
* https://shiny.rstudio.com/articles/running.html
* https://shiny.rstudio.com/deploy/

# Upload your own data files in csv (Comma Separated Values) format
Prepare data files in csv format and upload to Quickomics by using "Upload Files" tool from "Select Dataset" tab. The system will automatically process these files and create the R data files. You need the following files:
1. Sample Meta Data File. Sample Meta Data must have sampleid and group columns, with additional columns optional. The sample names in sampleid column must match those used in expression data file.
2. Expression Data File. Expression data should be a matrix of expression values with genes/proteins as rows, and samples as columns. The unique IDs for genes/proteins are in the first column. We recommend using log of normalized expression values (e.g. log2(TPM+1) for RNAseq).
3. Comparison Data File. Comparison data should have five columns, UniqueID, test, Adj.P.Value, P.Value and logFC. These are typically computed using R packages like limma, edgeR or DESeq2. Make sure to rename the column headers to Adj.P.Value, P.Value and logFC. The comparison names are listed in test column.
4. An optional Gene/Protein Name File.  Note the system will create gene/protein annotation based on the unique IDs from data files, so most users don't need to prepare the Gene/Protein Name File. Or user could prepare it with four required columns: id (sequential numbers), UniqueID (match with the IDs in the expression and comparison data file), Gene.Name (official gene symbols), Protein.ID (UniProt protein IDs, or enter empty values for RNA-Seq data). Additional columns (e.g. gene biotype) are optional.

After the data files are uploaded, the system will provide the user a link to explore the result in Quickomics. 

# Prepare R data files by computational biologist
We recommend uploading csv files, which is easier for most users. Nevertheless, Experienced R programmer can create R data files to be uploaded. 
Two files are needed for each data set, one contains the main data, the other contains network information. For the pre-loaded datasets, the main data files are located in the Quickomics/data/ folder, while the network files are located at Quickomics/networkdata/ folder.  One can check the Rdata files (e.g. Mouse_microglia_RNA-Seq.RData) in the data folder to see some examples. 
To prepare the main data files, you need to get the following R data frame objects. 
1. MetaData. Must have sampleid, group, Order, ComparePaires columns. Additional metadata columns about the samples can be added. sampleid values should match with expression data. The group is how you divide samples into biological meaning groups. Order is how you order the group. ComparePairs are the comparisons performed. 
2. ProteinGeneName. Must have UniqueID, Gene.Name and Protein.ID columns. UniqueID values should match with the gene IDs in the data files. Gene.Name should be official gene symbols. Protein.ID is UniProt protein IDs, or enter empty values for RNA-Seq data. Additional columns about the proteins or genes can be added.
3. data_wide. This is the expression matrix; rows are genes, columns are samples. Samples must match sampleid values from MetaData, gene IDs must match with UniqueID from ProteinGeneName. For RNA-Seq, we recommend using log2(TPM+1) as the expression values in data_wide. 
4. data_long. data_wide converted into long format, with the following four columns: UniqueID, sampleid, expr, group. The group values must match those listed in MetaData.
5. results_long. The comparison results in long format. The columns are: UniqueID, test, Adj.P.Value, P.Value and logFC. UniqueID must match with UniqueID values from ProteinGeneName; test column has the comparison names, must match values ComparePairs column of MetaData; the other values are typically computed from packages like DESeq2 or limma, but the data headers must be changed to Adj.P.Value, P.Value and logFC.
6. data_results. This is a summary table, it starts with UniqueID and Gene.Name columns, then the intensity (max or mean expression value from data_wide for each gene), followed by mean and SD expression values for each group, then followed by comparison data (test name added as prefix).

The network data file is computed from 'data_wide' expression matrix by using Hmisc package. See the example R code mentioned below. 

## Example R code to prepare R data files from RNA-Seq results
We have provided the example input files (TPM and count matrix files, sample grouping file, comparison list file) and the R code to generate the main data and network data files. These files are located at demo_files/Example_RNA_Seq_data/ directory. Please note that you may need to modify RNA_Seq_raw2quickomics.R to fit your input files.
* rsem_TPM.txt The TPM matrix. One can also use RPKM matrix if needed. 
* rsem_expected_count.txt  The gene count matrix. We used RSEM counts in this case, but gene count results from other packages (like kallisto, salmon, featureCount) can be used as well. 
* grpID.txt This file lists the group information for each sample
* comparison.txt This list lists the comparisons to perform (group 1 vs group 2 in each row).
* RNA_Seq_raw2quickomics.R The R code to read the data, run differential gene expression analysis using DESeq2, and create main and network data files.

## Example R code to prepare R data files from Proteomics results
We have provided the example input files (normalized protein expression, comparison data, sample information, protein and gene names) and the R code to generate the main data and network data files. These files are located at demo_files/Example_Proteomics_data/ directory. Please note that you may need to modify Proteomics2Quickomics.R to fit your input files.
* NormalizedExpression.csv Normalized protein expression (log2 transformed). 
* ComparisonData.csv  Comparison results. The statistic values are: logFC,  P.Value and Adj.P.Value. This can be created using R packages like limma.
* Sample.csv Sample information file.
* ProteinID_Symbol.csv This file lists the proteinIDs and associate gene symbols.
* Proteomics2Quickomics.R The R code to read the data and create main and network data files.

