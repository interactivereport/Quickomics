###########################################################################################################
## Proteomics Visualization R Shiny App
##
##This software belongs to Biogen Inc. All right reserved.
##
##@file: help.R
##@Developer : Benbo Gao (benbo.gao@Biogen.com)
##@Date : 5/16/2018
##@version 1.0
###########################################################################################################

tutorial_url <- "https://interactivereport.github.io/Quickomics/tutorial/docs/"
tutorial_link <- function(page) {
	paste0("<li><a href=\"", tutorial_url, page, ".html#", page, "\" target=\"_blank\">Tutorial</a></li>")
}


output$help_input <- renderText({
	htmlstr <- paste0("
	<br>
	<ul>
	<li>Select Dataset: open a saved project, upload an RData file (with an optional network file), or build a project from csv files in the Upload Files tab</li>
	<li>Project Overview: project information, groups, comparisons, and the expression data units used in plot labels</li>
	<li>Sample Table: sample groups and other sample attributes</li>
	<li>Result Table: comparison results, including log2 fold change, p value and adjusted p value (Benjamini-Hochberg)</li>
	<li>Data Table: normalized expression data</li>
	<li>Protein Gene Names: UniqueID, gene symbol and protein accession number</li>
	<li>Upload Files: create your own project from expression, sample and comparison files</li>",
	tutorial_link("dataset-module"), "
	</ul>")
	htmlstr
})


output$help_QC <- renderText({
	htmlstr <- paste0("
	<br>
	<ul>
	<li>Groups and samples are selected in the Groups and Samples tab. The changes apply here after you click Apply Changes to Project Data</li>
	<li>Click Save to output to keep a plot for download in the Output tab</li>
	<li>PCA Plot: select 2 of the first 10 principal components. You can color, shape and size dots by sample attributes. Ellipses show the 95% normal probability region of each group</li>
	<li>Covariates: tests which sample attributes are linked to the top principal components. Click Compute/Refresh</li>
	<li>AlignQC: alignment QC plots, shown when the project has these data</li>
	<li>Eigenvalues: percent of variance explained by the first 10 principal components</li>
	<li>PCA 3D Plot and PCA 3D Interactive: the first 3 principal components in 3D. The interactive plot can be saved from its upper-right menu</li>
	<li>Sample-sample Distance: heatmaps of the Euclidean distance and the correlation between samples, and the correlation table</li>
	<li>Dendrograms: sample clustering tree in tree, horizontal or circular layout. You can change the tree cut number and font size</li>
	<li>Box Plot: expression distribution of each sample (up to 2,000 random genes per sample)</li>
	<li>CV Distribution: histogram of the coefficient of variation (CV) of each group, with the median CV marked</li>",
	tutorial_link("qc-plots-module"), "
	</ul>")
	htmlstr
})


output$help_volcano <- renderText({
	htmlstr <- paste0("
	<br>
	<ul>
	<li>Set the fold change (not log2 FC) and p value or adjusted p value cutoffs in the left panel</li>
	<li>DEG Counts: number of up- and down-regulated genes in each comparison. Click a comparison to open its volcano plot</li>
	<li>Volcano Plot (Static): genes passing both cutoffs are red, and genes passing only the p value cutoff are green. You can label the top DEGs, uploaded genes or a gene set, and highlight a gene list in another color</li>
	<li>Volcano Plot (Interactive): hover to see genes, and lasso genes to list them in the table. The selected genes can be sent to the highlight list of the static plot</li>
	<li>DEGs in Two Comparisons: log2 fold changes of two comparisons plotted against each other</li>
	<li>Data Table: results of the genes passing the cutoffs</li>",
	tutorial_link("deg-module"), "
	</ul>")
	htmlstr
})


output$help_heatmap <- renderText({
	htmlstr <- paste0("
	<br>
	<ul>
	<li>Genes can be all genes (random or most variable), DEGs from a comparison, an uploaded gene list or annotated gene file, or a gene set</li>
	<li>Click Plot/Refresh after changing the settings</li>
	<li>Static Heatmap Layout 1: made with the 'ComplexHeatmap' package. Samples can be annotated by attributes</li>
	<li>Static Heatmap Layout 2: made with the 'gplots' package</li>
	<li>Static heatmaps show up to 5,000 genes (picked at random when there are more). Gene names are shown when the number of genes is at most Max Number of Genes to Label</li>
	<li>Interactive Heatmap: made with 'morpheus'. You can sort, cluster and zoom in the browser</li>
	<li>Major options: <br>
	Apply Clustering: cluster rows, columns, both or none<br>
	Apply Scaling: center and scale by row, by column, or none<br>
	Distance Metric and Linkage Algorithm (under Show More Options): how the hierarchical clustering is done</li>",
	tutorial_link("heatmap-module"), "
	</ul>")
	htmlstr
})


output$help_expression <- renderText({
	htmlstr <- paste0("
	<br>
	<ul>
	<li>Browsing: page through the genes that pass the cutoffs of one comparison</li>
	<li>Searched Expression Data: plot selected genes, an uploaded list or a gene set, as box, bar, violin or line plots</li>
	<li>Data Table: expression values of the searched genes</li>
	<li>Result Table: comparison results of the searched genes</li>
	<li>Rank Abundance Curve: genes ranked by expression level, with the searched genes labeled</li>
	<li>Click Plot/Refresh after changing the settings</li>",
	tutorial_link("expression-plot-module"), "
	</ul>")
	htmlstr
})


output$help_geneset <- renderText({
	htmlstr <- paste0("
	<br>
	<ul>
	<li>Gene sets are from the <a href=\"https://www.gsea-msigdb.org/gsea/msigdb\" target=\"_blank\">Molecular Signatures Database (MSigDB)</a> for human, mouse and rat. You can also add a custom gene set or a GMT file</li>
	<li>Select one or more comparisons</li>
	<li>Gene Set Enrichment Analysis (GSEA): uses all genes ranked by log2 fold change or p value</li>
	<li>Over-Representation Analysis (ORA): uses the DEGs that pass the fold change and p value cutoffs, or your own gene list</li>
	<li>Click a gene set in a result table to see its genes in Gene Expression and Gene Set Heatmap</li>
	<li>KEGG Pathway View and Wikipathways View: gene log2 fold changes shown on the pathway</li>
	<li>Dot Plot: top gene sets across comparisons</li>",
	tutorial_link("gene-set-enrichment-module"), "
	</ul>")
	htmlstr
})


output$help_pattern <- renderText({
	htmlstr <- paste0("
	<br>
	<ul>
	<li>Clustering analysis to visualize the expression in different groups</li>
	<li>Genes are the DEGs of a selected comparison, or your uploaded genes</li>
	<li>Cluster Method: <a href=\"https://en.wikipedia.org/wiki/Fuzzy_clustering\" target=\"_blank\"> soft (fuzzy) clustering</a> or
	 <a href=\"https://en.wikipedia.org/wiki/K-means_clustering\" target=\"_blank\"> k-means</a></li>
	<li>Each gene is summarised as its mean log expression in each selected group, in the order shown under Groups to Plot. Genes missing a group are left out.</li>
	<li>Cluster On decides what the clusters are based on:
	  <ul>
	  <li><b>Trend only (z-score)</b> (default): each gene's group means are centred and scaled to mean 0 and SD 1. Genes going up, down or peaking together cluster together whatever their expression level or size of change. Best for grouping genes by direction of change, e.g. Time only DEGs from Time Course Analysis.</li>
	  <li><b>Trend and size of change (centred)</b>: the gene's own mean is subtracted, so values are log2 fold changes vs the gene's average (e.g. +1 = twice its average). Same-trend genes with large and small changes fall in different clusters. Useful to split one trend cluster into strong and weak responders (paste its genes with upload genes and use 2-3 clusters).</li>
	  <li><b>Mean expression (no transform)</b>: previous behaviour. Clusters are driven mostly by expression level.</li>
	  </ul></li>
	<li>Plots show the values that were clustered. For z-score and centred, all K-means panels share one y scale and 0 is the gene's own average; Soft Clustering shows genes with membership &ge; 0.4, coloured by membership.</li>",
	tutorial_link("pattern-clustering-module"), "
	</ul>")
	htmlstr
})


output$help_network <- renderText({
	htmlstr <- paste0("
	<br>
	<ul>
	<li>The network is built from Pearson correlations between genes. Only pairs with |r| &gt; 0.7 and p &lt; 0.05 are kept, using the 10,000 most variable genes</li>
	<li>Pre-computed networks are used when all samples are selected. After samples are removed, the network is computed again, which can take a few minutes</li>
	<li>Select one or more genes, and set the r and p value cutoffs, to find the correlated genes</li>
	<li>visNetwork and networkD3 show the same network in two styles. Data Table lists the correlated pairs</li>",
	tutorial_link("correlation-network-module"), "
	</ul>")
	htmlstr
})


output$help_venn <- renderText({
	htmlstr <- paste0("
	<br>
	<ul>
	<li>Venn diagram of the DEGs of up to 5 comparisons in this project</li>
	<li>Intersection Output lists the genes in each part of the diagram</li>",
	tutorial_link("venn-diagram-module"), "
	</ul>")
	htmlstr
})


output$help_vennp <- renderText({
	htmlstr <- paste0("
	<br>
	<ul>
	<li>Venn diagram of the DEGs of up to 5 comparisons from different projects</li>
	<li>Projects may use different species or ID types, so only gene names are used. Check Upper Case Gene Symbols to match e.g. mouse and human genes</li>",
	tutorial_link("venn-diagram-module"), "
	</ul>")
	htmlstr
})


output$help_WGCNA <- renderText({
	htmlstr <- paste0("
	<br>
	<ul>
	<li>WGCNA finds modules of genes whose expression goes up and down together across samples</li>
	<li>Pre-computed results are shown when the project has them. Change Top N Genes or the Dendrogram Cut Height for Merging and click Run to compute new results (3-10 minutes). Run uses the samples selected in Groups and Samples</li>
	<li>Suggested order:
	  <ul>
	  <li><b>WGCNA QC</b>: check the soft-thresholding power (chosen as the lowest power with scale-free fit R&sup2; &gt; 0.8) and the sample tree. Remove outlier samples in Groups and Samples and click Run again</li>
	  <li><b>WGCNA Result</b>: the gene tree with module colors, and the genes of each module. Grey genes are not in any module. A lower cut height keeps more, smaller modules</li>
	  <li><b>Module Eigengenes</b>: the eigengene (summary profile) of each module in each sample, and how modules relate to each other</li>
	  <li><b>Module-Trait Relationships</b>: correlate module eigengenes with sample attributes. Click a heatmap cell to find the hub genes of that module and trait</li>
	  </ul></li>",
	tutorial_link("wgcna-module"), "
	</ul>")
	htmlstr
})


help_correlation_html <- paste0("
	<br>
	<ul>
	<li>Uses the samples selected in Groups and Samples</li>
	<li>Gene-Gene, Sample-Sample, Group-Group: correlate each pair of the selected genes, samples or groups. Group-Group uses group means</li>
	<li>Genes can be selected by name, from the DEGs of a comparison, an uploaded list or a gene set</li>
	<li>Pearson for linear relationships, Spearman for rank (monotonic) relationships</li>
	<li>Result Table: pairs ranked by correlation. Click View Correlation Plot to see a pair</li>
	<li>Gene-Target: correlate every gene with one target, a gene or a sample attribute. Use Ordered levels for attributes like Low, Mid, High (Spearman is then used)</li>
	<li>Compute separately for each level of: one set of correlations per level of another attribute, shown in the Correlation Heatmap</li>",
	tutorial_link("correlation-analysis-module"), "
	</ul>")

output$help_correlation <- renderText({ help_correlation_html })
output$help_correlation_target <- renderText({ help_correlation_html })


output$help_timecourse <- renderText({
	htmlstr <- paste0("
	<br>
	<ul>
	<li>Uses the samples selected in Groups and Samples. Any attribute with ordered levels can be the time variable</li>
	<li>Time-Series DE Analysis: DESeq2 LRT for RNA-seq or limma F-test for log values such as proteomics.
	Time only finds genes that change over time; Time x Condition finds genes whose time trend differs between conditions.
	DESeq2 uses expression values converted back to TPM-like counts, so its p values are approximate</li>
	<li>DEG Clustering Analysis: groups the significant genes by time pattern (DEGreport degPatterns), shown as cluster plots, tables and heatmaps</li>
	<li>Similar Pattern Search: find genes with the same (or opposite) time pattern as a seed gene. Run the DE analysis first</li>",
	tutorial_link("time-course-analysis-module"), "
	</ul>")
	htmlstr
})
