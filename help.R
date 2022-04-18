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


output$help_input <- renderText({
	htmlstr <- "
	<br>
	<ul>
	<li>Result Table: Statistics result, include log2 Fold Change, p value, p-values adjusted (Benjamini-Hochberg)</li>
	<li>Data Table: Normalized data</li>
	<li>Sample Table: sample group and comparison information</li>
	<li>Protein Gene Names: protein accession number and gene symbol</li>
    <li><a href=\"https://interactivereport.github.io/Quickomics/tutorial/docs/dataset-module.html#result-table\">Tutorial</a></li>
	</ul>"
	htmlstr
})


output$help_QC <- renderText({
	htmlstr <- "
	<br>
	<ul>
	<li>You can delete and add groups or samples. For some plots, you can re-arrange the order by delete, and add it back to particular position. This change will affect most plots</li>
	<li>If there is save to output option, you can click to save in the session and output as pdf file </li>
	<li>Dendrograms: 3 layouts: Tree, horizontal,  circular. You can change tree cut number and font size</li>
	<li>PCA Plot:  you can select principal components (up to 5) to visualize. The size of the concentration ellipse in normal probability is 95%</li>
	<li>PCA 3D Plot: This plot can be captured by print screen from keyboard</li>
	<li>PCA 3D Interactive: you can save this plot from upper-right menu</li>
	<li>Sample-sample Distance: the Euclidean distance computed from expression values of the samples</li>
	<li>Box Plot: The box plot (a.k.a. box and whisker diagram) is a standardized way of displaying the distribution of data based on the five number summary: minimum, first quartile, median, third quartile, and maximum</li>
	<li>CV Distribution:  show the histogram of CV and median of CV for each group</li>
    <li><a href=\"https://interactivereport.github.io/Quickomics/tutorial/docs/qc-plots-module.html#pca-plot\">Tutorial</a></li>
	</ul>"
	htmlstr
})



output$help_volcano <- renderText({
	htmlstr <- "
	<br>
	<ul>
	<li>Volcano Plot (Static): You can select the comparison group, fold change (not log2 FC) cutoff, p value cutoff. The selected proteins will be showed in red. Top 50 genes will be labeled </li>
	<li>Volcano Plot (Interactive)</li>
	<li>Data Output: selected genes data</li>
    <li><a href=\"https://interactivereport.github.io/Quickomics/tutorial/docs/deg-module.html#volcano-plot-static\">Tutorial</a></li>
	</ul>"
	htmlstr
})


output$help_heatmap <- renderText({
	htmlstr <- "
	<br>
	<ul>
	<li>You can select all data,subset, upload genes to visualize.The uploaded gene should be gene names, accession number</li>
	<li>If the selected gene or upload gene number is less than 50, it will show gene-accession names</li>
	<li>Static Heatmap 1: created by using package 'ComplexHeatmap' </li>
	<li>Static Heatmap 2: created by using package 'gplots'</li>
	<li>Interactive Heatmap: created by using package 'heatmaply'<br>
	This function is slow if you select large number of genes.</li>
	<li>Major options: <br>
	Apply Clustering: how the row and column dendrogram should be reordered<br>
	Distance Metric: : function used to compute the distance (dissimilarity) between both rows and columns. <br>
	Linkage Algorithm: methods of agglomerative hierarchical clustering<br>
	Scale: If the values should be centered and scaled in either the row direction or the column direction, or none. </li>
    <li><a href=\"https://interactivereport.github.io/Quickomics/tutorial/docs/heatmap-module.html#heatmap-module\">Tutorial</a></li>
	</ul>"
	htmlstr
})


output$help_expression <- renderText({
	htmlstr <- "
	<br>
	<ul>
	<li>Browsing: view all proteins pass the cutoff. </li>
	<li>Searched Protein Expression: You can search multiple proteins to plot </li>
	<li>Data Output: show the selected protein expression </li>
	<li>Rank Abundance Curve: helps interpret the distribution of relative abundance and expression levels of a set of genes.</li>
	<li>Options are straightforward </li>
    <li><a href=\"https://interactivereport.github.io/Quickomics/tutorial/docs/expression-plot-module.html#expression-plot-module\">Tutorial</a></li>
	</ul>"
	htmlstr
})



output$help_geneset <- renderText({
	htmlstr <- "
	<br>
	<ul>
	<li>The datasets used are from <a href=\"http://software.broadinstitute.org/gsea/msigdb/collections.jsp\" target=\"_blank\"> Molecular Signatures Database (MSigDB)</a></li>
	<li>Select genes for gene set enrichment analysis: by comparison, fold change, p value</li>
	<li>Select interested dataset</li>
	<li>Click one of the gene set, it will show the expression of genes in this dataset and the heatmap from all the studied groups </li>
	<li>If the gene sets are KEGG pathway, it will label the gene as well relative expression on the KEGG pathway. Save pathway by right click </li>
    <li><a href=\"https://interactivereport.github.io/Quickomics/tutorial/docs/gene-set-enrichment-module.html#gene-set-enrichment-module\">Tutorial</a></li>
	</ul>"
	htmlstr
})


output$help_pattern <- renderText({
	htmlstr <- "
	<br>
	<ul>
	<li>Clustering analysis to visualize the expresion in different groups</li>
	<li>The selection of genes used clustering analysis is based on any comparison in this experiment</li>
	<li>Cluster Mehtod: <a href=\"https://en.wikipedia.org/wiki/Fuzzy_clustering\" target=\"_blank\"> soft (fuzzy) clustering </a>,
	 <a href=\"https://en.wikipedia.org/wiki/K-means_clustering\" target=\"_blank\"> k-means</a>, 
	 <a href=\"https://en.wikipedia.org/wiki/K-medoids\" target=\"_blank\">  partitioning around medoids (PAM)</a></li>
	<li>program determined the optimal number of clusters</li>
	<li>If the gene sets are KEGG pathway, it will label the gene as well relative expression on the KEGG pathway. Save pathway by right click </li>
    <li><a href=\"https://interactivereport.github.io/Quickomics/tutorial/docs/pattern-clustering-module.html#pattern-clustering-module\">Tutorial</a></li>
	</ul>"
	htmlstr
})

output$help_network <- renderText({
	htmlstr <- "
	<br>
	<ul>
	<li>The co-expression network construction is based on protein-protein correlation matrix information. In order to speed up the response time, only keep correlation pairs with |r| > 0.6  </li>
	<li>You can search genes to find out expression correlated genes</li>
	<li>Two network visualization methods provided </li>
    <li><a href=\"https://interactivereport.github.io/Quickomics/tutorial/docs/pattern-clustering-module.html#pattern-clustering-module\">Tutorial</a></li>
	</ul>"
	htmlstr
})


output$help_venn <- renderText({
	htmlstr <- "
	<br>
	<ul>
	<li>Used to generate venn diagram among different comparsions in the experiment</li>
	<li>data output are intersection among different groups </li>
    <li><a href=\"https://interactivereport.github.io/Quickomics/tutorial/docs/venn-diagram-module.html#venn-diagram-module\">Tutorial</a></li>
	</ul>"
	htmlstr
})


output$help_vennp <- renderText({
	htmlstr <- "
	<br>
	<ul>
	<li>Used to generate venn diagram among different experiments in the list</li>
	<li>For different expriments, the samples maybe from different species, and used different database to search. Only gene names are used for intersection</li>
    <li><a href=\"https://interactivereport.github.io/Quickomics/tutorial/docs/venn-diagram-module.html#venn-diagram-module\">Tutorial</a></li>
	</ul>"
	htmlstr
})

#
