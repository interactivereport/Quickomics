#Use the QuickOmics.yml file to create a conda environment that has all the required packages for QuickOmics. 
#Load your conda, then run
conda env create -f QuickOmics.yml
#The above command may take a while, then you will have a new conda environment called QuickOmics

#Activate the environment, and you can run the QuickOmics Shiny App
conda activate QuickOmics
#Use your own QuickOmics path and port number below
nohup R -e 'shiny::runApp("/storage/share/ShinyApps/Quickomics/",port=8186,host="0.0.0.0")' >runGitHubcopy.out &

#######################
#Details
#The QuickOmics conda environment also includes a few packages for running RNA-Seq analysis (like DESeq2) using our ExpressionAnalysis package (to be published).
#If you don't have conda installed, you can do so using the example scripts below
#replace /share/miniconda to the directory where you want miniconda to be installed
cd /share/miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p /share/miniconda/Miniconda3
rm Miniconda3-latest-Linux-x86_64.sh

#The codes below shows how the conda environment was created from scratch. 
#Inspired by iamh2o  https://github.com/iamh2o/Quickomics/tree/master/environment
conda install -y -n base -c conda-forge mamba # Install mamba (an accelerated conda wrapper)

# Create conda environment
mamba create -y -n QuickOmics -c conda-forge -c bioconda -c anaconda -c r r-shiny r-shinythemes r-shinyjs plotly  r-plotly r-reshape2 r-tidyverse r-gplots r-ggpubr r-gridextra r-ggrepel  r-rcolorbrewer  r-pheatmap r-rgl r-car r-colourpicker r-venndiagram  r-factoextra r-openxlsx r-visnetwork r-cowplot r-circlize bioconductor-complexheatmap  r-svglite r-shinyjqui r-hmisc r-ggrastr r-ggextra r-networkd3 r-vctrs r-ragg r-textshaping r-biocmanager bioconductor-mfuzz bioconductor-pathview bioconductor-biomart r-psych r-broom r-png
#Add additional packages for ExpressionAnalysis
mamba install -c conda-forge -c bioconda -c anaconda r-r.utils bioconductor-genomicranges bioconductor-edger bioconductor-limma bioconductor-DESeq2 bioconductor-biocparallel bioconductor-apeglm r-caret
#######################

