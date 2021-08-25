#!/bin/bash
repo="QUICKO2"
################################################################################
# Sets up miniconda environment 'QUICKO' to run
################################################################################

if [[ $1 != "$repo" ]]; then

	echo """

    Hello! This is the  __QUICKOMICS___ installer. To allow this script to setup the proper conda environment, please run it with the following arguments: ./setup.sh QUICKO /path/to/your/conda/install (which for most people is ~/conda).  If you do not have conda installed at the location specified, this script will install it for you.

Once your environment is setup, you can source it with 'conda activate QUICKO'.
Then to start the app, from the main directory, you can run: R --vanilla -e "shiny::runApp(host='0.0.0.0',port=9804)"

    """
	exit 15
fi

if [[ $2 == "" ]]; then
	echo "

    .......Please specify a path to find/install conda as argument #2

    "
	exit 33
fi

export CONDA_DIR=$2
echo "(INSTALLING/LOKING FOR) CONDA IN "$CONDA_DIR
sleep 5

ENV_DIR=$(dirname $(realpath $_))

# Install miniconda if you don't have it
if [[ ! -d $CONDA_DIR ]]; then
	wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
	bash Miniconda3-latest-Linux-x86_64.sh -b -p $CONDA_DIR
	rm Miniconda3-latest-Linux-x86_64.sh
	source $CONDA_DIR/etc/profile.d/conda.sh # source the conda env
	conda init bash                          # init the shell
fi

conda deactivate                         # Clean start to be sure we can begin from a fresh bashrc source
source ~/.bashrc                         # source .bashrc
source $CONDA_DIR/etc/profile.d/conda.sh # re-source the conda profile

conda install -y -n base -c conda-forge mamba # Install mamba (an accelerated conda wrapper)

# Create conda environment
$CONDA_DIR/bin/mamba create -y -n $repo -c conda-forge -c bioconda -c anaconda -c r r-shiny r-shinythemes r-shinyjs plotly  r-plotly r-reshape2 r-tidyverse r-gplots r-ggpubr r-gridextra r-ggrepel  r-rcolorbrewer  r-pheatmap r-rgl r-car r-colourpicker r-venndiagram  r-factoextra r-openxlsx r-visnetwork r-cowplot r-circlize bioconductor-complexheatmap bioconductor-interactivecomplexheatmap r-svglite r-shinyjqui r-hmisc r-ggrastr r-ggextra r-networkd3 r-vctrs r-ragg r-textshaping r-biocmanager bioconductor-mfuzz

echo """

This part could take some time, and may throw a lot of warnings. You will be prompted to update packages after a long wait, select 'a' for all.

"""
sleep 6

R --vanilla -e 'BiocManager::install("pathview")'


echo """

       Attempting to launch Quickomics shiny app on port 9804....

"""

cd $ENV_DIR/../
R --vanilla -e "shiny::runApp(host='0.0.0.0',port=9804)"
