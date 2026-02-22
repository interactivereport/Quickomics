FROM rocker/r-base:4.3.3

ENV CRAN=https://cloud.r-project.org
ENV R_LIBS_USER=/opt/Rlibs

# Create custom R library
RUN mkdir -p $R_LIBS_USER && chmod -R 777 $R_LIBS_USER

# Install all system dependencies (complete set)

RUN apt-get update && apt-get install -y \
    pandoc \
    build-essential \
    gcc-12 g++-12 gfortran \
    binutils \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libgit2-dev \
    libyaml-dev \
    libx11-dev \
    libxt-dev \
    libpango1.0-dev \
    libreadline-dev \
    libblas-dev \
    liblapack-dev \
    libgfortran5 \
    libgfortran-12-dev \
    libfftw3-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Set gcc/g++ alternatives
RUN update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-12 100 && \
    update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-12 100

# Install CRAN packages
RUN R -e "install.packages(c( \
    'ggplot2', 'patchwork', 'gridExtra', 'gridtext', \
    'data.table', 'dplyr', 'tibble', 'yaml' \
), repos=Sys.getenv('CRAN'))"

# Install BiocManager
RUN R -e "install.packages('BiocManager', repos=Sys.getenv('CRAN'))"

# Install WGCNA dependencies from Bioconductor
RUN R -e "pkgs <- c('impute','preprocessCore'); \
          BiocManager::install(pkgs, ask=FALSE, update=FALSE); \
          if (!all(pkgs %in% rownames(installed.packages()))) quit(status=1)"

# Install WGCNA (Bioconductor version)
RUN R -e "BiocManager::install('WGCNA', ask=FALSE, update=FALSE); \
          if (!'WGCNA' %in% rownames(installed.packages())) quit(status=1)"
          
RUN R -e "install.packages('tinytex'); tinytex::install_tinytex()"
RUN R -e "tinytex::tlmgr_install(c('collection-fontsrecommended', 'xetex'))"

