###version 1.0, Nov 12 2020
## from the raw expression data to quickomics object
##
###########
## required: TPM and raw count expression matrix
## sample sheets two columns, first is the sample names matching the column names in TPM and raw counts matrics;
##                            second column is the group information (please order the group information as the display order)
## and comparison sheets two column, first is the numerator second is the denominator
##
require(tidyverse)
require(reshape2)
require(DESeq2)
rm(list=ls())
strTPM <- "rsem_TPM.txt"
strC <- "rsem_expected_count.txt"
strGrp <- "grpID.txt" #sample and group information
strCOM <- "comparison.txt" #the two groups used for comparison
strOut <- "mouse_microglia_RNA-Seq.RData" #output file
genome <- "mm10" #or use hg38 for human
minTPM <- 3  #cutoff for exprssed genes
minTPMsample <- 3 #number of samples that have the gene expressed
#########################
## process ....
#get genome annotation from BioMart. You can use your own annotation source if needed.
strGenome <- paste0(genome,".rds")
if(!file.exists(strGenome)){
  require(biomaRt)
  if(genome=="mm10"){
    ensembl <- useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl", version=91)
    gInfo <- getBM(c("ensembl_gene_id","mgi_symbol", "gene_biotype"),mart = ensembl,useCache=F)
    gInfo <- gInfo[!duplicated(gInfo$ensembl_gene_id),]
  }else if(genome=="hg38"){
    ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", version=97)
    gInfo <- getBM(c("ensembl_gene_id","hgnc_symbol","gene_biotype"),mart = ensembl)
  }else{
    stop("unknown genome!")
  }
  saveRDS(gInfo,strGenome)
}
gInfo <- readRDS(strGenome)
rownames(gInfo) <- gInfo[,1]

Grp <- read.table(strGrp,header=T,sep="\t",as.is=T,row.names=1)
Grp <- setNames(Grp[,1],rownames(Grp))
TPM <- read.table(strTPM,header=T,sep="\t",as.is=T,row.names=1,check.names=F)[,names(Grp)]
TPM <- TPM[apply(TPM,1,function(x)return(sum(x>minTPM)>=minTPMsample)),] #filter for expressed genes
rawC <- read.table(strC,header=T,sep="\t",as.is=T,row.names=1,check.names=F)[rownames(TPM),names(Grp)]
COM <- read.table(strCOM,header=T,sep="\t",as.is=T)

#Optional step, choose the subset of data and comparison to be uploaded
#In this example below, we only choose the age group (2mo, 1yr and 2yr) samples and compare between the genotypes.
Grp<-Grp[1:93] #only use first 93 samplses
TPM<-TPM[, 1:93]
rawC<-rawC[, 1:93]
COM<-COM[1:9, ] #only choose the between genotype comparisons.


rownames(TPM) <- rownames(rawC) <- sapply(strsplit(rownames(rawC),"\\."),head,1) #remove version in Enesmbl geneID
ProteinGeneName <- cbind(id=0:(nrow(TPM)-1),gInfo[rownames(TPM),])
colnames(ProteinGeneName) <- c("id","UniqueID","Gene.Name","Biotype")


TPM <- log2(1+TPM)
data_wide <- TPM
data_long <- melt(as.matrix(TPM))
colnames(data_long) <- c("UniqueID","sampleid","expr")
data_long <- cbind(data_long,group=Grp[data_long$sampleid])

MetaData <- list(sampleid=names(Grp),group=Grp,Order=unique(Grp),ComparePairs=apply(COM,1,paste,collapse="-"))
MetaData <- as.data.frame(lapply(MetaData,'length<-',max(sapply(MetaData,length))),stringsAsFactors=F)
MetaData[is.na(MetaData)] <- ""
#Optional step, add additional metaData columns like age, genotype, gender etc.
MetaData<-MetaData%>%mutate(Age=str_split_fixed(sampleid, "-", 3)[, 1], Genotype=str_split_fixed(sampleid, "-", 3)[, 2], Gender=str_replace(str_split_fixed(sampleid, "-", 3)[, 3], "\\d+", "") )


data_results <- cbind(ProteinGeneName[,c("UniqueID","Gene.Name","id")],Intensity=apply(TPM,1,max))
#compute mean and SD expression for each group
for(i in unique(Grp)){
  data_results <- cbind(data_results,t(apply(TPM[,Grp==i],1,function(x)return(setNames(c(mean(x),sd(x)),paste(i,c("Mean","sd"),sep="_"))))) )
}

#The code below uses DESeq to compute differentially expressed gene statistics.
#The DEG results are saved to two data frames: results_long and data_results.
#You can use results from limma or other packages instead, just prepare results_long and data_results accordingly.
results_long <- NULL
for(i in 1:nrow(COM)){
  strVS <- paste(COM[i,],collapse="vs")
  message(strVS)
  X <- rawC[,Grp%in%COM[i,]]
  Data <- DESeqDataSetFromMatrix(countData=matrix(as.integer(as.matrix(X)),nrow=nrow(X),dimnames=dimnames(X)),
                                 colData=data.frame(row.names=colnames(X),group=Grp[colnames(X)]),
                                 design=~group)
  dds <- DESeq(Data,betaPrior=TRUE,quiet=T)
  res <- as.data.frame(results(dds,c("group",COM[i,1],COM[i,2])))
  res <- res[data_results[,"UniqueID"],c("log2FoldChange","pvalue","padj")]
  colnames(res) <- paste0(strVS,"_DESeq.",c("logFC","P.value","Adj.P.value"))
  data_results <- cbind(data_results,res)

  colnames(res) <- c("logFC","P.Value","Adj.P.Value")
  res <- cbind(UniqueID=rownames(res),test=factor(strVS,levels=apply(COM,1,paste,collapse="vs")),res[,c("Adj.P.Value","P.Value","logFC")])
  if(is.null(results_long)){
    results_long <- res
  }else{
    results_long <- rbind(results_long,res)
  }
}

save(data_long,data_results,data_wide,MetaData,ProteinGeneName,results_long,file=strOut)




##Now prepare network data, it may take a while and you need a server with large memory (32 to 64GB for typical RNA-Seq data)
library(Hmisc)
system.time(cor_res <- Hmisc::rcorr(as.matrix(t(data_wide))) )
cormat <- cor_res$r
pmat <- cor_res$P
ut <- upper.tri(cormat)
network <- tibble (
from = rownames(cormat)[row(cormat)[ut]],
to = rownames(cormat)[col(cormat)[ut]],
cor  = signif(cormat[ut], 2),
p = signif(pmat[ut], 2),
direction = as.integer(sign(cormat[ut]))
)
dim(network)
network <- network %>% mutate_if(is.factor, as.character) %>%
dplyr::filter(!is.na(cor) & abs(cor) > 0.7 & p < 0.05)
dim(network)
#check network size. If it has>5 million rows, use higher cor and lower p to futher reduce size

#Not run for the example data choose the right cutoff to get fewer than 5 million data points
network_back=network
network <- network_back %>% mutate_if(is.factor, as.character) %>%
dplyr::filter(!is.na(cor) & abs(cor) > 0.85 & p < 0.005)
dim(network)
####
ProjectID <- "mouse_microglia_RNA-Seq"
save(network,file=str_c(ProjectID, "_network.RData") )
