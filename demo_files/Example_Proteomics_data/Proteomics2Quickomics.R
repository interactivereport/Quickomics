###version 1.0, Feb 24 2021
## From Proteomics Dat to Quickomics Rdata object
##
###########
## required:
## sample information
## Protein IDs and associated gene symbols
## Normalized Expression Levels
## Comparison Data (logFC, p-value, adj.p.value from packages like limma)
require(tidyverse)
require(reshape2)
library(Hmisc)

sample_file="Sample.csv"
exp_file="NormalizedExpression.csv"
comp_file="ComparisonData.csv"
ProteinID_file="ProteinID_Symbol.csv"
ProjectID="AD_PD_Brain_Proteomics" #give a name for the ouput file

MetaData=read.csv(sample_file, header=T, check.names=F)
data_wide=read.csv(exp_file, row.name=1, header=T, check.names=F)
results_long=read.csv(comp_file, header=T, check.names=F)
ProteinGeneName=read.csv("ProteinID_Symbol.csv")

#format files by changing column headers to those required by QuickOmics
#see https://github.com/interactivereport/Quickomics the README file for description of file format and headers
MetaData<-MetaData%>%transmute(sampleid, group=DiseaseState, PlatformName, Tissue)
ProteinGeneName<-ProteinGeneName%>%transmute(id=1:nrow(ProteinGeneName), UniqueID=Protein.ID, Gene.Name=Gene.Symbol, Protein.ID=Protein.ID )
results_long<-results_long%>%transmute(UniqueID=Protein.ID, test=ComparisonID, logFC, P.Value, Adj.P.Value)

#now process data
data_long <- reshape2::melt(as.matrix(data_wide))
colnames(data_long) <- c("UniqueID","sampleid","expr")
data_long<-data_long%>%mutate(sampleid=as.character(sampleid))
data_long <- data_long%>%left_join(MetaData%>%dplyr::select(sampleid, group) )

groups=sort(unique(MetaData$group))
tests=sort(unique(results_long$test))
MetaData$ComparePairs=""; MetaData$ComparePairs[1:length(tests)]=tests
MetaData$Order=""; MetaData$Order[1:length(groups)]=groups

data_results <- ProteinGeneName[,c("id", "UniqueID","Gene.Name","Protein.ID")]%>%
  left_join(data.frame(UniqueID=rownames(data_wide), Intensity=apply(data_wide,1,mean))%>%filter(!duplicated(UniqueID)))
sinfo1<-data.frame(sampleid=names(data_wide))%>%left_join(MetaData%>%dplyr::select(sampleid, group))
for(grp in unique(sinfo1$group) ){
  subdata<-data.frame(UniqueID=rownames(data_wide), t(apply(data_wide[,sinfo1$group==grp],1,function(x)return(setNames(c(mean(x),sd(x)),paste(grp,c("Mean","sd"),sep="_"))))), check.names=FALSE )
  data_results<-data_results%>%left_join(subdata)
}
for (ctr in tests) {
  subdata<-results_long%>%filter(test==ctr)%>%dplyr::select(UniqueID, logFC, P.Value, Adj.P.Value)
  names(subdata)[2:4]=str_c(ctr, "_", names(subdata)[2:4])
  data_results<-data_results%>%left_join(subdata)
}

strOut=str_c(ProjectID, ".RData")
save(data_long,data_results,data_wide,MetaData,ProteinGeneName,results_long,file=strOut)
cat("File ", strOut, " Saved\n" )

#network
system.time(cor_res <- Hmisc::rcorr(as.matrix(t(data_wide))) ) #120 seconds
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
cat(ProjectID," network size ", nrow(network), "\n" )
network <- network %>% mutate_if(is.factor, as.character) %>%
  dplyr::filter(!is.na(cor) & abs(cor) > 0.7 & p < 0.05)
if (nrow(network)>2e6) {
  network <- network_back %>% mutate_if(is.factor, as.character) %>%
    dplyr::filter(!is.na(cor) & abs(cor) > 0.8 & p < 0.005)
}
if (nrow(network)>2e6) {
  network <- network_back %>% mutate_if(is.factor, as.character) %>%
    dplyr::filter(!is.na(cor) & abs(cor) > 0.85 & p < 0.005)
}
cat(ProjectID," final network size ", nrow(network), "\n" )

save(network,file=str_c(ProjectID, "_network.RData") )
