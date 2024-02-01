output$upload.files.ui <- renderUI({
 tagList(tags$div(
  tags$p("Prepare your own data files in Excel, save them as csv files and upload here. The system will automatically process the files and create the R data files. You need sample metadata file, expression data file, and comparison data file. The system can create gene/protein annotation based on the IDs from data files, or you can upload your own Gene/Protein Name file.")),
  tags$a(href="RNA_Seq_Demo.zip", "Download RNA-Seq example csv files (200 genes from mouse microglia dataset)"),
  tags$br(),
  tags$a(href="Proteomics_Demo.zip", "Download Proteomics example csv files (200 proteins from AD PD dataset)"),
  tags$hr(),
  textInput("F_project_name", label="Project Name", value=""),
  radioButtons("Fspecies",label="Select Species", choices=c("human","mouse", "rat"), inline = T, selected="human"),
  tags$hr(),
  tags$p("Sample MetaData must have sampleid and group columns, with additional columns optional. The sample names in sampleid column must match the expression data file."),
  fileInput("F_sample", "Sample MetaData File"),
  tags$hr(),
  tags$p("Expression data should be matrix of expression values with genes/proteins as rows, and samples as columns.  The unique IDs for genes/proteins are in the first column. We recommend using log of normalized expression values (e.g. log2(TPM+1). Upload csv file, can be compressed as .gz or .zip file."),
  fileInput("F_exp", "Expression Data File"),
  tags$hr(),
  tags$p("Comparison data should have five columns, UniqueID, test, Adj.P.Value, P.Value and logFC. The comparison names are listed in test column. Upload csv file, can be compressed as .gz or .zip file."),
  fileInput("F_comp", "Comparison Data File"),
  tags$hr(),
  checkboxInput("F_annot_auto", "Create Gene/Protein Name File automatically (or uncheck to upload your own file)", TRUE, width="90%"),
  #radioButtons("F_annot_auto", label="Create Gene/Protein Names automatically:", inline = TRUE, choices = c("Yes","No"), selected = "Yes"),
  conditionalPanel(condition="input.F_annot_auto==1",
                   radioButtons("F_ID_type",label="Unique ID Type in the Data Files", choices=c("Ensembl Gene ID", "Gene Symbol", "NCBI GeneID","UniProtKB Protein ID", "UniProt Protein Name"), inline = T, selected="Ensembl Gene ID"),
  checkboxInput("F_ID_info", "Show ID type examples", FALSE, width="90%"),
  conditionalPanel(condition="input.F_ID_info==1",
                   tags$p("The system will extract gene/protein names based on the unique IDs from your data using BioMart or UniProt database. The unique IDs from the expression and comparison data can be one of the following formats:"),
                   tags$ol(
                     tags$li("Ensembl Gene ID (e.g. ENSG00000118260, or with version number ENSG00000118260.14)"), 
                     tags$li("Gene Symbol (e.g. CREB1)"), 
                     tags$li("NCBI GeneID (e.g. 1385)"),
                     tags$li("UniProtKB Protein ID (e.g. P16220, C9J4L5, P16220-2)"), 
                     tags$li("UniProt Protein Name (e.g. CREB1_HUMAN, C9J276_HUMAN)")
                   )),
                   checkboxInput("F_fillName", "Fill in uniqueID when Gene.Name not found", TRUE),
                   checkboxInput("F_description", "Add gene/protein description", FALSE)
  ),
  
   # radioButtons("F_description",label="Add gene/protein description?", choices=c("Yes", "No"), inline = T, selected="No"),
  conditionalPanel(condition="input.F_annot_auto==0",
                   tags$p("The Gene/Protein Name csv file must have four columns: id (sequential numbers), UniqueID (match with the IDs in the expression and comparison data file), Gene.Name (official gene symbols), Protein.ID (UniProt protein IDs, or enter empty values for RNA-Seq data). Additional columns (e.g. gene biotype) can be added."),                  
                   fileInput("F_annot", "Gene/Protein Name File")),
  tags$hr(),
  actionButton("uploadData", "Submit Data")
 )
})
up_message1="The URL for the uploaded dataset will be displayed here once the files are processed."
upload_message(up_message1) 
output$upload.message <- renderText({upload_message()})

showAlert("None")

observeEvent(input$uploadData, {  
  upload_message(up_message1) 
  URL_host <-(session$clientData$url_hostname)
  URL_port <- (session$clientData$url_port)
  URL_protocol<- (session$clientData$url_protocol)
  url_pathname<- (session$clientData$url_pathname)
  if (URL_port=="") {
    URL=str_c(URL_protocol,"//", URL_host) 
  } else {
    URL=str_c(URL_protocol,"//", URL_host, ":", URL_port)
  }
  if (url_pathname!="") {
    URL=str_c(URL, str_replace(url_pathname, "/$", "") ) 
  }
  cleanup_empty<-function(df) {
    df.empty=(is.na(df) | df=="")
    selCol=!(colSums( df.empty)==nrow(df))
    selRow=!(rowSums( df.empty)==ncol(df))
    return(df[selRow, selCol])
  }
  cat(URL, "\n")
  #create unique project ID
  Project_name=input$F_project_name
  ProjectID=str_c("PRJ_",  make.names(Project_name) )
  if (length(ProjectID>45) ) {ProjectID=substr(ProjectID, 1, 45)}
  ProjectID=str_c(ProjectID,"_", stri_rand_strings(1,6) )
  #get expression data
  withProgress(message = 'Processing...', value = 0, {
  MetaData=read.csv(input$F_sample$datapath, header=T, check.names=F)
  MetaData=cleanup_empty(MetaData)
  exp_file=input$F_exp$datapath
  if (str_detect(exp_file, "gz$") ) {exp_file=gzfile(exp_file, "rt")}
  if (str_detect(exp_file, "zip$") ) {
    fnames = as.character(unzip(exp_file, list = TRUE)$Name)
    exp_file=unz(exp_file, fnames[1])
  }
  data_wide=read.csv(exp_file, row.name=1, header=T, check.names=F)
  data_wide=cleanup_empty(data_wide)
  
  comp_file=input$F_comp$datapath
  if (str_detect(comp_file, "gz$") ) {comp_file=gzfile(comp_file, "rt")}
  if (str_detect(comp_file, "zip$") ) {
    fnames = as.character(unzip(comp_file, list = TRUE)$Name)
    comp_file=unz(comp_file, fnames[1])
  }
  results_long=read.csv(comp_file, header=T, check.names=F)
  results_long=cleanup_empty(results_long)
  species=input$Fspecies
  IDs=rownames(data_wide);
  IDs2=results_long$UniqueID
  IDall=unique(c(IDs, IDs2)); 
  IDall=IDall[!is.na(IDall)] 
  IDall=IDall[!(IDall=="")] 
  setProgress(0.1, detail = "Data files loaded. Work on Gene/Protein Name..."); Sys.sleep(0.1)
  #get gene or protein name
  #browser() #bebug
  if (input$F_annot_auto==0) {
    ProteinGeneName=read.csv(input$F_annot$datapath)
    ProteinGeneName=cleanup_empty(ProteinGeneName)
    ProteinGeneName<-ProteinGeneName%>%dplyr::filter(UniqueID %in% IDall)
    if (!"Protein.ID" %in% names(ProteinGeneName)) {ProteinGeneName$Protein.ID=NA} #Add Protein.ID column as it is required for certain tools.
    setProgress(0.3, detail = "Loaded Gene Names. Generate RData file..."); Sys.sleep(0.1)
  } else {
    if (str_detect(input$F_ID_type, "UniProt") ) { #protein name match
      ProteinInfo<-readRDS('db/ProteinInfo.rds')%>%mutate(Gene_Name=str_replace(Gene_Name, " .+", "")) #replace space, as UniProt put alias here
     if (input$F_ID_type=="UniProtKB Protein ID" ) {
      ProteinGeneName<-data.frame(id=1:length(IDall), UniqueID=IDall)%>%left_join(ProteinInfo%>%
                transmute(UniqueID=UniProtKB.AC, Gene.Name=Gene_Name, Protein.ID=UniProtKB.AC, Description=Protein_Name)%>%dplyr::filter(!duplicated(UniqueID)) )
     } else {
       ProteinGeneName<-data.frame(id=1:length(IDall), UniqueID=IDall)%>%left_join(ProteinInfo%>%
                transmute(UniqueID=UniProtKB.ID, Gene.Name=Gene_Name, Protein.ID=UniProtKB.AC, Description=Protein_Name)%>%dplyr::filter(!duplicated(UniqueID)) )
     }
    if (input$F_description==0) {ProteinGeneName<-ProteinGeneName%>%dplyr::select(-Description)} 
    if (input$F_fillName==1) {ProteinGeneName<-ProteinGeneName%>%mutate(Gene.Name=ifelse(is.na(Gene.Name), UniqueID, Gene.Name) ) } 
    setProgress(0.2, detail = "Loaded Protein Names"); Sys.sleep(0.1)
      
      
    } else { #gene
      cat("working on ",species," genes for project", ProjectID, "\n")
      if (species=="rat") {
       ensembl <- useEnsembl(biomart = "ensembl", dataset="rnorvegicus_gene_ensembl")
      } else if (species=="mouse") {
        ensembl <- useEnsembl(biomart = "ensembl", dataset="mmusculus_gene_ensembl")
      } else {
        ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
      }
      setProgress(0.2, detail = "Connected to Biomart, converting IDs to gene names..."); Sys.sleep(0.1)
      #browser() #bebug
      if (input$F_ID_type=="Ensembl Gene ID" ) {
        filter_type="ensembl_gene_id"
        IDall_old=IDall
        IDall=str_replace(IDall, "\\.\\d+$", "")
        EID=data.frame(IDall_old, IDall)
      } else if (input$F_ID_type=="NCBI GeneID" ) { filter_type="entrezgene_id"
      } else if (input$F_ID_type=="Gene Symbol" ) { filter_type="external_gene_name"}
      E_attributes<-c( 'ensembl_gene_id', "external_gene_name", "gene_biotype","entrezgene_id")
      if (input$F_description==1) {E_attributes<-c(E_attributes, "description") }
      system.time( output<-getBM(attributes = E_attributes,filters = filter_type, values = IDall, mart = ensembl, useCache=FALSE) )
      if (nrow(output)==0) {
        error_message="No gene annotation extracted from Biomart. Did you select the correct Species and Unique ID Type?"
        setProgress(0.9, detail = error_message); Sys.sleep(3)
        upload_message(error_message)}
      validate(need(nrow(output)>0, message = "No gene annotation extracted from Biomart. Did you select the correct Species and Unique ID Type?"))
      
      output<-output%>%arrange(entrezgene_id, ensembl_gene_id) #Favor IDs with smaller numbers
      if (input$F_description==0) {output$description=NA}
      F_TYPE=sym(filter_type)

      ProteinGeneName<-data.frame(id=1:length(IDall), UniqueID=IDall)%>%left_join(output%>%
                          transmute(UniqueID=!!F_TYPE, Gene.Name=external_gene_name, GeneType=gene_biotype, Description=description)%>%dplyr::filter(!duplicated(UniqueID)) )
      if (input$F_ID_type=="Ensembl Gene ID" ) {
        ProteinGeneName<-EID%>%transmute(UniqueID=IDall_old, Unique1=IDall)%>%left_join(ProteinGeneName%>%mutate(Unique1=UniqueID)%>%dplyr::select(-UniqueID))%>%dplyr::select(-Unique1)
      }
      
      ProteinGeneName$Protein.ID=NA
      ProteinGeneName<-ProteinGeneName%>%dplyr::select(id,UniqueID, Gene.Name, Protein.ID, GeneType, Description)
      if (input$F_description==0) {ProteinGeneName<-ProteinGeneName%>%dplyr::select(-Description)}  
      if (input$F_fillName==1) {ProteinGeneName<-ProteinGeneName%>%mutate(Gene.Name=ifelse(is.na(Gene.Name), UniqueID, Gene.Name) ) } 
     # browser() #bebug
      setProgress(0.3, detail = "Loaded Gene Names. Generate RData file..."); Sys.sleep(0.1)
      
    
    }
  }
  #now process data
  data_long <- melt(as.matrix(data_wide))
  colnames(data_long) <- c("UniqueID","sampleid","expr")
  data_long<-data_long%>%mutate(sampleid=as.character(sampleid))
  data_long <- data_long%>%left_join(MetaData%>%dplyr::select(sampleid, group) )
  
  groups=sort(unique(MetaData$group))
  tests=sort(unique(results_long$test))
  if (!("ComparePairs" %in% names(MetaData))) {
    MetaData$ComparePairs=""; MetaData$ComparePairs[1:length(tests)]=tests
  }
  if (!("Order" %in% names(MetaData))) {
   MetaData$Order=""; MetaData$Order[1:length(groups)]=groups
  }
  #browser() #bebug
  data_results <- ProteinGeneName[,c("id", "UniqueID","Gene.Name","Protein.ID")]%>%
    left_join(data.frame(UniqueID=rownames(data_wide), Intensity=apply(data_wide,1,mean))%>%dplyr::filter(!duplicated(UniqueID)))
  sinfo1<-data.frame(sampleid=names(data_wide))%>%left_join(MetaData%>%dplyr::select(sampleid, group))
  for(grp in unique(sinfo1$group) ){
    subdata<-data.frame(UniqueID=rownames(data_wide), t(apply(data_wide[,sinfo1$group==grp, drop=FALSE],1,function(x)return(setNames(c(mean(x),sd(x)),paste(grp,c("Mean","sd"),sep="_"))))), check.names=FALSE )
    data_results<-data_results%>%left_join(subdata%>%dplyr::filter(!duplicated(UniqueID)))
  }
  for (ctr in tests) {
    subdata<-results_long%>%dplyr::filter(test==ctr)%>%dplyr::select(UniqueID, logFC, P.Value, Adj.P.Value)
    names(subdata)[2:4]=str_c(ctr, "_", names(subdata)[2:4])
    data_results<-data_results%>%left_join(subdata%>%dplyr::filter(!duplicated(UniqueID)))
  }
    
  strOut=str_c("unlisted/", ProjectID, ".RData")
  save(data_long,data_results,data_wide,MetaData,ProteinGeneName,results_long,file=strOut)
  cat("File ", strOut, " Saved\n" )
  setProgress(0.5, detail ="Rdata files created. Now working on network, this may take a while...")
  
  #network
  #if data_wide has many genes, trim down to 10K
  if (nrow(data_wide)>10000 ) {
    dataSD=apply(data_wide, 1, function(x) sd(x,na.rm=T))
    dataM=rowMeans(data_wide)
    diff=dataSD/(dataM+median(dataM))
    data_wide=data_wide[order(diff, decreasing=TRUE)[1:10000], ]	 
    cat("reduce gene size to 10K for project ", ProjectID, "\n")
  }
  
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
    network <- network %>% mutate_if(is.factor, as.character) %>%
      dplyr::filter(!is.na(cor) & abs(cor) > 0.8 & p < 0.005)
  }
  if (nrow(network)>2e6) {
    network <- network %>% mutate_if(is.factor, as.character) %>%
      dplyr::filter(!is.na(cor) & abs(cor) > 0.85 & p < 0.005)
  }
  cat(ProjectID," final network size ", nrow(network), "\n" )
  
  save(network,file=str_c("unlisted/", ProjectID, "_network.RData") )
  setProgress(1.6, detail =str_c("Finished computing network. Final nodes: ", nrow(network)))
  
  #Load the data as an unlisted project
  #browser()
  unlisted_project<-data.frame(Name=Project_name,  ShortName=Project_name, 
                               ProjectID, Species=species)
  write.csv(unlisted_project, str_c("unlisted/", ProjectID, ".csv") )
  unlisted_project=read.csv(str_c("unlisted/", ProjectID, ".csv"))
  ProjectInfo$ProjectID=ProjectID
  ProjectInfo$Name=unlisted_project$Name
  ProjectInfo$Species=unlisted_project$Species
  ProjectInfo$ShortName=unlisted_project$ShortName
  ProjectInfo$file1= paste("unlisted/",  ProjectID, ".RData", sep = "")  #data file
  ProjectInfo$file2= paste("unlisted/", ProjectID, "_network.RData", sep = "") #Correlation results

  })
  cat("Finished processing data files for ", ProjectID, ".\n")
  up_message2=str_c("The direct URL for the uploaded dataset is: ", URL, "/?unlisted=", ProjectID)
  upload_message(up_message2) 
  showAlert(str_c(URL, "/?unlisted=", ProjectID))
  
#  alert(str_c("Data files processed successfully; please go to another tab to start using the system. In the future, you can use the following URL to access this dataset: ",URL, "/?unlisted=", ProjectID, " (this URL is also listed under the submit button of this page. Copy the URL and save it for future use." ) )
})

observe({
  Sys.sleep(0.5)
  if (str_detect(showAlert(), "unlisted") ){
    alert(str_c("Data files processed successfully; please go to another tab to start using the system. A direct URL to access this dataset will be shown below the Submit Data button of this page. Copy the URL and save it for future use." ) )
    showAlert("None") 
  }
})

