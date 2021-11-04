observe({
  DataIn = DataQCReactive()
  MetaData=DataIn$MetaData
  graph_width=max(800, nrow(MetaData)*15)
  updateSliderInput(session, "alignQC_width", value = graph_width)
  #select Numeric variables
  meta_num<-select_if(MetaData, is.numeric)
  num_col=colnames(meta_num)
  if (length(num_col)>0) {
    if ("RIN" %in% num_col) {sel_col="RIN"} else {sel_col=num_col[1]}
    updateSelectizeInput(session,'alignQC_var',choices=num_col, selected=sel_col)
  }
})

## intergenic, intronic and exonic -----
alignQC_RA_p <- reactive({
  DataIn = DataQCReactive()
  MetaData=DataIn$MetaData
  selN <- c("Exonic_Rate","Intronic_Rate","Intergenic_Rate")
  validate(need(sum(selN%in%colnames(MetaData))==length(selN), message = "Sample MetaData must have Exonic_Rate, Intronic_Rate and Intergenic_Rate columns in order to make Read Allocation graph." ))
  p=NULL
  if(sum(selN%in%colnames(MetaData))==length(selN)){
    D = melt(as.matrix(MetaData[,colnames(MetaData)%in%selN]))
    p<-ggplot(D,aes(x=Var1,y=value,fill=Var2))+
      geom_bar(position="stack",stat="identity")+
      ylab("Fraction of Reads")+xlab("")+
#      ggtitle("Mapped reads allocation")+
      ylim(0,1)+theme_minimal_grid(font_size=input$alignQC_fontsize)+
      scale_fill_manual(values=c("#66c2a5","#fc8d62","#8da0cb"))+
      theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
            legend.position = "top")+
      guides(fill=guide_legend(title=""))
    # MetaData <- MetaData[,!colnames(MetaData)%in%selN,drop=F]
  }
  p
})

output$alignQC_RA_pp<-renderPlot({alignQC_RA_p()})

output$alignQC_RA_plot <- renderUI({
  plotOutput("alignQC_RA_pp",height = input$alignQC_height, width=input$alignQC_width)
})


## top genes ratio----
observe({
  expU<-exp_unit()
  small_value=as.numeric(str_replace(str_split_fixed(expU, "\\+", 2)[2], "\\)", ""))
  if (!is.na(small_value)) {
    updateNumericInput(session, "convert_small", value=small_value)}
  N_log=as.numeric(str_replace(str_extract(expU, "log\\d"), "log", ""))
  if (!is.na(N_log)) {
    updateNumericInput(session, "convert_logbase", value=N_log)}
})
alignQC_TGR_p <- reactive({
  DataIn = DataQCReactive()
  data_wide=DataIn$tmp_data_wide
  #MetaData=DataIn$MetaData
  if (input$convert_exp=="None") {
    estT=data_wide
  } else {
    small_value=input$convert_small
    N_log=input$convert_logbase
    #default data_wide=log2(TPM+1)
    estT=N_log^data_wide - small_value
  }
  topN=c(1,10,30, 100)
  topN <- setNames(topN,paste0("Top",topN))
   D <- t(apply(estT,2,function(x){
    x <- sort(x,decreasing=T)
    return(sapply(topN,function(i)return(sum(x[1:i])/sum(x)*100)))
    }))
  D <- cbind(sID=rownames(D),data.frame(D))
  D2=melt(D, id.vars="sID")
  D2$sID=factor(D2$sID, levels=colnames(data_wide))
  p<-ggplot(D2, aes(x=sID, y=value, color=variable, group=variable))+geom_line(size=2)
  p<-p+ labs(x="", y=input$alignQC_TGR_Y, color="Top Genes")+
    theme_minimal_grid(font_size=input$alignQC_fontsize)+
    theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))
  p
})

output$alignQC_TGR_pp<-renderPlot({alignQC_TGR_p()})
output$alignQC_TGR_plot <- renderUI({
  plotOutput("alignQC_TGR_pp", height = input$alignQC_height, width=input$alignQC_width)
})



##Top gene list
inputReactive<-reactive({
  alignQC_TG_Ngene<-input$alignQC_TG_Ngene
  alignQC_TG_Ntotal<-input$alignQC_TG_Ntotal
  alignQC_fontsize<-input$alignQC_fontsize
  return(list(alignQC_TG_Ngene=alignQC_TG_Ngene, 
              alignQC_TG_Ntotal=alignQC_TG_Ntotal, alignQC_fontsize= alignQC_fontsize ))
})

alignQC_TG_data <- reactive({
  DataIn = DataQCReactive()
  data_wide=DataIn$tmp_data_wide
  inputData<-inputReactive%>%debounce(1000)
  if (input$convert_exp=="None") {
    estT=data_wide
  } else {
    small_value=input$convert_small
    N_log=input$convert_logbase
    #default data_wide=log2(TPM+1)
    estT=N_log^data_wide - small_value
  }
 # browser() #debug
  topUnion <- inputData()$alignQC_TG_Ngene
  topG <- unique(as.vector(apply(estT,2,function(x)return(names(sort(x,decreasing=T)[1:topUnion])))))
  while(length(topG)>inputData()$alignQC_TG_Ntotal){
    topUnion <- topUnion-1
    topG <- unique(as.vector(apply(estT,2,function(x)return(names(sort(x,decreasing=T)[1:topUnion])))))
  }
  estTsum <- apply(estT,2,sum)
  D <- apply(estT[topG,],1,function(x)return(100*x/estTsum))
  graph_height=max(700, ncol(D)*15)
  updateNumericInput(session, "alignQC_TG_height", value=graph_height)
  return(list(D=D, topUnion=topUnion, alignQC_fontsize=inputData()$alignQC_fontsize))
})

alignQC_TG_p <- reactive({
  D=alignQC_TG_data()$D
  topUnion=alignQC_TG_data()$topUnion
  ProteinGeneName = DataReactive()$ProteinGeneName
  D2 <- melt(D[,order(apply(D,2,median))])
  D2<-D2%>%mutate(Var2=as.character(Var2))%>%left_join(ProteinGeneName%>%transmute(Var2=UniqueID, Gene.Name))%>%mutate(ID=str_c(Var2, " ", Gene.Name))
 # browser() #debug
  D2$ID=factor(D2$ID, levels=unique(D2$ID))
  #write.csv(D,file=gsub("pdf","unionTop.csv",strPDF),row.names=F)
  p <- ggplot(D2,aes(x=value,y=ID))+
    geom_point(color="grey50",alpha=0.4,size=1)+
    geom_boxplot(color="#ff7f00",outlier.shape = NA,alpha=0)+
    xlab(input$alignQC_TGR_Y)+ylab("")+
    ggtitle(paste("Top", ncol(D), "Expressed Genes (Union of Top",topUnion,"Genes Per Sample)"))+
    theme_minimal_grid(font_size=alignQC_TG_data()$alignQC_fontsize)
 return(p)
})

output$alignQC_TG_pp<-renderPlot({alignQC_TG_p()})
output$alignQC_TG_plot <- renderUI({
  if (input$alignQC_TG_list==0) {
    plotOutput("alignQC_TG_pp", height = input$alignQC_TG_height, width=input$alignQC_TG_width)
  } else {
    output$TG_table <- DT::renderDT(server=FALSE,{
      D<-alignQC_TG_data()$D
      ProteinGeneName = DataReactive()$ProteinGeneName
      D<-round(t(D)*100)/100
      annot<-data.frame(UniqueID=rownames(D))%>%left_join(ProteinGeneName%>%dplyr::select(UniqueID, Gene.Name))%>%
          mutate(Median_Percentage=apply(D, 1, median, na.rm=T))
      results=cbind(annot, D)
      DT::datatable(results,  extensions = 'Buttons',
                    options = list(
                      dom = 'lBfrtip', pageLength = 50,
                      buttons = list(
                        list(extend = "csv", text = "Download Page", filename = "Page_results",
                             exportOptions = list(modifier = list(page = "current"))),
                        list(extend = "csv", text = "Download All", filename = "All_Results",
                             exportOptions = list(modifier = list(page = "all")))
                      )
                    ),rownames= F)
    })
    dataTableOutput('TG_table')
  }
})

#Plot Other Variables
alignQC_OV_p <- reactive({
  req(input$alignQC_var)
  DataIn = DataQCReactive()
  MetaData=DataIn$MetaData
  MetaData$sampleid=factor(MetaData$sampleid, levels=MetaData$sampleid)
  i=input$alignQC_var
  p<-ggplot(MetaData,aes_string(x="sampleid",y=i))+
    geom_bar(stat="identity")+
    xlab("")+ylab("")+
    ggtitle(i)+
    theme_minimal_grid(font_size=input$alignQC_fontsize)+
    theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
          legend.position = "none")
  p
})
output$alignQC_OV_pp<-renderPlot({alignQC_OV_p()})
output$alignQC_OV_plot <- renderUI({
  plotOutput("alignQC_OV_pp",height = input$alignQC_height, width=input$alignQC_width)
})


