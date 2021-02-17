###########################################################################################################
## Proteomics Visualization R Shiny App
##
##This software belongs to Biogen Inc. All right reserved.
##update 09/04/2020 to add DEGs of TWo Comparisons, control on range of logFC and -log10(P-value)
##@file: volcano.R
##@Developer : Benbo Gao (benbo.gao@Biogen.com)
##@Date : 5/16/2018
##@version 1.0
###########################################################################################################

observe({
  DataIn = DataReactive()
  tests = DataIn$tests
  ProteinGeneName = DataIn$ProteinGeneName
  updateRadioButtons(session,'volcano_genelabel', inline = TRUE, choices=colnames(ProteinGeneName)[-1], selected="Gene.Name")
  updateSelectizeInput(session,'volcano_test',choices=tests, selected=tests[1])
  updateSelectizeInput(session,'volcano_test1',choices=tests, selected=tests[1])
  if (length(tests)>1) {	updateSelectizeInput(session,'volcano_test2',choices=tests, selected=tests[2])}
  else {	updateSelectizeInput(session,'volcano_test2',choices=tests, selected=tests[1])}
})


observe({
  DataIn = DataReactive()
  results_long = DataIn$results_long
  test_sel = input$volcano_test
  FCcut = log2(as.numeric(input$volcano_FCcut))
  pvalcut = as.numeric(input$volcano_pvalcut)
  if (input$volcano_psel == "Padj") {
    tmpdat = results_long %>% filter(test==test_sel & Adj.P.Value < pvalcut & abs(logFC) > FCcut) 
  } else {
    tmpdat = results_long %>% filter(test==test_sel & P.Value < pvalcut & abs(logFC) > FCcut) 
  }
  output$volcano_filteredgene <- renderText({paste("Genes pass cutoff (DEGs):",nrow(tmpdat),sep="")})
  #browser()#debug
  DEGs=tmpdat$UniqueID
  if (nrow(tmpdat)>input$Ngenes) {
    DEGs=sample(DEGs, input$Ngenes)
  }
  updateTextAreaInput(session, "volcano_gene_list", value=paste(DEGs, collapse="\n"))
})


DatavolcanoReactive <- reactive({
  DataIn = DataReactive()
  results_long = DataIn$results_long
  
  test_sel = input$volcano_test
  FCcut = log2(as.numeric(input$volcano_FCcut))
  FCcut_rd=round(FCcut*1000)/1000
  pvalcut = as.numeric(input$volcano_pvalcut)
  volcano_genelabel = input$volcano_genelabel
  
  res = results_long %>% filter(test==test_sel) %>%
    filter(!is.na(P.Value)) %>%
    dplyr::mutate (color="Not Significant") %>% as.data.frame() 
  
  
  res$labelgeneid = res[,match(volcano_genelabel,colnames(res))]
  
  if (input$volcano_psel == "Padj") {
    res$color[which((abs(res$logFC)>FCcut)*(res$Adj.P.Value<pvalcut)==1)] = paste0("Padj","<",pvalcut," & abs(log2FC)>",FCcut_rd)
    res$color[which((abs(res$logFC)<FCcut)*(res$Adj.P.Value<pvalcut)==1)] =  paste0("Padj","<",pvalcut, " & abs(log2FC)<",FCcut_rd)
    res$color = factor(res$color,levels = unique(c("Not Significant",	paste0("Padj","<",pvalcut, " & abs(log2FC)<",FCcut_rd),	paste0("Padj","<",pvalcut, " & abs(log2FC)>",FCcut_rd))))
    if (input$Max_Pvalue>0) {
      res<-res%>%mutate(Adj.P.Value=pmax(Adj.P.Value, 10^(0-input$Max_Pvalue) ))
    }
  } else { 
    res$color[which((abs(res$logFC)>FCcut)*(res$P.Value<pvalcut)==1)] = paste0("pval","<",pvalcut," & abs(log2FC)>",FCcut_rd)
    res$color[which((abs(res$logFC)<FCcut)*(res$P.Value<pvalcut)==1)] =  paste0("pval","<",pvalcut, " & abs(log2FC)<",FCcut_rd)
    res$color = factor(res$color,levels = unique(c("Not Significant",	paste0("pval","<",pvalcut, " & abs(log2FC)<",FCcut_rd),	paste0("pval","<",pvalcut, " & abs(log2FC)>",FCcut_rd))))
    if (input$Max_Pvalue>0) {
      res<-res%>%mutate(P.Value=pmax(P.Value, 10^(0-input$Max_Pvalue) ))
    }
  }
  
  res$logFC_ori=res$logFC
  if (input$Max_logFC>0) {
    res<-res%>%mutate(logFC=ifelse(logFC>=0, pmin(input$Max_logFC, logFC), pmax(0-input$Max_logFC, logFC) ) )
  }
  
  return(res)
})

DatavolcanoReactive1 <- reactive({
  DataIn = DataReactive()
  results_long = DataIn$results_long
  
  test_sel = input$volcano_test1
  FCcut = log2(as.numeric(input$volcano_FCcut))
  FCcut_rd=round(FCcut*1000)/1000
  pvalcut = as.numeric(input$volcano_pvalcut)
  volcano_genelabel = input$volcano_genelabel
  
  res = results_long %>% filter(test==test_sel) %>%
    filter(!is.na(P.Value)) %>%
    dplyr::mutate (color="Not Significant") %>% as.data.frame() 
  
  res$labelgeneid = res[,match(volcano_genelabel,colnames(res))]
  res$Sig="X_notsig"
  
  if (input$volcano_psel == "Padj") {
    res$Sig[which((abs(res$logFC)>FCcut)*(res$Adj.P.Value<pvalcut)==1)] = "X_sig"
    res$color[which((abs(res$logFC)>FCcut)*(res$Adj.P.Value<pvalcut)==1)] = paste0("Padj","<",pvalcut," & abs(log2FC)>",FCcut_rd)
    res$color[which((abs(res$logFC)<FCcut)*(res$Adj.P.Value<pvalcut)==1)] =  paste0("Padj","<",pvalcut, " & abs(log2FC)<",FCcut_rd)
    res$color = factor(res$color,levels = unique(c("Not Significant",	paste0("Padj","<",pvalcut, " & abs(log2FC)<",FCcut),	paste0("Padj","<",pvalcut, " & abs(log2FC)>",FCcut_rd))))
    if (input$Max_Pvalue>0) {
      res<-res%>%mutate(Adj.P.Value=pmax(Adj.P.Value, 10^(0-input$Max_Pvalue) ))
    }
  } else { 
    res$Sig[which((abs(res$logFC)>FCcut)*(res$P.Value<pvalcut)==1)] = "X_sig"
    res$color[which((abs(res$logFC)>FCcut)*(res$P.Value<pvalcut)==1)] = paste0("pval","<",pvalcut," & abs(log2FC)>",FCcut_rd)
    res$color[which((abs(res$logFC)<FCcut)*(res$P.Value<pvalcut)==1)] =  paste0("pval","<",pvalcut, " & abs(log2FC)<",FCcut_rd)
    res$color = factor(res$color,levels = unique(c("Not Significant",	paste0("pval","<",pvalcut, " & abs(log2FC)<",FCcut_rd),	paste0("pval","<",pvalcut, " & abs(log2FC)>",FCcut_rd))))
    if (input$Max_Pvalue>0) {
      res<-res%>%mutate(P.Value=pmax(P.Value, 10^(0-input$Max_Pvalue) ))
    }
  }
  if (input$Max_logFC>0) {
    res<-res%>%mutate(logFC=ifelse(logFC>=0, pmin(input$Max_logFC, logFC), pmax(0-input$Max_logFC, logFC) ) )
  }
  return(res)
})

DatavolcanoReactive2 <- reactive({
  DataIn = DataReactive()
  results_long = DataIn$results_long
  
  test_sel = input$volcano_test2
  FCcut = log2(as.numeric(input$volcano_FCcut))
  FCcut_rd=round(FCcut * 1000)/1000
  pvalcut = as.numeric(input$volcano_pvalcut)
  volcano_genelabel = input$volcano_genelabel
  
  res = results_long %>% filter(test==test_sel) %>%
    filter(!is.na(P.Value)) %>%
    dplyr::mutate (color="Not Significant") %>% as.data.frame() 
  if (input$Max_logFC>0) {
    res<-res%>%mutate(logFC=ifelse(logFC>=0, pmin(input$Max_logFC, logFC), pmax(0-input$Max_logFC, logFC) ) )
  }
  res$labelgeneid = res[,match(volcano_genelabel,colnames(res))]
  res$Sig="Y_notsig"
  
  if (input$volcano_psel == "Padj") {
    res$Sig[which((abs(res$logFC)>FCcut)*(res$Adj.P.Value<pvalcut)==1)] = "Y_sig"
    res$color[which((abs(res$logFC)>FCcut)*(res$Adj.P.Value<pvalcut)==1)] = paste0("Padj","<",pvalcut," & abs(log2FC)>",FCcut_rd)
    res$color[which((abs(res$logFC)<FCcut)*(res$Adj.P.Value<pvalcut)==1)] =  paste0("Padj","<",pvalcut, " & abs(log2FC)<",FCcut_rd)
    res$color = factor(res$color,levels = unique(c("Not Significant",	paste0("Padj","<",pvalcut, " & abs(log2FC)<",FCcut_rd),	paste0("Padj","<",pvalcut, " & abs(log2FC)>",FCcut_rd))))
    if (input$Max_Pvalue>0) {
      res<-res%>%mutate(Adj.P.Value=pmax(Adj.P.Value, 10^(0-input$Max_Pvalue) ))
    }
  } else { 
    res$Sig[which((abs(res$logFC)>FCcut)*(res$P.Value<pvalcut)==1)] = "Y_sig"
    res$color[which((abs(res$logFC)>FCcut)*(res$P.Value<pvalcut)==1)] = paste0("pval","<",pvalcut," & abs(log2FC)>",FCcut_rd)
    res$color[which((abs(res$logFC)<FCcut)*(res$P.Value<pvalcut)==1)] =  paste0("pval","<",pvalcut, " & abs(log2FC)<",FCcut_rd)
    res$color = factor(res$color,levels = unique(c("Not Significant",	paste0("pval","<",pvalcut, " & abs(log2FC)<",FCcut_rd),	paste0("pval","<",pvalcut, " & abs(log2FC)>",FCcut_rd))))
    if (input$Max_Pvalue>0) {
      res<-res%>%mutate(P.Value=pmax(P.Value, 10^(0-input$Max_Pvalue) ))
    }
  }
  
  if (input$Max_logFC>0) {
    res<-res%>%mutate(logFC=ifelse(logFC>=0, pmin(input$Max_logFC, logFC), pmax(0-input$Max_logFC, logFC) ) )
  }
  
  return(res)
})



output$volcanoplot <- renderPlotly({
  res = DatavolcanoReactive()
  test_sel = input$volcano_test
  FCcut = log2(as.numeric(input$volcano_FCcut))
  pvalcut = as.numeric(input$volcano_pvalcut)
  if (input$volcano_psel == "Padj") {
    p <- ggplot(res, aes(x = logFC, y =-log10(Adj.P.Value), text=UniqueID))
    ylab <- "-log10(Padj.Value)"
  } else {
    p <- ggplot(res, aes(x = logFC, y =-log10(P.Value), text=UniqueID))
    ylab <- "-log10(P.Value)"
  }
  p <- p	+
    scale_color_manual(values = c("grey", "green2","red2")) +
    geom_point(aes(color = color)) +
    theme_bw(base_size = 20) +
    geom_hline(yintercept = -log10(pvalcut), colour="grey") +
    geom_vline(xintercept = c(-FCcut,0,FCcut), colour="grey") +
    ylab(ylab) + xlab("log2 Fold Change") +
    ggtitle(test_sel)+theme(plot.title = element_text(size = 20),
                            axis.title.x = element_text(size = 14),
                            axis.title.y = element_text(size = 14),
                            legend.text=element_text(size=12))
  p$elementId <- NULL
  
  p <- ggplotly(p) %>% layout(legend = list(orientation = 'h', y=-0.2))
  p
  
})

volcanoplotstatic_out <- reactive({
  res = DatavolcanoReactive()
  DataIn = DataReactive()
  ProteinGeneName = DataIn$ProteinGeneName
  test_sel = input$volcano_test
  FCcut = log2(as.numeric(input$volcano_FCcut))
  FCcut_rd=round(FCcut*1000)/1000
  pvalcut = as.numeric(input$volcano_pvalcut)
  
  if (input$volcano_psel == "Padj") {
    p <- ggplot(res, aes(x = logFC, y = -log10(Adj.P.Value)))
    ylab <- "-log10(Padj.Value)"
    
    filterSig <- paste0("Padj", "<", pvalcut, " & abs(log2FC)>", FCcut_rd)
    data.label <- filter(res, color == filterSig)
    if (nrow(data.label) > input$Ngenes) {
      data.label <- top_n(data.label, input$Ngenes, abs(logFC_ori))
    }
    
  } else {
    filterSig <- paste0("pval", "<", pvalcut, " & abs(log2FC)>",FCcut_rd)
    data.label <- filter(res, color == filterSig)
    if (nrow(data.label) > input$Ngenes) {
      data.label <- top_n(data.label, input$Ngenes, abs(logFC_ori))
    }
    p <- ggplot(res, aes(x = logFC, y = -log10(P.Value)))
    ylab <- "-log10(P.Value)"
  }
  
  if (input$volcano_label=="Upload") {
    req(input$volcano_gene_list)
    volcano_gene_list <- input$volcano_gene_list
    if(grepl("\n",volcano_gene_list)) {
      volcano_gene_list <-  stringr::str_split(volcano_gene_list, "\n")[[1]]
    } else if(grepl(",",volcano_gene_list)) {
      volcano_gene_list <-  stringr::str_split(volcano_gene_list, ",")[[1]]
    }
    volcano_gene_list <- gsub(" ", "", volcano_gene_list, fixed = TRUE)
    volcano_gene_list <- unique(volcano_gene_list[volcano_gene_list != ""])
    uploadlist <- dplyr::filter(ProteinGeneName, (UniqueID %in% volcano_gene_list) | (Protein.ID %in% volcano_gene_list) | (Gene.Name %in% volcano_gene_list))  %>%
      dplyr::select(UniqueID) %>% 	collect %>%	.[["UniqueID"]] %>%	as.character()
    validate(need(length(uploadlist)>0, message = "input gene list"))
    if (length(uploadlist)>input$Ngenes) {uploadlist=uploadlist[1:input$Ngenes]}
    data.label<-res%>%filter(UniqueID %in% uploadlist)
    
  }
  
  if (input$volcano_label=="Geneset") {
    req(input$geneset_list)
    volcano_gene_list <- input$geneset_list
    if(grepl("\n",volcano_gene_list)) {
      volcano_gene_list <-  stringr::str_split(volcano_gene_list, "\n")[[1]]
    } else if(grepl(",",volcano_gene_list)) {
      volcano_gene_list <-  stringr::str_split(volcano_gene_list, ",")[[1]]
    }
    volcano_gene_list <- gsub(" ", "", volcano_gene_list, fixed = TRUE)
    volcano_gene_list <- unique(volcano_gene_list[volcano_gene_list != ""])
    
    uploadlist <- dplyr::filter(ProteinGeneName, (UniqueID %in% volcano_gene_list) | (Protein.ID %in% volcano_gene_list) | (toupper(Gene.Name) %in% toupper(volcano_gene_list)) )  %>%
      dplyr::select(UniqueID) %>% 	collect %>%	.[["UniqueID"]] %>%	as.character()
    validate(need(length(uploadlist)>0, message = "Please select at least one valid gene."))
    if (length(uploadlist)>input$Ngenes) {uploadlist=uploadlist[1:input$Ngenes]}
    data.label<-res%>%filter(UniqueID %in% uploadlist)
  }
  
  
  p <- p	+
    scale_color_manual(values = c("grey", "green2","red2"))
  if (input$rasterize=="Yes") { p<-p+geom_point_rast(aes(color = color), size=0.7, alpha=0.6, na.rm=TRUE, dev="ragg") 
  } else {p<-p+geom_point(aes(color = color), size=0.7) }
  p <- p+
    theme_bw(base_size = 20) +
    geom_hline(yintercept = -log10(pvalcut), colour="grey") +
    geom_vline(xintercept = c(-FCcut,0,FCcut), colour="grey") +
    ylab(ylab) + xlab("log2 Fold Change") +
    ggtitle(test_sel) +
    theme(legend.position = input$vlegendpos, legend.text=element_text(size=input$yfontsize))
  if (input$volcano_label!="None") {
    p=p+geom_text_repel(data = data.label,  aes(label=labelgeneid),	size = input$lfontsize,	box.padding = unit(0.35, "lines"),	point.padding = unit(0.3, "lines"))
  }
  p <- p + guides(color = guide_legend(override.aes = list(alpha = 1, size = 4)))
  return(p)
})

output$volcanoplotstatic <- renderPlot({
  volcanoplotstatic_out()
})

DEG_Compare <- reactive({
  res = DatavolcanoReactive1()
  res2=DatavolcanoReactive2()
  DataIn = DataReactive()
  ProteinGeneName = DataIn$ProteinGeneName  
  test_sel = input$volcano_test1
  test_sel2 = input$volcano_test2	
  FCcut = log2(as.numeric(input$volcano_FCcut))
  pvalcut = as.numeric(input$volcano_pvalcut)
  plotdata=merge(res, res2, by="UniqueID")
  plotdata<-plotdata%>%mutate(color1=paste(Sig.x, Sig.y))
  data.label <-plotdata%>%filter(str_detect(color1, "_sig"))%>%mutate(H_logFC=pmax(abs(logFC.x), abs(logFC.y))) #at least one is sig
  c.res<-cor.test(plotdata$logFC.x, plotdata$logFC.y, method="pearson",use = "complete.obs" )
  cor_string=paste("Pearson corr:", format(c.res$estimate, digits=3), "; p-value:", format(c.res$p.value, digits=3), sep="")
  #browser() #debug
  if (nrow(data.label) > input$Ngenes) {data.label <- top_n(data.label,input$Ngenes,H_logFC) }
  
  if (input$volcano_psel == "Padj") {
    p<-ggplot(plotdata, aes(x=logFC.x, y=logFC.y, color=color1,
                            size=-pmin(log10(Adj.P.Value.x),log10(Adj.P.Value.y))))
    if (input$rasterize=="Yes") { p<-p+geom_point_rast(na.rm=TRUE, dev="ragg") 
    } else {p<-p+geom_point()}
    p<-p+
      theme_bw(base_size = 20) + ylab(str_c("log2FC in ", test_sel2)) + xlab(str_c("log2FC in ", test_sel))+
      labs(color='Significance',size='-log10 min Adj.P.Value',title=cor_string) 
    
  } else {
    p<-ggplot(plotdata, aes(x=logFC.x, y=logFC.y, color=color1,
                            size=-pmin(log10(P.Value.x),log10(P.Value.y)))) 
    if (input$rasterize=="Yes") { p<-p+geom_point_rast(na.rm=TRUE, dev="ragg") 
    } else {p<-p+geom_point()}
    p<-p+
      geom_point_rast(na.rm=TRUE, dev="ragg")+
      theme_bw(base_size = 20) + ylab(str_c("log2FC in ", test_sel2)) + xlab(str_c("log2FC in ", test_sel))+
      labs(color='Significance',size='-log10 min P.Value', title=cor_string) 
  }
  if (input$volcano_label=="Upload") {
    req(input$volcano_gene_list)
    volcano_gene_list <- input$volcano_gene_list
    if(grepl("\n",volcano_gene_list)) {
      volcano_gene_list <-  stringr::str_split(volcano_gene_list, "\n")[[1]]
    } else if(grepl(",",volcano_gene_list)) {
      volcano_gene_list <-  stringr::str_split(volcano_gene_list, ",")[[1]]
    }
    volcano_gene_list <- gsub(" ", "", volcano_gene_list, fixed = TRUE)
    volcano_gene_list <- unique(volcano_gene_list[volcano_gene_list != ""])
    uploadlist <- dplyr::filter(ProteinGeneName, (UniqueID %in% volcano_gene_list) | (Protein.ID %in% volcano_gene_list) | (Gene.Name %in% volcano_gene_list))  %>%
      dplyr::select(UniqueID) %>% 	collect %>%	.[["UniqueID"]] %>%	as.character()
    validate(need(length(uploadlist)>0, message = "Please enter at least one valid gene."))
    if (length(uploadlist)>input$Ngenes) {uploadlist=uploadlist[1:input$Ngenes]}
    data.label<-plotdata%>%filter(UniqueID %in% uploadlist)
    
  }
  
  if (input$volcano_label=="Geneset") {
    req(input$geneset_list)
    volcano_gene_list <- input$geneset_list
    if(grepl("\n",volcano_gene_list)) {
      volcano_gene_list <-  stringr::str_split(volcano_gene_list, "\n")[[1]]
    } else if(grepl(",",volcano_gene_list)) {
      volcano_gene_list <-  stringr::str_split(volcano_gene_list, ",")[[1]]
    }
    volcano_gene_list <- gsub(" ", "", volcano_gene_list, fixed = TRUE)
    volcano_gene_list <- unique(volcano_gene_list[volcano_gene_list != ""])
    uploadlist <- dplyr::filter(ProteinGeneName, (UniqueID %in% volcano_gene_list) | (Protein.ID %in% volcano_gene_list) | 
                                  (toupper(Gene.Name) %in% toupper(volcano_gene_list)))  %>%
      dplyr::select(UniqueID) %>% 	collect %>%	.[["UniqueID"]] %>%	as.character()
    validate(need(length(uploadlist)>0, message = "Please select at least one valid gene."))
    if (length(uploadlist)>input$Ngenes) {uploadlist=uploadlist[1:input$Ngenes]}
    data.label<-plotdata%>%filter(UniqueID %in% uploadlist)
  }
  
  
  p<-p+ scale_color_manual(values=c('X_sig Y_sig'='blue3','X_sig Y_notsig'='green3',
                                    'X_notsig Y_sig'='orange','X_notsig Y_notsig'='#00000022')) + 
    theme(legend.position=input$vlegendpos, legend.text=element_text(size=input$yfontsize), legend.title=element_text(size=input$yfontsize+1))
  
  if (input$volcano_label=="Upload" || input$volcano_label=="Geneset" ) {
    p=p+ geom_text_repel(data = data.label,  aes(label=labelgeneid.x),	size = input$lfontsize,	box.padding = unit(0.35, "lines"),	
                         color="coral3",  point.padding = unit(0.3, "lines"))
  } #uploaded list use a different color. The DEGs colors are hard to see for un-sig genes.
  if (input$volcano_label=="DEGs") {
    if (input$DEG_comp_color=="No") {
      p=p+ geom_text_repel(data = data.label, color="coral3",  aes(label=labelgeneid.x),	size = input$lfontsize,	box.padding = unit(0.35, "lines"),	
                           point.padding = unit(0.3, "lines"))
    } else {
      p=p+ geom_text_repel(data = data.label, aes(label=labelgeneid.x),	size = input$lfontsize,	box.padding = unit(0.35, "lines"),	
                           point.padding = unit(0.3, "lines"))
    }
  }  
  #browser() #debug
  if (input$DEG_comp_XY=="Yes"){
    XY_min=min(min(plotdata$logFC.x), min(plotdata$logFC.y))
    XY_max=max(max(plotdata$logFC.x), max(plotdata$logFC.y))
    p<-p+xlim(XY_min, XY_max)+ylim(XY_min, XY_max)
  }
  # Remove strange "a" from legend
  p <- p + guides(color = guide_legend(override.aes = list(alpha = 1, size = 4)))
  return(p)
})

output$DEG_Compare <- renderPlot({
  DEG_Compare()
})

observeEvent(input$volcano, {
  test_sel = input$volcano_test
  saved_plots$volcano[[test_sel]] <- volcanoplotstatic_out()
})

observeEvent(input$DEG_comp, {
  test_sel = paste(input$volcano_test1, "vs", input$volcano_test2)
  saved_plots$volcano[[test_sel]] <- DEG_Compare()
})

DEG_data <-reactive ({
  DataIn = DataReactive()
  results_long = DataIn$results_long
  test_sel = input$volcano_test
  FCcut = log2(as.numeric(input$volcano_FCcut))
  pvalcut = as.numeric(input$volcano_pvalcut)
  if (input$volcano_psel == "Padj") {
    tmpdat = results_long %>% filter(test==test_sel & Adj.P.Value < pvalcut & abs(logFC) > FCcut) 
  } else {
    tmpdat = results_long %>% filter(test==test_sel & P.Value < pvalcut & abs(logFC) > FCcut) 
  }
  tmpdat[,sapply(tmpdat,is.numeric)] <- signif(tmpdat[,sapply(tmpdat,is.numeric)],3)
  return(tmpdat)
})
output$volcanoData <- DT::renderDataTable({
  DT::datatable(DEG_data(),extensions = 'Buttons',  options = list(
    dom = 'lBfrtip', buttons = c('csv', 'excel', 'print'), pageLength = 20), rownames= FALSE)
})

observeEvent(input$DEG_data, {
  saved_table$DEG_data <- DEG_data()
})

