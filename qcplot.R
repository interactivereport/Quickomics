###########################################################################################################
## Proteomics Visualization R Shiny App
##
##This software belongs to Biogen Inc. All right reserved.
##
##@file: qcplot.R
##@Developer : Benbo Gao (benbo.gao@Biogen.com)
##@Date : 5/16/2018
##@version 1.0
###########################################################################################################


# observe({
#   #DataIn = DataReactive()
#   MetaData=all_metadata()
#   attributes=sort(setdiff(colnames(MetaData), c("sampleid", "Order", "ComparePairs") ))
#   updateSelectInput(session, "PCAcolorby", choices=attributes, selected="group")
#   updateSelectInput(session, "PCAshapeby", choices=c("none", attributes), selected="none")
#   updateSelectInput(session, "PCAsizeby", choices=c("none", attributes), selected="none")
#   attrs=sort(setdiff(colnames(MetaData), c("Order", "ComparePairs") ))
#   updateRadioButtons(session,'PCA_label', inline = TRUE, choices=attrs, selected="sampleid")
#   updateTextInput(session, "Ylab", value=exp_unit())
#   if (!is.null(MetaData)) {
#     if (nrow(MetaData)>100) {updateRadioButtons(session,'PCA_subsample', selected="None")} #when there are too many samples, don't show  labels
#   }
# })

observeEvent(all_metadata(), {
  MetaData <- all_metadata()
  attributes <- sort(setdiff(colnames(MetaData), c("Order", "ComparePairs")))
  updateSelectInput(session, "PCAcolorby", choices=attributes, selected="group")
  updateSelectInput(session, "PCAshapeby", choices=c("none", attributes), selected="none")
  updateSelectInput(session, "PCAsizeby", choices=c("none", attributes), selected="none")
  attrs=sort(setdiff(colnames(MetaData), c("Order", "ComparePairs") ))
  updateRadioButtons(session,'PCA_label', inline = TRUE, choices=attrs, selected="sampleid")
  updateTextInput(session, "Ylab", value=exp_unit())
  if (!is.null(MetaData)) {
    if (nrow(MetaData)>100) {updateRadioButtons(session,'PCA_subsample', selected="None")} #when there are too many samples, don't show  labels
  }
})

observe({
  samples <- sample_order()
  updateTextAreaInput(session, "PCA_list", value=paste(samples, collapse="\n"))
})


####################################################
output$selectGroupSampleQC <- renderUI(shared_header_content())

DataPCAReactive <- reactive({
  browser()
  DataIn <-  DataQCReactive()
  tmp_sampleid <- DataIn$tmp_sampleid
  validate(need(length(tmp_sampleid)>1, message = "Please select at least two samples (please note samples are filtered by group selection as well)."))
  
  MetaData=DataIn$MetaData
  tmp_group = MetaData$group
  
  tmp_data_wide <- DataIn$tmp_data_wide
  #clean up data
  tmp_data_wide=na.omit(tmp_data_wide)
  tmp_data_wide <- tmp_data_wide[apply(tmp_data_wide, 1, sd) != 0, ] #remove rows with all 0s, or the same value across all samples   
  
  pca <- 	prcomp(t(tmp_data_wide),rank. = 10, scale = TRUE)
  percentVar <- 	round((pca$sdev)^2/sum(pca$sdev^2), 3) * 100
  scores <- as.data.frame(pca$x)
  rownames(scores) <- tmp_sampleid
  all_groups <- DataIn$tmp_group$group
  scores$group <- factor(tmp_group, levels = all_groups)
  attributes=setdiff(colnames(MetaData), c("Order", "ComparePairs", "group") )
  
  colsel=match(attributes, colnames(MetaData) )
  scores=cbind(scores, MetaData[, colsel, drop=F])
  for (attr in attributes) {
    if (attr %in% names(DataIn$tmp_group)) {
      scores[[attr]] <- factor(
        scores[[attr]],
        levels = DataIn$tmp_group[[attr]]
      )
    }
  }
  #browser() #debug
  return(list('scores'=scores,'percentVar'=percentVar))
})

#Eigenvalue bar chart
Eigenvalues_plot<-reactive({
  req(DataPCAReactive())
  PCAlist <- DataPCAReactive()
  scores <- PCAlist$scores
  percentVar <- PCAlist$percentVar
  plotdata<-data.frame(PC=names(scores)[1:10], perVar=percentVar[1:10])
  plotdata$PC=factor(plotdata$PC, levels=plotdata$PC)
  plotdata<-plotdata%>%mutate(TotalVar=cumsum(perVar))
  adj.factor=max(plotdata$TotalVar)/max(plotdata$perVar)*0.9
  p<-ggplot(plotdata, aes(x=PC) )+geom_bar(aes(y=perVar), stat="identity", fill="blue4")+
    geom_line(aes(y=TotalVar/adj.factor), size=1.5, color="red4", group=1)+geom_point(aes(y=TotalVar/adj.factor), size=3, color="red4")+
    labs(x="Principal Components")+scale_y_continuous(name="Percentage of Variance", sec.axis=sec_axis(~.*adj.factor, name="Total Variance") ) +theme_cowplot()
  return(p)
})
output$Eigenvalues <- renderPlot({
  Eigenvalues_plot()
})


########## boxplot
QCboxplot_out <- reactive({
	withProgress(message = 'Making box plot', value = 0, {
	  
		DataQC <-  DataQCReactive()
		tmp_sampleid <- DataQC$tmp_sampleid
		#tmp_data_long <- DataQC$tmp_data_long %>% dplyr::filter(expr !=0) %>% sample_n(1000)
		tmp_data_long<-DataQC$tmp_data_long%>% group_by(sampleid) %>% dplyr::slice_sample(n=2000) %>% ungroup #max 2K genes/proteins per sample	
		tmp_data_long <- tmp_data_long %>%
		  mutate(sampleid = factor(sampleid, levels = tmp_sampleid))
		
		p <- ggplot(tmp_data_long, aes(x=sampleid, y=expr)) +
		geom_boxplot(aes(color=factor(sampleid)), outlier.colour = NA) +

		#scale_fill_manual(values=rep("Dark2", length(tmp_sampleid)))+

		coord_cartesian(ylim = range(boxplot(tmp_data_long$expr, plot=FALSE)$stats)*c(.9, 1.2)) +
		labs(x = "Sample", y = exp_unit()) +
		theme_bw(base_size = 20) +
		theme(legend.position = "bottom",	legend.title=element_blank(),	axis.text.x = element_blank(), plot.margin=unit(c(1,1,1,1),"mm")) +
		guides(col = guide_legend(ncol = 8))
		return(p)
	}
)
})

output$QCboxplot <- renderPlot({
	QCboxplot_out()
})

observeEvent(input$QCboxplot, {
	saved_plots$QCboxplot <- QCboxplot_out()
})


######## PCA
observeEvent(input$plot_PCA, {
  plot_pca_control(plot_pca_control()+1)
})

pcaplot_out <- eventReactive (plot_pca_control(), {
  browser()
  ptm <- proc.time()
  req(DataPCAReactive())
  req(input$PCA_label != "")
  pcnum=as.numeric(input$pcnum)
	validate(need(length(pcnum)==2, message = "Select 2 Prinical Components."))

	#DataQC <-  DataQCReactive()
	#tmp_group = DataQC$tmp_group
	
	PCAlist <- DataPCAReactive()
	scores <- PCAlist$scores
	percentVar <- PCAlist$percentVar
	samples=scores$sampleid

	xlabel <- paste("PC",pcnum[1],"(",round(percentVar[pcnum[1]]),"%)",sep="")
	ylabel <- paste("PC",pcnum[2],"(",round(percentVar[pcnum[2]]),"%)",sep="")

	PC1 <- paste("PC",pcnum[1],sep="")
	PC2 <- paste("PC",pcnum[2],sep="")
	n <- length(unique(as.character(unlist(scores[, colnames(scores)==input$PCAcolorby]))))
	colorpal=get_pal_ramp(input$PCAcolpalette, n)

	#if (all(table(tmp_group))<4)
	#  ellipsoid = FALSE 

	if (input$PCA_subsample=="None" ) {
	  labels=NULL
	} else {
	  label_sel=match(input$PCA_label, names(scores))
	 # browser() #debug
	  labels=unlist(scores[, label_sel])	
	  if (input$PCA_subsample=="Subset") {
	    PCA_list=str_split(input$PCA_list, "\n")[[1]]
	    N_sel=match(PCA_list, samples)
	    N_sel=N_sel[!is.na(N_sel)]
	    validate(need(length(N_sel)>0, message = "Enter at least one valid sampleid to label"))
	    keep_s=rep(FALSE, length(labels))
	    keep_s[N_sel]=TRUE
	    labels[!keep_s]=""
	    #browser() #debug
	  }
	}

	if (input$PCAshapeby=="none") {shape_by=19} else {shape_by=input$PCAshapeby}
	if (input$PCAsizeby=="none") {size_by=input$PCAdotsize} else {size_by=input$PCAsizeby}	
	if (is.numeric(scores[[input$PCAcolorby]])) {  #when colorby is numeric, don't use color palette
	  p <- ggpubr::ggscatter(scores,x =PC1, y=PC2, color =input$PCAcolorby, shape=shape_by, size =size_by , ellipse = input$ellipsoid, mean.point = input$mean_point, rug = input$rug,
	                         label =labels, font.label = input$PCAfontsize, repel = TRUE,  ggtheme = theme_bw(base_size = 20) )
	} else {
	  p <- ggpubr::ggscatter(scores,x =PC1, y=PC2, color =input$PCAcolorby, shape=shape_by, size =size_by , palette= colorpal, ellipse = input$ellipsoid, mean.point = input$mean_point, rug = input$rug,
	                         label =labels, font.label = input$PCAfontsize, repel = TRUE,  ggtheme = theme_bw(base_size = 20) )
	}

	p <- ggpubr::ggpar(p, xlab = xlabel, ylab = ylabel)
	#	browser() #debug	
	#	p <- ggpubr::ggpar(p, legend.title ="", xlab = xlabel, ylab = ylabel, legend = "bottom") #works only when use color by. 
	p <- p + guides(color = guide_legend(override.aes = list(label="")))
	#cat("generated pca plot p",(proc.time() - ptm)[["elapsed"]], "\n")
	return(p)
})

output$pcaplot <- renderPlot({
  ptm <- proc.time()
  withProgress(message = 'Drawing PCA Plot...', value = 0, {
	print(pcaplot_out()) })
	cat("plotted PCA",(proc.time() - ptm)[["elapsed"]], "\n")
})

observeEvent(input$pcaplot, {
	saved_plots$pcaplot <- pcaplot_out()
}
)

output$pca_legend <- renderPlot({
  PCAlist <- DataPCAReactive()
  scores <- PCAlist$scores
  color_by=input$PCAcolorby
  tmp_group=as.character(unlist(scores[, colnames(scores)==color_by]))
  n <- length(unique(tmp_group))
  colorpal = get_pal_ramp(input$PCAcolpalette, n)
  tmp_plot<-ggplot(scores, aes_string(x="PC1", y="PC2", color=color_by))+geom_point()+scale_color_manual(values=colorpal)+ theme_cowplot(12)
  legend_only <- get_legend(tmp_plot +theme(legend.position = "bottom",  legend.title = element_text(size = 16),
                                            legend.text = element_text(size = 14))+guides(color = guide_legend(override.aes = list(size=8))))
  plot_grid(legend_only)
})

######## PCA 3D
output$plot3d <- renderRglwidget({
 	PCAlist <- DataPCAReactive()
	scores <- PCAlist$scores
	percentVar <- PCAlist$percentVar

	xlabel <- paste("PC1(",round(percentVar[1]),"%)",sep="")
	ylabel <- paste("PC2(",round(percentVar[2]),"%)",sep="")
	zlabel <- paste("PC3(",round(percentVar[3]),"%)",sep="")

	sampleid <- rownames(scores)
	
	tmp_group=as.character(unlist(scores[, colnames(scores)==input$PCAcolorby]))
	n <- length(unique(tmp_group))
	#colorpal = topo.colors(n, alpha = 1)
	#colorpal = get_palette("Dark2", n)
	colorpal = get_pal_ramp(input$PCAcolpalette, n)
	scores$tmp_group=unlist(scores[, colnames(scores)==input$PCAcolorby])

	
	#rgl.open(useNULL=T)
	options(rgl.useNULL=TRUE)
	if (input$ellipsoid3d == "Yes") {
	  ellipsoid3d = TRUE
	} else {
	  ellipsoid3d = FALSE 
	}


	if (any(table(tmp_group) <= 3))
	  ellipsoid3d = FALSE 
	

	
	if (input$dotlabel == "Yes") {
	  dotlabel=TRUE
	} else {
	  dotlabel=FALSE
	}

	scatter3d(PC3 ~ PC1 + PC2 | tmp_group, data= scores,
	          axis.col= c("black", "black", "black"),
	          xlab=xlabel, ylab=ylabel,  zlab=zlabel, labels = as.factor(sampleid), id=dotlabel, id.n=length(sampleid),
	          axis.scales=FALSE,  axis.ticks=FALSE,
	          ellipsoid = ellipsoid3d,
	          surface=FALSE, grid = FALSE,
	          cex.lab=3,
	          surface.col = colorpal)
	rglwidget(width = 800, height = 800)
})

output$plotly3d <- renderPlotly({
	PCAlist <- DataPCAReactive()
	scores <- PCAlist$scores
	scores<-scores%>%mutate_if(is_character, as.factor)
	percentVar <- PCAlist$percentVar
	symbol_list=rep(c('circle', 'square',  'diamond',  'circle-open','square-open','diamond-open'), 2) #symbols which work with plotly scatter3d
	plot_symbols=symbol_list[unique(as.numeric(unlist(scores[, colnames(scores)==input$PCAshapeby])))]
	
	xlabel <- paste("PC1(",round(percentVar[1]),"%)",sep="")
	ylabel <- paste("PC2(",round(percentVar[2]),"%)",sep="")
	zlabel <- paste("PC3(",round(percentVar[3]),"%)",sep="")

	sampleid <- str_c(scores$sampleid, "\n", scores$group)
	n <- length(unique(as.character(unlist(scores[, colnames(scores)==input$PCAcolorby]))))
	colorpal = get_pal_ramp(input$PCAcolpalette, n)
	if (input$PCAshapeby=="none"){
	  p <- plot_ly(scores, x = ~PC1, y = ~PC2, z = ~PC3, color = as.formula(paste0("~", input$PCAcolorby)), 
	               colors = colorpal,text = sampleid) %>%
	    add_markers() %>%
	    layout(scene = list(xaxis = list(title = xlabel), yaxis = list(title = ylabel),  zaxis = list(title = zlabel)))
	  
	} else{
	p <- plot_ly(scores, x = ~PC1, y = ~PC2, z = ~PC3, color = as.formula(paste0("~", input$PCAcolorby)), 
	              symbol=as.formula(paste0("~", input$PCAshapeby)),symbols=plot_symbols, 
	             colors = colorpal,text = sampleid) %>%
	add_markers() %>%
	layout(scene = list(xaxis = list(title = xlabel), yaxis = list(title = ylabel),  zaxis = list(title = zlabel)))
	}
	p$elementId <- NULL
	p
})

############heatmap
pheatmap_out <- reactive({
	DataQC <-  DataQCReactive()
	tmp_sampleid <- DataQC$tmp_sampleid
	tmp_data_wide <- DataQC$tmp_data_wide
	MetaData=DataQC$MetaData
	
	selCol=which(names(MetaData)==input$PCAcolorby)
	annotation=MetaData[, selCol, drop=F]
	#annotation = data.frame("group" = tmp_group)
	rownames(annotation) <- tmp_sampleid

	sampleDistMatrix <- as.matrix(dist(t(tmp_data_wide)))
	rownames(sampleDistMatrix) <- tmp_sampleid
	colors <- colorRampPalette(rev(brewer.pal(9, "Blues")) )(32)

	p <- pheatmap::pheatmap(sampleDistMatrix, annotation_row=annotation,	annotation_col=annotation, col=colors)
	return(p)
})

output$pheatmap <- renderPlot({
	grid.draw(pheatmap_out()$gtable)
})

observeEvent(input$SampleDistance, {
	saved_plots$SampleDistance <- pheatmap_out()$gtable
}
)

############Dendrograms
Dendrograms_out <- reactive({
	hc <- pheatmap_out()$tree_row
	if (input$dendroformat=="tree") {
		p <- fviz_dend(hc, k = input$DendroCut, cex = input$DendroFont, k_colors = "jco",	color_labels_by_k = TRUE, rect = TRUE, rect_border = "jco",	rect_fill = TRUE)
	} else 	if (input$dendroformat=="horiz") {
		p <- fviz_dend(hc, k = input$DendroCut, cex = input$DendroFont, k_colors = "jco",	 horiz = TRUE, color_labels_by_k = TRUE, rect = TRUE, rect_border = "jco", rect_fill = TRUE)
	} else if (input$dendroformat=="circular") {
		p <- fviz_dend(hc, k = input$DendroCut, cex = input$DendroFont, k_colors = "jco", type = "circular")
	}
	return(p)
})

output$Dendrograms <- renderPlot({
	Dendrograms_out()
})

observeEvent(input$Dendrograms, {
	saved_plots$Dendrograms <- Dendrograms_out()
})

############histplot
histplot_out <- reactive({
  withProgress(message = 'Calculating.',  detail = 'This may take a while...', value = 0, {
	DataQC <-  DataQCReactive()
	tmp_sampleid <- DataQC$tmp_sampleid
	tmp_data_long <- DataQC$tmp_data_long
	# tmp_group = DataQC$tmp_group
	if (!"id" %in% names(tmp_data_long)) {tmp_data_long$id=tmp_data_long$UniqueID} 
	#browser()
	# CV.df <- tmp_data_long %>%
	# group_by(.,  group, id) %>%
	# dplyr::summarise( mean=mean(expr, na.rm = TRUE), sd=sd(expr, na.rm = TRUE)) %>%
	# dplyr::mutate(CV=100*(sd/mean))
	
	CV.df <- tmp_data_long %>%
	  dplyr::group_by(group, id) %>%
	  dplyr::summarise(mean = mean(expr, na.rm = TRUE), sd = sd(expr, na.rm = TRUE)) %>%
	  mutate(
	    sd = ifelse(is.na(sd), 0, sd),
	    CV = 100 * sd / mean
	  )
	
	mu <- group_by(CV.df,group) %>%
	dplyr::summarise(median = round(median(CV, na.rm = TRUE),1))

	interval <- seq.int(0, 100, 5)
	xlimmin <- interval[cut(min(mu$median), interval, include.lowest = TRUE, labels = FALSE)]
	xlimmax <- interval[cut(max(mu$median), interval, include.lowest = TRUE, labels = FALSE) +1]
	p <- ggplot(CV.df, aes(x=CV, color=group)) +
	geom_freqpoly (position="dodge", na.rm = TRUE, bins = 10) +
	geom_vline(data=mu, aes(xintercept=median, color=group), linetype="dashed") +
	geom_text(data=mu, mapping=aes(x=median, y=0, label=paste(median,"(",group,")", sep="")), size=4, angle=90, vjust=-0.4, hjust=0) +
	scale_x_continuous(breaks = seq(xlimmin, xlimmax, by=5), limits=c(xlimmin,xlimmax)) +
	theme_bw(base_size = 20) +
	theme(legend.position = "bottom")
	return(p)
  })
})

output$histplot <- renderPlot({
	histplot_out()
})

observeEvent(input$histplot, {
	saved_plots$histplot <- histplot_out()
})


############PC_covariates QC Plots
PC_covariates_out <-  eventReactive(input$compute_PC,{
  DataQC <-  DataQCReactive()
  tmp_data_wide <- DataQC$tmp_data_wide
  MetaData=DataQC$MetaData
  meta=MetaData[, !(colnames(MetaData) %in% c("sampleid", "Order", "ComparePairs")), drop=FALSE]
  meta=meta[, (colnames(meta) %in% input$covar_variates), drop=FALSE]
  rownames(meta)=MetaData$sampleid
  #browser() #debug
  #clean up data
	tmp_data_wide=na.omit(tmp_data_wide)
	tmp_data_wide <- tmp_data_wide[apply(tmp_data_wide, 1, sd) != 0, ] #remove rows with all 0s, or the same value across all samples 
  res<-Covariate_PC_Analysis(tmp_data_wide, meta, out_prefix=NULL, PC_cutoff=input$covar_PC_cutoff, 
            FDR_cutoff=input$covar_FDR_cutoff, N_col=input$covar_ncol)
  #print(res$PC_info) #debug, should match main PCA results
  return(res)
})

#output$covar_table=renderTable(PC_covariates_out()$selVar_All, colnames=T)

output$covar_table <- DT::renderDT(server=FALSE,{
  results<-PC_covariates_out()$selVar_All
  if (!is.null(results)) {
    results["P-value"]=as.numeric(formatC(unlist(results["P-value"]), format="e", digits=2))
    results["FDR"]=as.numeric(formatC(unlist(results["FDR"]), format="e", digits=2))
  }
  DT::datatable(results,  extensions = 'Buttons',
                options = list(
                  dom = 'lBfrtip', pageLength = 25,
                  buttons = list(
                    list(extend = "csv", text = "Download Page", filename = "Page_results",
                         exportOptions = list(modifier = list(page = "current"))),
                    list(extend = "csv", text = "Download All", filename = "All_Results",
                         exportOptions = list(modifier = list(page = "all")))
                  )
                ),rownames= T)
})


output$PC_covariatesC <- renderPlot({
  data=PC_covariates_out()$sel_dataC
  if (!is.null(data)) {
    data$plot
  }
})
output$plot.PC_covariatesC=renderUI({
  tagList(
    textOutput("N_pairs_C"),
    plotOutput("PC_covariatesC",height = input$covar_cat_height)
  )
}) 


output$PC_covariatesN <- renderPlot({
  data=PC_covariates_out()$sel_dataN
  if (!is.null(data)) {
    data$plot
  }
})
output$plot.PC_covariatesN=renderUI({
  tagList(
  textOutput("N_pairs_N"),
  plotOutput("PC_covariatesN",height = input$covar_num_height)
  )
})

Npairs_cov<-reactive({
  res<-PC_covariates_out()
  C=res$sel_dataC$selVar
  if (is.null(C)) {N1=0} else {N1=nrow(C)}
  N=res$sel_dataN$selVar
  if (is.null(N)) {N2=0} else {N2=nrow(N)}
  return(c(N1, N2))
})

observe({
  H_C=ceiling(Npairs_cov()[1]/PC_covariates_out()$ncol)*400
  if (H_C>0)  { updateSliderInput(session, "covar_cat_height", value = H_C)}
  H_N=ceiling(Npairs_cov()[2]/PC_covariates_out()$ncol)*400
  if (H_N>0)  { updateSliderInput(session, "covar_num_height", value = H_N)}  
})

output$N_pairs_C<-renderText({str_c("There are ", Npairs_cov()[1], " significant categorical covariate-PC pairs.")})
output$N_pairs_N<-renderText({str_c("There are ", Npairs_cov()[2], " significant numeric covariate-PC pairs.")})
output$N_pairs<-renderText({str_c("There are ", Npairs_cov()[1]+Npairs_cov()[2], " significant covariate-PC pairs.")})


observeEvent(input$covar_cat, {
  data=PC_covariates_out()$sel_dataC
  saved_plots$covar_cat <- data$plot
})

observeEvent(input$covar_num, {
  data=PC_covariates_out()$sel_dataN
  saved_plots$covar_num<- data$plot
})
