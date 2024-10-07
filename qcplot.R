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


observe({
	#DataIn = DataReactive()
	MetaData=all_metadata()
	attributes=sort(setdiff(colnames(MetaData), c("sampleid", "Order", "ComparePairs") ))
	updateSelectInput(session, "PCAcolorby", choices=attributes, selected="group")
	updateSelectInput(session, "PCAshapeby", choices=c("none", attributes), selected="none")
	updateSelectInput(session, "PCAsizeby", choices=c("none", attributes), selected="none")
	sampleIDs=sort(setdiff(colnames(MetaData), c("Order", "ComparePairs") ))
	updateRadioButtons(session,'PCA_label', inline = TRUE, choices=sampleIDs, selected="sampleid")
	updateSelectInput(session, "covar_variates", choices=attributes, selected=attributes)
	updateTextInput(session, "Ylab", value=exp_unit())
	if (!is.null(MetaData)) {
	  if (nrow(MetaData)>100) {updateRadioButtons(session,'PCA_subsample', selected="None")} #when there are too many samples, don't show  labels
	}
})

observe({
  #DataIn = DataReactive()
  allsamples = all_samples()
  allgroups = all_groups()
  samples <- sample_order()
  groups = group_order()
  updateSelectizeInput(session,'QC_groups', choices=allgroups, selected=groups)
  updateSelectizeInput(session,'QC_samples', choices=allsamples, selected=samples)
  updateTextAreaInput(session, "PCA_list", value=paste(samples, collapse="\n"))
})

output$reorder_group=renderUI({
  req(group_order())
  orderInput(inputId = 'order_groups', label = 'Drag and Drop to Reorder Groups. (Use Select Groups at left menu to delete or add groups.)', items =group_order(), width="90%", item_class = 'primary', legacy =TRUE )
})

attribute_filters_text=reactive({
  if (!is.null(attribute_filters())) {
    text=""
    attr=attribute_filters()
    #browser() #debug
    for (i in 1:length(attr) ) {
      m1=attr[[i]]
      name1=names(attr[i])
      text=str_c(text," (", name1, ":  ", paste(m1, collapse=","), ")" )
    }
  } else (text="None (Filters that select 0 samples will be reset)")
  return(text)  
})

output$sample_choose_order=renderUI({
  req(group_order())
  req(sample_order())
  group_exclude<-setdiff(all_groups(), group_order())
  sample_exclude<- setdiff(all_samples(), sample_order())
  MetaData=all_metadata()
  attributes=c("None",  setdiff(colnames(MetaData), c("sampleid", "Order", "ComparePairs", "group") ) )
  attributes=sort(attributes)
  #browser() #debug
  tagList(
    tags$div(
      tags$p("Groups excluded: ", paste(group_exclude, collapse=", "), " (samples from the excluded groups are removed)")),
      tags$p(tags$em("Add/remove/re-roder groups will reset sample selection based on slected groups.")),
    tags$hr(style="border-color: RoyalBlue;"),
    tags$p(tags$strong("Please finalized group selection/order before working on further sample selection.")),
    radioButtons("Select_Sample", label="Method to Select Samples:", inline = TRUE, choices = c("Sample Filter", "Upload Sample List", "From Comparison", "None"), selected = "Sample Filter"),
    conditionalPanel("input.Select_Sample=='Sample Filter'",
      tags$p("Sample attribute filters already applied: ", attribute_filters_text() ),
      tags$p("Manualy removed samples: ", paste(samples_excludeM(), collapse=",") ),
      selectizeInput("meta_col_sel", label="Filter Samples by the Attribute (Covariate) Below:", choices=attributes, selected="None", multiple=FALSE),
      conditionalPanel(condition="input.meta_col_sel!='None'",
                     uiOutput('filter_meta')),
      checkboxInput("remove_samples", "Remove additonal samples?", FALSE, width="90%"),
      conditionalPanel(condition="input.remove_samples==1",
                     textAreaInput("sample_exclude_list", "Enter Samples to Exclude:", "",  width="500px", height="50px"),
                     actionButton("remove_sample", "Remove Samples in the Box Above")) ),
    conditionalPanel("input.Select_Sample=='Upload Sample List'",
      textAreaInput("sample_upload_list", "Enter Samples to Use (IDs separated by comma or line break):", "",  width="500px", height="150px"),
      actionButton("upload_sample", "Upload Samples in the Box Above")),
    conditionalPanel("input.Select_Sample=='From Comparison'",
      tags$p("This tool is only available when there is comparison information (comp_info in RData) with valid Subsetting_group."),
      tags$p("To select multiple comparisons, use the left menu tool 'Use Samples in Subset/Comparison'."),
      tags$hr(),
      uiOutput("samples_from_comp")),  
    tags$hr(),
    tags$p("Selected samples: ", paste(sample_order(), collapse=", ")),
    tags$p("Excluded samples: ", paste(sample_exclude, collapse=", ")),
    tags$br(),
    tags$hr(style="border-color: RoyalBlue;"),
    checkboxInput("show_samples", "Reorder Selected samples?", TRUE, width="90%"),
    conditionalPanel(condition="input.show_samples==1",
                     orderInput(inputId = 'order_samples', label = 'Drag and Drop to Reorder Samples.', items =sample_order(), width="90%", item_class = 'success', legacy =TRUE ))
  )
})

output$samples_from_comp=renderUI({
  sel_comp=DataReactive()$sel_comp
  if (!is.null(sel_comp)) {
    comp_all=c("None", "All_Samples",sel_comp$Comparison)
    sel_comp=sel_comp[, 1:(ncol(sel_comp)-2)]
    sel_comp$N_samples=as.integer(sel_comp$N_samples)
    sel_comp=sel_comp[, c(1, ncol(sel_comp), 2:(ncol(sel_comp)-1))]
    subsets=unique(sel_comp$Subsetting_group)
    sel_sub=str_detect(subsets, ":"); subsets=subsets[sel_sub]
    if (length(subsets)>0) {comp_all=c(comp_all, str_c("(Subset) ", subsets))}
    output$table_sel_comp=renderTable(sel_comp, rownames=F, colnames=T)
    tagList(
      selectizeInput("comp_4_samples", label="Subset Samples Based on a Comparison:", choices=comp_all, selected="None", multiple=FALSE),
      tableOutput("table_sel_comp"),
      HTML("<hr>")
    )
  }
})

output$QC_samples_from_comp<-renderUI({
  sel_comp=DataReactive()$sel_comp
  if (!is.null(sel_comp)) {
    comp_all=c("None","All_Samples", sel_comp$Comparison)
    subsets=unique(sel_comp$Subsetting_group)
    sel_sub=str_detect(subsets, ":"); subsets=subsets[sel_sub]
    if (length(subsets)>0) {comp_all=c(comp_all, str_c("(Subset) ", subsets))}
    tagList(
      selectizeInput("QC_comp_4_samples", label="Use Samples in Comparison:", choices=comp_all, selected="None", multiple=TRUE),
    )
  } else {
    tagList(
      tags$p("No comparison-sample information. (Need comp_info in RData file)."))
  }
})


output$filter_meta=renderUI({
  MetaData=all_metadata()
  selCol=input$meta_col_sel
  if (selCol!="None") {
  values=as.character(MetaData[[selCol]])
  UniqueValues=sort(unique(values))
  UniqueValues[UniqueValues==""]="Empty_Value" #change so empty values can be displayed and selected
  group_info<-data.frame(values)%>%group_by(values)%>%dplyr::count()%>%t()
  output$table_selCol=renderTable(group_info, colnames=F)
  sel_row=which(MetaData$sampleid %in% sample_order())
  values2=values[sel_row]
  group_info2<-data.frame(values2)%>%group_by(values2)%>%dplyr::count()%>%t()
  output$table_selCol2=renderTable(group_info2, colnames=F)
  menu_label=str_c("Select Values from ", selCol, " (ctrl & click to select multiple values to remove)")
  #browser() ##debug
  #check if existing filter cover the selected columns
  selValues=UniqueValues
  attr=attribute_filters()
  if (!is.null(attr)) {
    for (i in 1:length(attr)) {
      m1=attr[[i]]
      name1=names(attr[i])
      if (name1==selCol) {
        selValues=m1
      }
    }
  }
  
  tagList(
    selectizeInput("sel_values_meta_col", label=menu_label, choices=UniqueValues, selected=selValues, multiple=TRUE, width="80%"),
    actionButton("clear_values_meta_col", "Deselect All Values"),
    tags$p("Number of samples for each value from all samples"),
    tableOutput("table_selCol"),
    tags$p("Number of samples for each value from selected samples"),
    tableOutput("table_selCol2"),
    HTML("<hr>")
  )
  }
})

observeEvent(input$clear_values_meta_col, {  
  #cat("reset values!\n")
  updateSelectizeInput(session,'sel_values_meta_col', selected="")
})

observeEvent(input$order_groups_order, {  
  group_order(input$order_groups_order)
})
observeEvent(input$order_samples_order, {  
  sample_order(input$order_samples_order)
})
observeEvent(input$QC_groups, {  
  group_order(input$QC_groups)
})

output$selectGroupSample <- renderText({ paste("Selected ",length(group_order()), " out of ", length(all_groups()), " Groups, ", 
                                               length(sample_order()), " out of ", length(all_samples()), " Samples.", sep="")})


#update sample list and sample order if group changes
observe({
  req(group_order())
  MetaData=all_metadata()
  groups = group_order()
  MetaData1<-MetaData%>%filter(group %in% groups)
  samplesG <- as.character( MetaData1$sampleid[order(match(MetaData1$group,groups))])
  sample_order(samplesG)
  attribute_filters(NULL)
  samples_excludeM(""); samples_excludeF("")
  resetComp2Sample(TRUE)
})

observe({
 if ( resetComp2Sample() ) {
   updateCheckboxInput(session, "QC_comp2sample", value = FALSE)
   updateSelectizeInput(session,'QC_comp_4_samples',  selected="None")
   resetComp2Sample(FALSE)
 }
})

#update sample list based on filter and manual list
observe({
  req(input$Select_Sample)
  if (input$Select_Sample=='Sample Filter') {
    #get samplese from group first
    MetaData=all_metadata()
    groups = group_order()
    MetaData1<-MetaData%>%filter(group %in% groups)
    sampleG <- as.character( MetaData1$sampleid[order(match(MetaData1$group,groups))])
    samples=sampleG
    #samples=sample_order()
    sample_R= unique(c(samples_excludeM(), samples_excludeF())) #extra samples to remove
    ToRemove=( toupper(samples) %in% toupper(sample_R) )
   # browser() #bebug
    if  (sum(ToRemove)>0) {samples=samples[!ToRemove]; 
      if (length(samples)>0) {
        sample_order(samples)
      } else {
        cat("No samples left after applying sample filter, reset!\n")
        sample_order(sampleG)
        attribute_filters(NULL)
        samples_excludeM(""); samples_excludeF("")
      }
     
    }
  }
})

#update sample list based on comparison (left menu)
observe({
  req(input$QC_comp_4_samples)
  comp1<-input$QC_comp_4_samples
  sel_comp=DataReactive()$sel_comp
  subsets<-sel_comp%>%filter(!duplicated(Subsetting_group), str_detect(Subsetting_group, ":"))
  if (nrow(subsets)>0){
    subsets<-subsets%>%mutate(Comparison=str_c("(Subset) ", Subsetting_group), sample_list=subset_list, Comparison=NA)
    sel_comp=rbind(sel_comp, subsets)
  }
  #browser() #bebug
  if (length(comp1)==1) {
  if (comp1!="None"){
    if (comp1=="All_Samples"){
      samples=all_samples()
    } else {
      samples=sel_comp$sample_list[sel_comp$Comparison==comp1]
      samples=str_split(samples, ",")[[1]]
    }
    sample_order(samples)
    attribute_filters(NULL)
    samples_excludeM(""); samples_excludeF("")
    updateRadioButtons(session, "Select_Sample", selected = "None")
  }
  } else {
    if ("All_Samples" %in% comp1) {
      samples=all_samples()
    } else {
      samples<-sel_comp%>%dplyr::filter(Comparison %in% comp1)%>%dplyr::select(sample_list)%>%unlist%>%unname%>%paste(collapse=",")
      samples=str_split(samples, ",")[[1]]
    }
    sample_order(samples)
    attribute_filters(NULL)
    samples_excludeM(""); samples_excludeF("")
    updateRadioButtons(session, "Select_Sample", selected = "None")
  }
})


#update sample list based on comparison (right panel)
observe({
  req(input$comp_4_samples)
  comp1<-input$comp_4_samples
  sel_comp=DataReactive()$sel_comp
  subsets<-sel_comp%>%filter(!duplicated(Subsetting_group), str_detect(Subsetting_group, ":"))
  if (nrow(subsets)>0){
    subsets<-subsets%>%mutate(Comparison=str_c("(Subset) ", Subsetting_group), sample_list=subset_list, Comparison=NA)
    sel_comp=rbind(sel_comp, subsets)
  }
  
  if (length(comp1)==1) {
    if (comp1!="None"){
      if (comp1=="All_Samples"){
        samples=all_samples()
      } else {
        samples=sel_comp$sample_list[sel_comp$Comparison==comp1]
        samples=str_split(samples, ",")[[1]]
      }
      sample_order(samples)
      attribute_filters(NULL)
      samples_excludeM(""); samples_excludeF("")
      resetComp2Sample(TRUE)
    }
  } else {
    if ("All_Samples" %in% comp1) {
      samples=all_samples()
    } else {
      samples<-sel_comp%>%dplyr::filter(Comparison %in% comp1)%>%dplyr::select(sample_list)%>%unlist%>%unname%>%paste(collapse=",")
      samples=str_split(samples, ",")[[1]]
    }
    sample_order(samples)
    attribute_filters(NULL)
    samples_excludeM(""); samples_excludeF("")
    resetComp2Sample(TRUE)
  }
})

observeEvent(input$QC_samples, {  
  sample_order(input$QC_samples)
})

observeEvent(input$remove_sample, {  
  sample_list=input$sample_exclude_list
  if(grepl("\n",sample_list)) {
    sample_list <-  stringr::str_split(sample_list, "\n")[[1]]
  } else if(grepl(",",sample_list)) {
    sample_list <-  stringr::str_split(sample_list, ",")[[1]]
  }
  sample_list <- gsub(" ", "", sample_list, fixed = TRUE)
  sample_list <- unique(sample_list[sample_list != ""])
  samples=all_samples()
  ToRemove=( toupper(samples) %in% toupper(sample_list) )
  #browser() #debug
  if  (sum(ToRemove)>0) {
    samples_excludeM(samples[ToRemove])
  }
})

observeEvent(input$upload_sample, {  
  sample_list=input$sample_upload_list
  if(grepl("\n",sample_list)) {
    sample_list <-  stringr::str_split(sample_list, "\n")[[1]]
  }  else if(grepl(",",sample_list)) {
    sample_list <-  stringr::str_split(sample_list, ",")[[1]]
  }
  sample_list <- gsub(" ", "", sample_list, fixed = TRUE)
  sample_list <- unique(sample_list[sample_list != ""])
  samples=all_samples()
  ToAdd=( toupper(sample_list) %in% toupper(samples) )
  #browser() #debug
  if  (sum(ToAdd)>0) {
    sample_order(sample_list[ToAdd])
    attribute_filters(NULL)
    samples_excludeM(""); samples_excludeF("")
    resetComp2Sample(TRUE)
  }
})

#remove samples from filter on attributes
observeEvent(input$sel_values_meta_col, {  
  MetaData=all_metadata()
  selCol=input$meta_col_sel
  values=as.character(MetaData[[selCol]])
  samples=MetaData$sampleid
  UniqueValues=sort(unique(values))
  selValues=input$sel_values_meta_col
  selValues1=selValues
  selValues[selValues=="Empty_Value"]=""
  RemoveValues=setdiff(UniqueValues, selValues)
  #check if selCol in filters already
  attr=attribute_filters()
  InFilter=0
  if (!is.null(attr)) {
    for (i in 1:length(attr)) {
      name1=names(attr[i])
      if (name1==selCol ) {
        InFilter=i
      }
    }
  }
  
  if (length(RemoveValues)>0 || InFilter>0)  {
    Fnew=list(selValues1)
    names(Fnew)=selCol
    #add selection to attribute_filters
    if (is.null(attr) & length(RemoveValues)>0) {
      attr=Fnew
    } else if (InFilter>0) {
      attr[InFilter]=Fnew
    } else if  (!is.null(attr) & length(RemoveValues)>0) {
      attr=c(attr, Fnew)
    }
  
    attribute_filters(attr)
    resetComp2Sample(TRUE)
    # browser() #debug
    #now loop through attributes
    ToRemove=rep(F, length(samples))
    for (i in 1:length(attr)) {
      selValues=attr[[i]]
      selCol=names(attr[i])
      values=as.character(MetaData[[selCol]])
      UniqueValues=sort(unique(values))
      selValues[selValues=="Empty_Value"]=""
      RemoveValues=setdiff(UniqueValues, selValues)
      NRemove=which(values %in% RemoveValues)
      ToRemove[NRemove]=T
    }

    if  (sum(ToRemove)>0) {
      samples_excludeF(samples[ToRemove])
    }
  }
})




observeEvent(input$reset_group, {
  allgroups = all_groups()
  group_order(allgroups)
  samples_excludeM("")
  attribute_filters(NULL)
  samples_excludeF("")
  samples=all_samples()
  sample_order(samples)
})

DataQCReactive <- reactive({
	DataIn = DataReactive()
	results_long = DataIn$results_long
	ProteinGeneName = DataIn$ProteinGeneName
	MetaData = DataIn$MetaData
	data_long = DataIn$data_long
	data_wide = DataIn$data_wide

	input_groups = input$QC_groups
	#group_order(input$QC_groups)
	input_samples = input$QC_samples
	tmp_data_long = dplyr::filter(data_long, (group %in% input_groups) & (sampleid %in% input_samples))
	
	sel_sample_order=match(input_samples, MetaData$sampleid)
	sel_sample_order=sel_sample_order[!is.na(sel_sample_order)]
	MetaData=MetaData[sel_sample_order, ] #user input sample order
	input_keep = which(MetaData$group %in% input_groups) 
	MetaData=MetaData[input_keep, ]
	tmp_group = MetaData$group
	tmp_sampleid = MetaData$sampleid
	data_wide  <- data_wide[apply(data_wide, 1, function(x) sum(length(which(x==0 | is.na(x)))) < 3),]
	sel_sample_order2=match(tmp_sampleid, colnames(data_wide))
	sel_sample_order2=sel_sample_order2[!is.na(sel_sample_order2)]
	tmp_data_wide = data_wide[, sel_sample_order2] %>% as.matrix()
	
	return(list('tmp_data_wide'=tmp_data_wide,'tmp_data_long'=tmp_data_long,'tmp_group' = tmp_group, 'tmp_sampleid'=tmp_sampleid, "MetaData"=MetaData ))
})

DataPCAReactive <- reactive({
	DataQC <-  DataQCReactive()
	tmp_sampleid <- DataQC$tmp_sampleid
	validate(need(length(tmp_sampleid)>1, message = "Please select at least two samples (please note samples are filtered by group selection as well)."))
	tmp_data_wide <- DataQC$tmp_data_wide
	tmp_group = DataQC$tmp_group

	tmp_data_wide[is.na(tmp_data_wide)] <- 0 
	pca <- 	prcomp(t(tmp_data_wide),rank. = 10, scale = FALSE)
	percentVar <- 	round((pca$sdev)^2/sum(pca$sdev^2), 3) * 100
	scores <- as.data.frame(pca$x)
	rownames(scores) <- tmp_sampleid
	scores$group <- factor(tmp_group, levels = group_order())
	attributes=setdiff(colnames(DataQC$MetaData), c("Order", "ComparePairs", "group") )
	MetaData=DataQC$MetaData
	colsel=match(attributes, colnames(MetaData) )
	scores=cbind(scores, MetaData[, colsel, drop=F])
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
		tmp_data_long <- DataQC$tmp_data_long %>% dplyr::filter(expr !=0) %>% sample_n(1000)
			
		tmp_group = DataQC$tmp_group
		#colorpal = get_palette("Dark2", length(tmp_sampleid))
	
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
  ptm <- proc.time()
  req(DataPCAReactive())
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
	#colorpal = topo.colors(n, alpha = 1)
	#colorpal = get_palette("Dark2", n)
	colorpal = colorRampPalette(brewer.pal(8, input$PCAcolpalette))(n)
	
	#if (all(table(tmp_group))<4)
	#  ellipsoid = FALSE 

	if (input$PCA_subsample=="None" ) {labels=NULL
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

observe({
  
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
  colorpal = colorRampPalette(brewer.pal(8, input$PCAcolpalette))(n)
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
	colorpal = colorRampPalette(brewer.pal(8, input$PCAcolpalette))(n)
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
	colorpal = colorRampPalette(brewer.pal(8, input$PCAcolpalette))(n)
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
	tmp_group = DataQC$tmp_group
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
	tmp_group = DataQC$tmp_group
	if (!"id" %in% names(tmp_data_long)) {tmp_data_long$id=tmp_data_long$UniqueID} 
	#browser()
	CV.df <- tmp_data_long %>%
	group_by(.,  group, id) %>%
	dplyr::summarise( mean=mean(expr, na.rm = TRUE), sd=sd(expr, na.rm = TRUE)) %>%
	dplyr::mutate(CV=100*(sd/mean))

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
  res<-Covariate_PC_Analysis(tmp_data_wide, meta, out_prefix=NULL, PC_cutoff=input$covar_PC_cutoff, 
            FDR_cutoff=input$covar_FDR_cutoff, N_col=input$covar_ncol)
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
