###########################################################################################################
## Proteomics Visualization R Shiny App
##
##This software belongs to Biogen Inc. All right reserved.
##
##@file: ouput.R
##@Developer : Benbo Gao (benbo.gao@Biogen.com)
##@Date : 3/18/2021
##@version 2.0
###########################################################################################################
saved_plot_list<- reactiveVal()

output$downloadPDF <- downloadHandler(
	filename = function() {
		paste("output_",Sys.Date(),".pdf", sep="")
	},
	content = function(file) { withProgress(message = 'This takes a minute or two',{
		plots_checked=input$plots_checked
		validate(need(length(plots_checked)>0, message = "Please choose at least one saved plot."))
		Np=0
		pdf(file = file, width=input$pdf_width, height=input$pdf_height)

		### Boxplot
		if (!is.null(saved_plots$QCboxplot) & ("QC Boxplot" %in% plots_checked)){
			print(saved_plots$QCboxplot)
			Np=Np+1
		}

		## PCA
		if (!is.null(saved_plots$pcaplot) & ("PCA Plot" %in% plots_checked)){
			print(saved_plots$pcaplot)
			Np=Np+1
		}
		
		#sample-to-sample distances
		if (!is.null(saved_plots$SampleDistance) & ("Sample Distance Plot" %in% plots_checked) ){
			grid.newpage()
			pushViewport(viewport(width=0.8,height=0.8))
			grid.draw(saved_plots$SampleDistance)
			Np=Np+1
		}
		
		#Dendrograms
		if (!is.null(saved_plots$Dendrograms) & ("Sample Dendrogram" %in% plots_checked)){
			print(saved_plots$Dendrograms)
			Np=Np+1
		}
		
		###CV distribution
		if (!is.null(saved_plots$histplot) & ("CV Distribution" %in% plots_checked)){
			print(saved_plots$histplot)
			Np=Np+1
		}
		
		
		#Covariates
		if (!is.null(saved_plots$covar_cat) & ("Categorical Covariates vs PCs" %in% plots_checked)){
		  print(saved_plots$covar_cat)
		  Np=Np+1
		}
		
		if (!is.null(saved_plots$covar_num) & ("Numeric Covariates vs PCs" %in% plots_checked)){
		  print(saved_plots$covar_num)
		  Np=Np+1
		}
		
		
		### volcano plot
		if (!is.null(saved_plots$volcano)){
			for (i in 1:length(names(saved_plots$volcano))) {
				p.name <- names(saved_plots$volcano)[i]
				title_tmp <- paste("Volcano Plot (",p.name,")",sep="")
				if (title_tmp %in% plots_checked) {
					print(saved_plots$volcano[[p.name]])
					Np=Np+1
				}
			}
		}
		
		## Heatmap
		if (!is.null(saved_plots$pheatmap2) & ("Heatmap" %in% plots_checked)){
			draw(saved_plots$pheatmap2, merge_legend=T,auto_adjust = FALSE)
			Np=Np+1
		}
		
		## Heatmap
		if (!is.null(saved_plots$staticheatmap)  & ("Heatmap layout2" %in% plots_checked)){
			replayPlot(saved_plots$staticheatmap)
			Np=Np+1
		}
		
		## Geneset Heatmap
		if (!is.null(saved_plots$genesetheatmap) & ("Geneset Heatmap" %in% plots_checked)){
	    draw(saved_plots$genesetheatmap, merge_legend=T,auto_adjust = FALSE)
		  Np=Np+1
		}

		## Geneset KEGG View
		if (!is.null(saved_plots$keggSave)){
			for (i in 1:length(names(saved_plots$keggSave))) {
				p.name <- names(saved_plots$keggSave)[i]
				title_tmp <- paste("KEGG View (",p.name,")",sep="")
				if (title_tmp %in% plots_checked) {
					par(mai=c(0,0,0,0))
					plot(c(0,1),c(0,1),type="n")
					rasterImage(saved_plots$keggSave[[i]],0,0,1,1)
					Np=Np+1
				}
			}
		}

		### boxplot plot
		if (!is.null(saved_plots$boxplot)){
			for (i in 1:length(saved_plots$boxplot)) {
				title_tmp <- paste("boxplot (",i,")",sep="")
				if (title_tmp %in% plots_checked) {
					print(saved_plots$boxplot[[i]])
					Np=Np+1
				}
			}
		}
		
		### browing plot
		if (!is.null(saved_plots$browsing)){
			for (i in 1:length(saved_plots$browsing)) {
				title_tmp <- paste("Browsing Plot (",i,")",sep="")
				if (title_tmp %in% plots_checked) {
					print(saved_plots$browsing[[i]])
					Np=Np+1
				}
			}
		}
		
		###vennDiagram
		if (!is.null(saved_plots$vennDiagram)){
			for (i in 1:length(saved_plots$vennDiagram)) {
				title_tmp <- paste("vennDiagram (",i,")",sep="")
				if (title_tmp %in% plots_checked) {
					grid.newpage()
					pushViewport(viewport(width=0.8,height=0.8))
					grid.draw(saved_plots$vennDiagram[[i]])
					Np=Np+1
				}
			}
		}
		
		## pattern
		if (!is.null(saved_plots$patternkmeans) & ("Pattern Clustering (K-means)" %in% plots_checked)){
			print(saved_plots$patternkmeans)
			Np=Np+1
		}
		if (!is.null(saved_plots$patternpam) & ("Pattern Clustering (PAM)" %in% plots_checked)){
			print(saved_plots$patternpam)
			Np=Np+1
		}
		if (!is.null(saved_plots$patternmfuzz) & ("Pattern Clustering (Soft Clustering)" %in% plots_checked)){
			replayPlot(saved_plots$patternmfuzz)
			Np=Np+1
		}
		dev.off()
		cat("Saved to PDF", Np, "graphs.\n")
	})
},contentType = "application/pdf"

)

observeEvent(input$clear_saved_plots, {
	cat("clear saved plots\n")
	lapply(X = names(saved_plots), FUN = function(x) { saved_plots[[x]] <- NULL })
})

#list all saved plots
observe({
	req(saved_plots)
	summary=NULL

	if (!is.null(saved_plots$QCboxplot)){
		summary=c(summary, 'QC Boxplot')
	}
	
	## PCA
	if (!is.null(saved_plots$pcaplot)){
		summary=c(summary, 'PCA Plot')
	}
	
	#sample-to-sample distances
	if (!is.null(saved_plots$SampleDistance)){
		summary=c(summary, 'Sample Distance Plot')
	}
	
	#Dendrograms
	if (!is.null(saved_plots$Dendrograms)){
		summary=c(summary, 'Sample Dendrogram')
	}
	
	###CV distribution
	if (!is.null(saved_plots$histplot)){
		summary=c(summary, 'CV Distribution')
	}
	
	###Covariates
	if (!is.null(saved_plots$covar_cat)){
	  summary=c(summary, 'Categorical Covariates vs PCs')
	}
	if (!is.null(saved_plots$covar_num)){
	  summary=c(summary, 'Numeric Covariates vs PCs')
	}
	
	### volcano plot
	if (!is.null(saved_plots$volcano)){
		for (i in 1:length(names(saved_plots$volcano))) {
			p.name <- names(saved_plots$volcano)[i]
			title_tmp <- paste("Volcano Plot (",p.name,")",sep="")
			summary=c(summary,  title_tmp)
		}
	}
	
	## Heatmap
	if (!is.null(saved_plots$pheatmap2)){
		summary=c(summary, 'Heatmap')
	}
	
	## Heatmap
	if (!is.null(saved_plots$staticheatmap)){
		summary=c(summary, 'Heatmap layout2')
	}
	
	## Geneset Heatmap
	if (!is.null(saved_plots$genesetheatmap)){
	  summary=c(summary, 'Geneset Heatmap')
	}
	
	## Geneset KeggView
	if (!is.null(saved_plots$keggSave)){
		for (i in 1:length(names(saved_plots$keggSave))) {
			p.name <- names(saved_plots$keggSave)[i]
			title_tmp <- paste("KEGG View (",p.name,")",sep="")
			summary=c(summary,  title_tmp)
		}
	}


	### boxplot plot
	if (!is.null(saved_plots$boxplot)){
		for (i in 1:length(saved_plots$boxplot)) {
			title_tmp <- paste("boxplot (",i,")",sep="")
			summary=c(summary, title_tmp)
		}
	}
	
	### browing plot
	if (!is.null(saved_plots$browsing)){
		for (i in 1:length(saved_plots$browsing)) {
			title_tmp <- paste("Browsing Plot (",i,")",sep="")
			summary=c(summary, title_tmp)
		}
	}
	
	###vennDiagram
	if (!is.null(saved_plots$vennDiagram)){
		for (i in 1:length(saved_plots$vennDiagram)) {
			title_tmp <- paste("vennDiagram (",i,")",sep="")
			summary=c(summary, title_tmp)
		}
	}
	
	## pattern
	if (!is.null(saved_plots$patternkmeans)){
		summary=c(summary, 'Pattern Clustering (K-means)')
	}
	if (!is.null(saved_plots$patternpam)){
		summary=c(summary, 'Pattern Clustering (PAM)')

	}
	if (!is.null(saved_plots$patternmfuzz)){
		summary=c(summary, 'Pattern Clustering (Soft Clustering)')
	}
	
	#cat("saved plots are:", summary, "\n")
	saved_plot_list(summary)
	if (is.null(summary)) {saved_plot_list(character(0))}
})

observe({
	updateCheckboxGroupInput(session, "plots_checked", choices=saved_plot_list(), selected=saved_plot_list() )
})


##Product SVG for first selected Graph
output$downloadSVG <- downloadHandler(
  filename = function() {
    paste("output_",Sys.Date(),".svg", sep="")
  },
  content = function(file) { withProgress(message = 'This takes a minute or two',{ 
    plots_checked=input$plots_checked
    validate(need(length(plots_checked)>0, message = "Please choose at least one saved plot."))
    Np=0
    svglite(file = file, width=input$pdf_width, height=input$pdf_height)
    Title <- paste("","Data Analysis Report",sep=" ")
    today <- Sys.Date() %>%	format(format="%m/%d/%Y")
    ### Boxplot
    if (!is.null(saved_plots$QCboxplot) & ("QC Boxplot" %in% plots_checked)){
      print(saved_plots$QCboxplot)
      Np=Np+1
    }
    
    ## PCA
    if (Np==0 & !is.null(saved_plots$pcaplot) & ("PCA Plot" %in% plots_checked)){
      print(saved_plots$pcaplot)
      Np=Np+1
    }
    #sample-to-sample distances
    if (Np==0 & !is.null(saved_plots$SampleDistance) & ("Sample Distance Plot" %in% plots_checked) ){
      plot.new()
      grid.draw(saved_plots$SampleDistance)
      Np=Np+1
    }
    #Dendrograms
    if (Np==0 & !is.null(saved_plots$Dendrograms) & ("Sample Dendrogram" %in% plots_checked)){
      print(saved_plots$Dendrograms)
      Np=Np+1
    }
    ###CV distribution
    if (Np==0 & !is.null(saved_plots$histplot) & ("CV Distribution" %in% plots_checked)){
      print(saved_plots$histplot)
      Np=Np+1
    }
    ###Covariates
    if (!is.null(saved_plots$covar_cat) & ("Categorical Covariates vs PCs" %in% plots_checked)){
      print(saved_plots$covar_cat)
      Np=Np+1
    }
    if (!is.null(saved_plots$covar_num) & ("Numeric Covariates vs PCs" %in% plots_checked)){
      print(saved_plots$covar_num)
      Np=Np+1
    }
    
    
    ### volcano plot
    if (Np==0 & !is.null(saved_plots$volcano)){
      for (i in 1:length(names(saved_plots$volcano))) {
        p.name <- names(saved_plots$volcano)[i]
        title_tmp <- paste("Volcano Plot (",p.name,")",sep="")
        if (Np==0 & title_tmp %in% plots_checked) {
          print(saved_plots$volcano[[p.name]])
          Np=Np+1
        }
      }
    }
    ## Heatmap
    if (Np==0 & !is.null(saved_plots$pheatmap2) & ("Heatmap" %in% plots_checked)){
      draw(saved_plots$pheatmap2, merge_legend=T,auto_adjust = FALSE)
      Np=Np+1
    }
    ## Heatmap
    if (Np==0 & !is.null(saved_plots$staticheatmap)  & ("Heatmap layout2" %in% plots_checked)){
      replayPlot(saved_plots$staticheatmap)
      Np=Np+1
    }
    
    
    ## Geneset Heatmap
    if (Np==0 & !is.null(saved_plots$genesetheatmap) & ("Geneset Heatmap" %in% plots_checked)){
      draw(saved_plots$genesetheatmap, merge_legend=T,auto_adjust = FALSE)
      Np=Np+1
    }
    
    ## Geneset KEGG View
    if (Np==0 &!is.null(saved_plots$keggSave)){
      for (i in 1:length(names(saved_plots$keggSave))) {
        p.name <- names(saved_plots$keggSave)[i]
        title_tmp <- paste("KEGG View (",p.name,")",sep="")
        if (title_tmp %in% plots_checked) {
          par(mai=c(0,0,0,0))
          plot(c(0,1),c(0,1),type="n")
          rasterImage(saved_plots$keggSave[[i]],0,0,1,1)
          Np=Np+1
        }
      }
    }
    
    ### boxplot plot
    if (Np==0 & !is.null(saved_plots$boxplot)){
      for (i in 1:length(saved_plots$boxplot)) {
        title_tmp <- paste("boxplot (",i,")",sep="")
        if (Np==0 & title_tmp %in% plots_checked) {
          print(saved_plots$boxplot[[i]])
          Np=Np+1
        }
      }
    }
    ### browing plot
    if (Np==0 & !is.null(saved_plots$browsing)){
      for (i in 1:length(saved_plots$browsing)) {
        title_tmp <- paste("Browsing Plot (",i,")",sep="")
        if (Np==0 & title_tmp %in% plots_checked) {
          print(saved_plots$browsing[[i]])
          Np=Np+1
        }
      }
    }
    ###vennDiagram
    if (Np==0 & !is.null(saved_plots$vennDiagram)){
      for (i in 1:length(saved_plots$vennDiagram)) {
        title_tmp <- paste("vennDiagram (",i,")",sep="")
        if (Np==0 & title_tmp %in% plots_checked) {
          plot.new()
          grid.draw(saved_plots$vennDiagram[[i]])
          Np=Np+1
        }
      }
    }
    ## pattern
    if (Np==0 & !is.null(saved_plots$patternkmeans) & ("Pattern Clustering (K-means)" %in% plots_checked)){
      print(saved_plots$patternkmeans)
      Np=Np+1
    }
    if (Np==0 & !is.null(saved_plots$patternpam) & ("Pattern Clustering (PAM)" %in% plots_checked)){
      print(saved_plots$patternpam)
      Np=Np+1
    }
    if (Np==0 & !is.null(saved_plots$patternmfuzz) & ("Pattern Clustering (Soft Clustering)" %in% plots_checked)){
      replayPlot(saved_plots$patternmfuzz)
      Np=Np+1
    }
    dev.off()
    
 
  })
  },contentType = "application/svg"
  
)


output$downloadXLSX <- downloadHandler( 

	filename = function() {
		paste("output_",Sys.Date(),".xlsx", sep="")
	},
	content = function(file) { withProgress(message = 'This takes a minute or two',{
		wb <- createWorkbook(file)
		for (i in names(saved_table)) {
			addWorksheet(wb, i)
			res_table <- saved_table[[i]]
			writeData(wb, i, res_table,rowNames = TRUE)
		}
		saveWorkbook(wb, file = file, overwrite = TRUE)
	})
},contentType = "application/vnd.ms-excel"
)

