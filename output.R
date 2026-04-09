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
saved_gct_list <- reactiveVal()

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

		## Geneset dot plot
		if (!is.null(saved_plots$dotplot)){
		  for (i in 1:length(names(saved_plots$dotplot))) {
		    title_tmp <- paste("KEGG View (", i,")",sep="")
		    if (title_tmp %in% plots_checked) {
		      par(mai=c(0,0,0,0))
		      plot(c(0,1),c(0,1),type="n")
		      rasterImage(saved_plots$dotplot[[i]],0,0,1,1)
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
		
		### Time Series Cluster
		if (!is.null(saved_plots$ts_cluster_plot)){
		  for (i in 1:length(saved_plots$ts_cluster_plot)) {
		    title_tmp <- paste("Time Series Cluster Plot (",i,")",sep="")
		    if (title_tmp %in% plots_checked) {
		      print(saved_plots$ts_cluster_plot[[i]])
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
		
		## correlation analysis
		if (!is.null(saved_plots$CorrPlot) ){
		  for (i in 1:length(names(saved_plots$CorrPlot))) {
		    p.name <- names(saved_plots$CorrPlot)[i]
		    title_tmp <- paste("Correlation Plot (",p.name,")",sep="")
		    if (title_tmp %in% plots_checked) {
		      print(saved_plots$CorrPlot[[p.name]])
		      Np=Np+1
		    }
		  }
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

	## Geneset dotplot
	if (!is.null(saved_plots$dotplot)){
	  for (i in 1:length(names(saved_plots$dotplot))) {
	    title_tmp <- paste("pathway dot plot (", i,")",sep="")
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

	## correlation analysis
	if (!is.null(saved_plots$CorrPlot)){
	  for (i in 1:length(names(saved_plots$CorrPlot))) {
	    p.name <- names(saved_plots$CorrPlot)[i]
	    title_tmp <- paste("Correlation Plot (",p.name,")",sep="")
	    summary=c(summary,  title_tmp)
	  }
	}
	
	## Time Series Cluster
	if (!is.null(saved_plots$ts_cluster_plot)){
	  for (i in 1:length(saved_plots$ts_cluster_plot)) {
	    title_tmp <- paste("Time Series Cluster Plot (",i,")",sep="")
	    summary=c(summary, title_tmp)
	  }
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
    ## correlation analysis
    if (Np==0 & !is.null(saved_plots$CorrPlot) ) {
      for (i in 1:length(names(saved_plots$CorrPlot))) {
        p.name <- names(saved_plots$CorrPlot)[i]
        title_tmp <- paste("Correlation Plot (",p.name,")",sep="")
        if (Np==0 & title_tmp %in% plots_checked) {
          print(saved_plots$CorrPlot[[p.name]])
          Np=Np+1
        }
      }
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

observe({
  req(saved_gcts)
  summary=NULL
  
  ## Heatmap
  if (!is.null(saved_gcts$heatmap_gct)){
    summary=c(summary, 'Heatmap GCT file')
  }
  
  ## WGCNA Cluster Heatmap
  if (!is.null(saved_gcts$wgcna_gct)){
    summary=c(summary, 'WGCNA cluster GCT file')
  }
  
  #Time Series Cluster Heatmap
  if (!is.null(saved_gcts$ts_tpm_gct)){
    summary=c(summary, 'Time Series cluster sample normalized log2TPM GCT file')
  }
  
  if (!is.null(saved_gcts$ts_zscore_gct)){
    summary=c(summary, 'Time Series cluster group z-score GCT file')
  }

  #cat("saved plots are:", summary, "\n")
  saved_gct_list(summary)
  if (is.null(summary)) {saved_gct_list(character(0))}
})

observe({
  updateCheckboxGroupInput(session, "GCT_table_checked", choices=saved_gct_list(), selected=saved_gct_list() )
})


output$downloadGCT <- downloadHandler(
  filename = function() {
    paste("GCT_Export_", Sys.Date(), ".zip", sep = "")
  },
  content = function(file) {
    # Create a temporary directory to hold the individual GCT files
    cwd <- getwd()
    tmpdir <- tempdir()
    setwd(tempdir())
    files_to_zip <- c()

    GCT_table_checked=input$GCT_table_checked
    
    withProgress(message = 'Preparing GCT files...', value = 0, {
      # Example: Check for Heatmap GCT
      if (!is.null(saved_gcts$heatmap_gct) && ("Heatmap GCT file" %in% GCT_table_checked)) {
        fname <- "heatmap_data.gct"
        cmapR::write_gct(saved_gcts$heatmap_gct, fname, appenddim = FALSE)
        files_to_zip <- c(files_to_zip, fname)
      }
      
      if (!is.null(saved_gcts$wgcna_gct) && ("WGCNA cluster GCT file" %in% GCT_table_checked)) {
        fname <- "wgcna_cluster_heatmap_data.gct"
        cmapR::write_gct(saved_gcts$wgcna_gct, fname, appenddim = FALSE)
        files_to_zip <- c(files_to_zip, fname)
      }
      
      if (!is.null(saved_gcts$ts_tpm_gct) && ("Time Series cluster sample normalized log2TPM GCT file" %in% GCT_table_checked)) {
        fname <- "ts_cluster_sample_normalized_log2TPM_data.gct"
        cmapR::write_gct(saved_gcts$ts_tpm_gct, fname, appenddim = FALSE)
        files_to_zip <- c(files_to_zip, fname)
      }
      if (!is.null(saved_gcts$ts_zscore_gct) && ("Time Series cluster group z-score GCT file" %in% GCT_table_checked)) {
        fname <- "ts_cluster_group_zscore_data.gct"
        cmapR::write_gct(saved_gcts$ts_zscore_gct, fname, appenddim = FALSE)
        files_to_zip <- c(files_to_zip, fname)
      }
      # Verify if any files were actually created
      if (length(files_to_zip) == 0) {
        writeLines("No GCT files were selected for download.", "README.txt")
        files_to_zip <- "README.txt"
      }
      # Zip the files together
      zip::zip(zipfile = file, files = files_to_zip)
    })
    setwd(cwd)
  },
  contentType = "application/zip"
)
