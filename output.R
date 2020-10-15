###########################################################################################################
## Proteomics Visualization R Shiny App
##
##This software belongs to Biogen Inc. All right reserved.
##
##@file: ouput.R
##@Developer : Benbo Gao (benbo.gao@Biogen.com)
##@Date : 5/16/2018
##@version 1.0
###########################################################################################################


output$downloadPDF <- downloadHandler(
	filename = function() {
		paste("output_",Sys.Date(),".pdf", sep="")
	},
	content = function(file) { withProgress(message = 'This takes a minute or two',{ 
		pdf(file = file, width=input$pdf_width, height=input$pdf_height)
		Title <- paste("","Data Analysis Report",sep=" ")
		today <- Sys.Date() %>%	format(format="%m/%d/%Y")

		
		### Boxplot
		if (!is.null(saved_plots$QCboxplot)){
		print(saved_plots$QCboxplot)
		}

		## PCA
		if (!is.null(saved_plots$pcaplot)){
			print(saved_plots$pcaplot)
		}
		#sample-to-sample distances
		if (!is.null(saved_plots$SampleDistance)){
		  plot.new()
			grid.draw(saved_plots$SampleDistance)
		}
		#Dendrograms
		if (!is.null(saved_plots$Dendrograms)){
			print(saved_plots$Dendrograms)
		}
		###CV distribution
		if (!is.null(saved_plots$histplot)){
			print(saved_plots$histplot)
		}
		### volcano plot
		if (!is.null(saved_plots$volcano)){
			for (i in 1:length(names(saved_plots$volcano))) {
				p.name <- names(saved_plots$volcano)[i]
				title_tmp <- paste("Volcano Plot (",p.name,")",sep="")
			print(saved_plots$volcano[[p.name]])
			}
		}
		## Heatmap
		if (!is.null(saved_plots$pheatmap2)){
			draw(saved_plots$pheatmap2, merge_legend=T,auto_adjust = FALSE)
		}
		## Heatmap
		if (!is.null(saved_plots$staticheatmap)){
			replayPlot(saved_plots$staticheatmap)
		}
		### boxplot plot
		if (!is.null(saved_plots$boxplot)){
			for (i in 1:length(saved_plots$boxplot)) {
				title_tmp <- paste("boxplot (",i,")",sep="")
				print(saved_plots$boxplot[[i]])
			}
		}
		### browing plot
		if (!is.null(saved_plots$browsing)){
			for (i in 1:length(saved_plots$browsing)) {
				title_tmp <- paste("Browsing Plot (",i,")",sep="")
				print(saved_plots$browsing[[i]])
			}
		}
		###vennDiagram
		if (!is.null(saved_plots$vennDiagram)){
			for (i in 1:length(saved_plots$vennDiagram)) {
				title_tmp <- paste("vennDiagram (",i,")",sep="")
				plot.new()
				grid.draw(saved_plots$vennDiagram[[i]])
			}
		}
		## pattern
		if (!is.null(saved_plots$patternkmeans)){
			print(saved_plots$patternkmeans)
		}
		if (!is.null(saved_plots$patternpam)){
			print(saved_plots$patternpam)
		}
		if (!is.null(saved_plots$patternmfuzz)){
			replayPlot(saved_plots$patternmfuzz)
		}

		dev.off()
	})
	},contentType = "application/pdf"

)

observeEvent(input$clear_saved_plots, {  
  #cat("clear saved plots\n")
  lapply(X = names(saved_plots), FUN = function(x) { saved_plots[[x]] <- NULL }) 
})

#list all saved plots
saved_plot_list<-reactive({
  req(saved_plots)
 
  summary=str_c('<style type="text/css">
.disc {
 list-style-type: disc;
}
</style>',
"<h4>Saved Plots for Project: ", ProjectInfo$ShortName, "</h4><br>", '<ul class="disc">')
  
  ### Boxplot
  if (!is.null(saved_plots$QCboxplot)){
   summary=str_c(summary, '<li>QC Boxplot</li>')
  }
  ## PCA
  if (!is.null(saved_plots$pcaplot)){
    summary=str_c(summary, '<li>PCA Plot</li>')
  }
  #sample-to-sample distances
  if (!is.null(saved_plots$SampleDistance)){
    summary=str_c(summary, '<li>Sample Distance Plot</li>')
  }
  #Dendrograms
  if (!is.null(saved_plots$Dendrograms)){
    summary=str_c(summary, '<li>Sample Dendrogram</li>')
  }
  ###CV distribution
  if (!is.null(saved_plots$histplot)){
    summary=str_c(summary, '<li>CV Distribution</li>')
  }
  ### volcano plot
  if (!is.null(saved_plots$volcano)){
    for (i in 1:length(names(saved_plots$volcano))) {
      p.name <- names(saved_plots$volcano)[i]
      title_tmp <- paste("Volcano Plot (",p.name,")",sep="")
      summary=str_c(summary, '<li>', title_tmp, '</li>')
    }
  }
  ## Heatmap
  if (!is.null(saved_plots$pheatmap2)){
    summary=str_c(summary, '<li>Heatmap</li>')
  }
  ## Heatmap
  if (!is.null(saved_plots$staticheatmap)){
    summary=str_c(summary, '<li>Heatmap layout2</li>')
  }
  ### boxplot plot
  if (!is.null(saved_plots$boxplot)){
    for (i in 1:length(saved_plots$boxplot)) {
      title_tmp <- paste("boxplot (",i,")",sep="")
      summary=str_c(summary, '<li>', title_tmp, '</li>')
    }
  }
  ### browing plot
  if (!is.null(saved_plots$browsing)){
    for (i in 1:length(saved_plots$browsing)) {
      title_tmp <- paste("Browsing Plot (",i,")",sep="")
      summary=str_c(summary, '<li>', title_tmp, '</li>')
    }
  }
  ###vennDiagram
  if (!is.null(saved_plots$vennDiagram)){
    for (i in 1:length(saved_plots$vennDiagram)) {
      title_tmp <- paste("vennDiagram (",i,")",sep="")
      summary=str_c(summary, '<li>', title_tmp, '</li>')
    }
  }
  ## pattern
  if (!is.null(saved_plots$patternkmeans)){
    summary=str_c(summary, '<li>Pattern Clustering (K-means)</li>')
  }
  if (!is.null(saved_plots$patternpam)){
    summary=str_c(summary, '<li>Pattern Clustering (PAM)</li>')
    
  }
  if (!is.null(saved_plots$patternmfuzz)){
    summary=str_c(summary, '<li>Pattern Clustering (Soft Clustering)</li>')
  }
summary=str_c(summary, "</ul><br><hr>")
})
output$saved_plot_list=renderText(saved_plot_list())



##Product SVG for first Graph
output$downloadSVG <- downloadHandler(
  filename = function() {
    paste("output_",Sys.Date(),".svg", sep="")
  },
  content = function(file) { withProgress(message = 'This takes a minute or two',{ 
    svglite(file = file, width=input$pdf_width, height=input$pdf_height)
    Title <- paste("","Data Analysis Report",sep=" ")
    today <- Sys.Date() %>%	format(format="%m/%d/%Y")
    
    ### Boxplot
    if (!is.null(saved_plots$QCboxplot)){
      print(saved_plots$QCboxplot)
    } else   if (!is.null(saved_plots$pcaplot)){
      print(saved_plots$pcaplot)
    } else if (!is.null(saved_plots$SampleDistance)){
     # plot.new()
      grid.draw(saved_plots$SampleDistance)
    } else  if (!is.null(saved_plots$Dendrograms)){
      print(saved_plots$Dendrograms)
    } else if (!is.null(saved_plots$histplot)){
      print(saved_plots$histplot)
    } else if (!is.null(saved_plots$volcano)){
        i=1
        p.name <- names(saved_plots$volcano)[i]
        title_tmp <- paste("Volcano Plot (",p.name,")",sep="")
        print(saved_plots$volcano[[p.name]])
    } else if (!is.null(saved_plots$pheatmap2)){
      draw(saved_plots$pheatmap2, merge_legend=T,auto_adjust = FALSE)
    } else if (!is.null(saved_plots$staticheatmap)){
      replayPlot(saved_plots$staticheatmap)
    } else if (!is.null(saved_plots$boxplot)){
        i=1
        title_tmp <- paste("boxplot (",i,")",sep="")
        print(saved_plots$boxplot[[i]])
    } else if (!is.null(saved_plots$browsing)){
        i=1
        title_tmp <- paste("Browsing Plot (",i,")",sep="")
        print(saved_plots$browsing[[i]])
    } else if (!is.null(saved_plots$vennDiagram)){
        i=1
        title_tmp <- paste("vennDiagram (",i,")",sep="")
        #plot.new()
        grid.draw(saved_plots$vennDiagram[[i]])
    } else if (!is.null(saved_plots$patternkmeans)){
      print(saved_plots$patternkmeans)
    }  else if (!is.null(saved_plots$patternpam)){
      print(saved_plots$patternpam)
    } else if (!is.null(saved_plots$patternmfuzz)){
      replayPlot(saved_plots$patternmfuzz)
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

