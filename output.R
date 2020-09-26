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
		pdf(file = file, width=12, height=8)
		Title <- paste("","Data Analysis Report",sep=" ")
		today <- Sys.Date() %>%	format(format="%m/%d/%Y")

		
		### Boxplot
		if (!is.null(saved_plots$QCboxplot)){
		print(saved_plots$QCboxplot)
		}

		## PCA
		if (!is.null(saved_plots$pcaplot)){
			print(saved_plots$pcaplot)
			plot.new()
		}
		#sample-to-sample distances
		if (!is.null(saved_plots$SampleDistance)){
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
			plot.new()
			}
		}
		## Heatmap
		if (!is.null(saved_plots$pheatmap2)){
			grid.draw(saved_plots$pheatmap2)
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
				plot.new()
			}
		}
		###vennDiagram
		if (!is.null(saved_plots$vennDiagram)){
			for (i in 1:length(saved_plots$vennDiagram)) {
				title_tmp <- paste("vennDiagram (",i,")",sep="")
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

