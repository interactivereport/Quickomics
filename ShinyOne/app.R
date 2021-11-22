#This app enables users to create an online table listing all QuickOmics Projects with one click link to launch each project.
#Edit the project.csv file in this folder to add/remove projects. 
#The URLs for projects should be QuickOmics URL of your instance, followed by /?serverfile=projectID, or /?testfile=projectID. Use the config.csv file in QuickOmics directory to define server file and test file folders.

library(shiny)
library(DT)
library(dplyr)

createLink <- function(val) {
  sprintf('<a href="%s" target="_blank" class="btn btn-primary">View Project</a>',val)
}

ui <- fluidPage(  
  titlePanel("Datasets Loaded into Quickomics"),
  tags$hr(style="border-color: RoyalBlue;"),
      dataTableOutput('table1')

)

server <- function(input, output) {
  
  output$table1 <- renderDataTable({
    my_table=read.csv("projects.csv", check.names=F)
    my_table<-my_table%>%dplyr::mutate(Link=createLink(URL) )%>%dplyr::select(-URL)
    #browser() #debug
    return(my_table) }
   , escape = FALSE)
}

shinyApp(ui, server)
