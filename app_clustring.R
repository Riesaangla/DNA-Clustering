
#GD-10
rm(list=ls())
# source("3-functions.R")
library(rsconnect)
library(shiny)
library(markdown)
library(shinyWidgets)
library(dplyr)
library(ggplot2)
library(ROCit)
library(rlist)



# User interface ----

ui <- fluidPage(
  
  navbarPage(title="DNA Clustering",
             
            tabPanel("Information of DNA Clustering",
                     column(1),
                     column(5,br(),br(),br(),
                            withMathJax(p("Here is an introduction of the app."),
                                        p("Here is the second paragraph pf the introduction of this app."))
                            ),
                     column(1)),
  
            tabPanel("Data Exploration",
                     
                     tabsetPanel(
                       
                       tabPanel("Sample Data",
                                fluidRow(
                                  column(3,
                                         wellPanel(
                                           selectInput("Sample_DT", "Choose one sample Data:", choices=list("Sample-1"=1,"Sample-2"=2), selected=1),
                                           # bsPopover(id="sampdat", title="Data set information", content="The one-sample data set is collected from the Old Faithful geyser in Yellowstone National Park, Wyoming, USA.  The variable of interest is the duration of each eruption in minutes.  The two-sample data set is associated with the 1974 Motor Trend US magazine.  The explanatory variable is the type of transmission and the response variable is horsepower.",
                                           #           trigger="hover",placement="right"),
                                           checkboxInput("usesample_DT", "Get result by using sample data", TRUE),
                                           tags$hr(),
                                           helpText("Remember to uncheck this when not using sample data!"),
                                           br(),br(),br())),
                                  
                                  column(9,
                                         conditionalPanel(
                                           condition="input.usesample_DT",
                                           tableOutput("Sample_DT"))))),
                       
                       tabPanel("Upload Data",
                                
                                column(3,wellPanel(
                                  fileInput("file1", label = h3("Choose CSV File"),multiple = FALSE,
                                            accept = c("text/csv/txt", "text/comma-separated-values,text/plain",".csv"),
                                            placeholder = "No file selected"),
                                  tags$hr(),
                                  helpText("Please upload a csv file (make sure there are only one sheet in the file)."),
                                  )),
                                
                                column(9,conditionalPanel(
                                  condition="input.file1!='NULL'",
                                  tableOutput("input_DT")
                                ))),
                       
                       tabPanel("Result",column(3,tableOutput("output_DT")))
                     )
            )
  )
)




# Define server logic ----


server <- function(input, output) {
  
  output$input_DT <- renderTable({
    
    inFile <- input$file1
    
    if (is.null(inFile))
      return(NULL)
    
    Dat <- readxl::read_excel(inFile$datapath,1)[,-2]
    
    return(Dat)
    
  })
  
  output$output_DT <- renderTable({
    
    inFile <- input$file1
    
    if (is.null(inFile))
      return(NULL)
    
    Dat <- readxl::read_excel(inFile$datapath,1)[,-2]
    
    sort <- data.frame("Row"=c(1:nrow(Dat)),Dat$Sample, "group"= rep(0,nrow(Dat)))
    
    sort_new <- sort
    # for a few (<= 3) allele values, group it to Group 0;
    sort_new$group <- as.numeric( apply(Dat[,-1], 1, FUN=function(x) sum(!is.na(x))) > 3 )
    ################ Exclusion ################
    while (!identical(sort_new$group,sort$group)){
      sort <- sort_new
      ref_g <- max(sort$group)
      list_g <- which(sort$group==ref_g)
      ref_r <- list_g[which.max(apply(Dat[list_g,], 1, FUN=function(x) sum(!is.na(x))))]# reference sample row in group 
      for  (i in list_g){
        m<- Dat[c(i,ref_r),-1]
        m_ref <- m[,colSums(is.na(m)) == 0]
        if (length(m_ref) == 0){ # if no variables in common, then belongs to remaining group
          sort_new$group[i] = sort_new$group[i] +1
        } else{
          if (sum(m_ref[1,] != m_ref[2,]) != 0){# if 2 rows(after excluding NA columns) are diff, they belong to the remaining groups
            sort_new$group[i] = sort_new$group[i] +1
          }
        } 
      }
    }
    #sort_new$group <- ifelse(sort_new$group==0,'Too few allele values',sort_new$group)
    return(sort_new[order(sort_new$group),])
    
  })

  
}
  


# Run the app ----
shinyApp(ui, server)


