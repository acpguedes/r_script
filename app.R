#!/usr/bin/R
check_lib <- function(mypkg, install_pkg=1, ...){
  if(is.element(mypkg, installed.packages()[,1])){
    suppressMessages(library(mypkg, character.only=TRUE)) #, envir = "/usr/local/lib/R/site-library"))
  } else {
    if(install_pkg == 1){
      install.packages(mypkg, character.only=TRUE, repos="http://vps.fmvz.usp.br/CRAN/")
      suppressMessages(library(mypkg, character.only=TRUE)) #, envir = "/usr/local/lib/R/site-library"))
    } else {
      stop(paste("Missing package:", mypkg, sep=" "), call.=FALSE)
    }
  }
}

library(shiny)
#library(httpuv)
#library(xtable)
#library(htmltools)
library(shinythemes)
library(shinydashboard)
library(tidyverse)
library(DT)
#check_lib("plotly")

#myFile <- '/home/acpguedes/projects/sig_trans/work/SBP_5/20190730/test2/sbps_node_annotation_leiden.tsv'
myFile <-"sbps_node_annotation_leiden.tsv"
hh_raw_leiden <- dplyr::as_tibble(data.table::fread(myFile))

tibble(
  Cluster      = c("A", "B", "C", "D", "D G", "E", "F", "F E", "G", "B A"),
  ClusterColor = c("red", "blue", "purple", "black", "green", "yellow", "orange", "pink", "violet", "gray")
) -> mycolors

hh_raw_leiden %>% 
  filter(from_genomes) %>% 
  #left_join(mycolors) %>% 
  select(node.name, id_0cov_0.8evalue_0.001, 
         full, Cluster, 
         TMSP = `periprate1e-3`, cache_rate,  
         groupsize, hhArch, `archConsensus1e-3` #, ClusterColor 
  ) %>% 
  distinct() -> tmpdata

tmpdata_mapping <- data.frame(
  Species = tmpdata$node.name,
  row_id = 1:length(tmpdata$node.name),
  stringsAsFactors = FALSE
)

ui <- fluidPage(
  fluidRow(
      plotOutput("plot1", 
                 click = "plot1_click",
                 brush = brushOpts(
                   id = "plot1_brush")
      )
  ),
  
  fluidRow(
    div(
      p(
            'Selecione uma area com circulos para ver os detalhes abaixo'
      )
    
  ),
  
  dashboardBody(
    
    fluidRow(
      div( 
        style = 'overflow-x: scroll', 
        DT::dataTableOutput( "table", width = "100%" ) 
      )
    )
  )
  )
)


server <- function(input, output) {
  
  output$plot1 <- renderPlot({
    p <- tmpdata %>% 
      ggplot(
        aes(
          TMSP, 
          cache_rate, 
          group = as_factor(Cluster), 
          color = as_factor(Cluster),
          label1 = id_0cov_0.8evalue_0.001,
          label2 = groupsize
        )
      ) + 
      geom_point(aes(size  = log10(groupsize+1), alpha = 0.3)) + #aes(shape = as_factor(Cluster))
      geom_rug(col=rgb(.5,0,0,alpha=.2)) +
      scale_color_manual(values = deframe(mycolors)) +
      facet_wrap( . ~ full)  +
      labs(
        title = "Relation between SBP and Cache",
        subtitle = "Relationship between localization of a SBP and neighborhood with a Cache",
        y = "Frequency of SBP neighbors of a Cache",
        x = "Frequency of SBP with TM and/or SP"
      ) +
      scale_size_discrete(name = "Poolman cluster") +
      scale_size_continuous(name = "Cluster size (log scale)") +
      scale_y_continuous(expand = expansion(mult = .11)) +
      scale_x_continuous(expand = expansion(mult = .11)) +
      scale_alpha(guide = FALSE) +
      theme_classic() +
      theme(axis.text.x = element_text(angle=-45))
    #ggplotly(p)
    p
  })
  
  
  dat <- reactive({
    user_brush <- input$plot1_brush#
    
    bP <- brushedPoints(tmpdata, user_brush, xvar = "TMSP", yvar = "cache_rate")
    
    hh_raw_leiden %>% 
      filter(node.name %in% bP$node.name ) %>% 
      select(-e0.1, -e0.001, -e5, -e6, -e10, -e15, -e20) %>% 
      arrange(desc(is.close))
  })
  output$table <- DT::renderDataTable({DT::datatable(dat(), 
                                                     filter = 'top', options = list(autoWidth = T, scrollX = T))})
  
}

shiny::shinyApp(ui, server)


