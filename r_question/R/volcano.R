#'volcano.R 
#'
#'Shiny Module for building an interactive plotly volcano plot in a shiny 
#'dashboard box.

library(shiny)
library(manhattanly)
library(shinydashboard)
library(shinydashboardPlus)
library(plotly)
library(shinyjqui)



#'Helper function for formatting mappings between term ids and genes associated 
#'with them.
.format_intersection <- function(intersection_string){
    return(as.numeric(strsplit(intersection_string, split = ",")[[1]]))
}

volcano_ui  <- function(id, label = NULL, color = NULL, study_name = NULL){
    ns <- NS(id)
    
    div(tags$style(HTML(paste0( "#", id, " .box-header.with-border", "{ background-color:", color, " !important;}"))),
    
    #enabling draggable & resizable box using. 
    #Need to specify header as a handle by using the class id
    jqui_draggable(options = list(handle = ".box-header.with-border", containment = "#plot_space"), 
        jqui_resizable(options =  list(minWidth = 300, maxWidth = 1000, maxHeight = 461, minHeight = 461 ),
        box(
            title = paste0(study_name, " Gene Expression Volcano Plot"),
            width = 6,
            collapsible = T,
            closable = T,
            sidebar = boxSidebar(
                id = ns("myboxsidebar"),
                width = 40,
                selectInput(
                    ns("term_type"),
                    label = "Enrichment Type",
                    choices = list("Biological Process" = "GO:BP", 
                                   "Cell Component" = "GO:CC",
                                   "Molecular Function" = "GO:MF",
                                   "Transcription Factor Motif" = "TF",
                                   "KEGG" = "KEGG",
                                   "Reactome"= "REAC",
                                   "All" = "all"),
                    selected = "all")
                ),
            plotlyOutput(ns("plot"))
            )
        )
    )
)  
    
   
    
}

volcano_server <- function(id, gprofiler, dataset, label = NULL){
    moduleServer(id, 
                 function(input, output, session){
                    
                     ns <- session$ns
                     
                     #If I include these lines and reference the inputs outside of the reactive context, I get the desired result.
                     #data <- dataset
                     #gp <- gprofiler
                     
                     filtered_data <- reactive({
                                                if(input$term_type != "all"){
                         
                                                    genes <- dplyr::filter(gprofiler$result, 
                                                                       source == input$term_type)$intersection
                                                    genes <- unique(unlist(sapply(genes, .format_intersection, USE.NAMES = F)))
                                                
                                                    dataset[which(dataset$Entrez.ID %in% genes),] %>% na.omit
                                                    } else {
                                                        dataset %>% na.omit
                                                    }
                         
                                                })
                     
                     output$plot <- renderPlotly(
                         volcanoly(filtered_data(), gene = "Entrez.ID", 
                                   effect_size = "logFC", 
                                   p = "adj.P.Val", 
                                   genomewideline = -log10(1e-2), effect_size_line = FALSE, 
                                   snp = "Gene.symbol", annotation1 = "adj.P.Val", title = NULL))
                     
                     #placeholder return Value to make sure shiny knows its there when rendered
                     return(TRUE)
                    })
}


#' test ui and server
ui <- dashboardPage( title = "volcano-test",
                     #wrapping everything in a dashBoard page makes the box work                
                     dashboardHeader(),
                     dashboardSidebar(),
                     dashboardBody(
                         tags$style("body { background-color: ghostwhite}"),
                         volcano_ui("plot")
                     )
                     
)


server <- function(input, output) {
    #example data    
    data <- read.csv("datasets\\GSE4136_Scer\\GSE4136_25thGen.csv")
    fdata <- data %>% dplyr::filter(adj.P.Val <= 0.1) %>% dplyr::filter(logFC >= 1 || logFC <= -1)
    
    up <- fdata %>% filter(logFC >= 1) %>% na.omit 
    down <- fdata %>% filter(logFC <= -1) %>% na.omit
    
    gp <- gprofiler2::gost(list("up-regulated" =  as.numeric(up$Entrez.ID), "down-regulated" = as.numeric(down$Entrez.ID)), 
               organism = "scerevisiae", numeric_ns = "ENTREZGENE_ACC"
               ,ordered_query = T, multi_query = F, evcodes = T)
    
    
    volcano_server("plot", gp, data)
}

# Run the application 
shinyApp(ui = ui, server = server)
