
library(shiny)
library(shinydashboard)
library(shinydashboardPlus)
library(shinyjqui)
library(gprofiler2)
library(shinyWidgets)
library(dplyr)
library(shinyjs)
library(DT)
library(shinyjqui)
library(purrr)



#'building dynamic sidebar menu
#'
build_menu <- function(input_data){
    tags_out <- list()
    for(i in 1:length(names(input_data))){
        study <- names(input_data)[i]
        tags_out[[i]] <-  menuItem(text = study, 
                                           menuSubItem(tabName =  paste0(study,"_gexp_tab"), text = paste0("Gene Expression"), icon = icon("line-chart")),
                                           menuItem(text = paste0("Gene Functional Groups"), 
                                                    menuSubItem(tabName =  paste0(study, "_manhat_tab"), text  = "Manhattan Plot", icon = icon("line-chart")),
                                                    menuSubItem(tabName =  paste0(study, "_biofab_tab"), text = "BioFabric Plot", icon = icon("line-chart")))
                                    )
                                  
    }
    
    return(tags_out)
}



# Processing GO enrichment via gprofiler2
process_data <- function(dataset, GSEcode, Organism){
    require(gprofiler2, dpylr)
    
    filtered <- dataset %>% dplyr::filter(adj.P.Val <= 0.1) %>% 
        dplyr::filter(logFC >= 1 || logFC <= -1)
    up <- filtered %>% dplyr::filter(logFC >= 1) %>% na.omit 
    down <- filtered %>% dplyr::filter(logFC <= -1) %>% na.omit 
    
    gp <- gprofiler2::gost(list("up-regulated" =  as.numeric(up$Entrez.ID), "down-regulated" = as.numeric(down$Entrez.ID)), 
                           organism = Organism, numeric_ns = "ENTREZGENE_ACC"
                           ,ordered_query = T, multi_query = F, evcodes = T)
    outname <- paste0(GSEcode, ": ", Organism)
    return(list( "data" = list("gene_expression" = dataset,"gp_result" = gp), "name" = outname))
    
    
}

ui <- dashboardPage( title = "dash-test",
    dashboardHeader(),
    dashboardSidebar(
        sidebarMenu(menuItemOutput("menuItem"))
        ),
    dashboardBody(
        uiOutput(outputId = "tabs"))

   
)


server <- function(input, output) {
    
    data_scer <- read.csv("test_data\\GSE4136_25thGen.csv")
    data_cele <- read.csv("test_data\\GSE27338.csv")
    
    
    
    scer <- process_data(data_scer, "GSE4136", "scerevisiae")
    cele <- process_data(data_cele, "GSE27338", "celegans")
    
    input_data <- list("GSE4136" = scer,"GSE27338" = cele)
    
    #DEFINING COLOR PALETTE AND ALL POSSIBLE PLOT IDS
    palette <- suppressWarnings(RColorBrewer::brewer.pal(n = length(input_data), name = "Pastel2"))
    names(palette) <- names(input_data)
    
    
    #BUILDING SIDEBAR MENU
    output$menuItem <- renderMenu({sidebarMenu(.list =  build_menu(input_data))})
    
    #BUILDING TABSETS
    output$tabs <- renderUI({
        all_tabs <- list()
        
        for(i in 1:length(input_data)){
            study <- names(input_data)[i]
            #'volcano plot definitions
            all_tabs[[paste0(study,"_gexp_tabi")]] <- tabItem(tabName = paste0(study, "_gexp_tab"),
                                                            volcano_ui(id = paste0(study, "_gexp"), 
                                                                       color = palette[[study]], 
                                                                       study_name = study))
            
            #all_tabs[[paste0(study,"_manhat_tabi")]] <- tabItem(tabName = paste0(study, "_manhat_tab"),
            #                                                   manhat_ui(id = paste0(study, "_manhat"), 
            #                                                              color = palette[[study]], 
            #                                                              study_name = study))
            
            #all_tabs[[paste0(study,"_biofab_tabi")]] <- tabItem(tabName = paste0(study, "_biofab_tab"),
            #                                                   go_fabric_ui(id = paste0(study, "_biofab"), 
            #                                                             color = palette[[study]], 
            #                                                             study_name = study))
        }
        return(tags$div(class = "tab-content", all_tabs))
    })
    
    #MODULE DEFINITIONS
    module_outputs <- reactiveValues()
    
    for(i in 1:length(input_data)){
        study <- names(input_data)[i]
        module_outputs[[paste0(study,"_gexpi")]] <- volcano_server(id = paste0(study,"_gexp"), 
                                                                  gprofiler = input_data[[study]]$data$gp_result,
                                                                  dataset = input_data[[study]]$data$gene_expression
                                                                  )
        #module_outputs[[paste0(study,"_manhat")]] <- manhat_server(id = paste0(study,"_manhat"), 
          #                                                        gprofiler = input_data[[study]]$data$gp_result
                                                                  
         #                                                        )
        #module_outputs[[paste0(study,"_biofab")]] <- go_fabric_server(id = paste0(study,"_biofab"), 
         #                                                          gprofiler = input_data[[study]]$data$gp_result,
         #                                                          dataset = input_data[[study]]$data$gene_expression
         #                                                        )
        
    }
    

}

# Run the application 
shinyApp(ui = ui, server = server)
