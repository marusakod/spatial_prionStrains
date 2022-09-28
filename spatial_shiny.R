

library(Seurat)
library(tidyverse)
library(shiny)

all_spatial <- readRDS("all_spatial_seurat_normalized.rds")
allen_cortex <-  readRDS("allen_cortex.rds")
all_spatial_gene_symbols <- rownames(GetAssayData(all_spatial$RML6, assay = "Spatial", slot ="counts"))
all_reference_symbols <- rownames(allen_cortex)

ui <- fluidPage(tabsetPanel(
  tabPanel("Gene search", 
           value = "gene_search", 
           sidebarLayout(
             sidebarPanel = sidebarPanel(width = 2, 
        
                                         selectizeInput("geneSymbolSpatial",
                                                        label = "Select Gene Symbol", 
                                                        multiple = FALSE,
                                                        size = 20,
                                                        choices = NULL,
                                                        width = "200px",
                                                        #options = list(
                                                        # placeholder = 'eg. Agap2',
                                                        #   onInitialize = I('function() { this.setValue("");}'))
                                         ), 
                                         numericInput("alphaSpatialImage",
                                                      label = "Select image transparency", 
                                                      value = 1, 
                                                      min = 0,
                                                      max = 1,
                                                      step = 0.1, 
                                                      width = "200px"), 
                                         numericInput("alphaSpatialSpot",
                                                      label = "Select spot transparency", 
                                                      value = 1, 
                                                      min = 0,
                                                      max = 1,
                                                      step = 0.1, 
                                                      width = "200px"), 
                                         numericInput("sizeSpatialSpot",
                                                      label = "Select spot size", 
                                                      value = 1.2, 
                                                      min = 0.1,
                                                      max = 1.5,
                                                      step = 0.1, 
                                                      width = "200px")),
             mainPanel = mainPanel(width = 10, 
                                   uiOutput("get_plot"), 
                                   plotOutput("spatial_plot", width = "1300px", height = "3000px"))
             )), 
  tabPanel("scRNAseq-Spatial integration", 
           value = "integration", 
           sidebarLayout(
             sidebarPanel = sidebarPanel(width = 2, 
                                         radioButtons("strainToIntegrate", 
                                                      "Select strain", 
                                                      choices = names(all_spatial), 
                                                      selected = "NBH", 
                                                      width = "200px"), 
                                         selectInput("cellToDisplay", 
                                                     "Select celltype", 
                                                     choices = unique(allen_cortex@meta.data$subclass), 
                                                     selected = "Macrophage", 
                                                     multiple = FALSE, 
                                                     width = "200px"), 
                                         selectizeInput("geneInRef", 
                                                        "Select gene to display in reference data", 
                                                        choices = NULL, 
                                                        multiple = FALSE, 
                                                        size = 20, 
                                                        width = "200px")
                                         ), 
             mainPanel = mainPanel(width = 10, 
                                   plotOutput("refPlot", width = "800px", height = "700px"), 
                                   plotOutput("cellInSpatial", width = "1200px", height = "1200px"))
           ))
  )
)



# update select imput for spatial

 server <- function(input, output, session){

   #### GENE SEARCH TAB ###############################################

updateSelectizeInput(session  = session,
                     inputId  = "geneSymbolSpatial",
                     server   = TRUE,
                     choices  = all_spatial_gene_symbols,
                     selected = character(0),
                     options = list(
                     placeholder = 'eg. Gfap',
                     onInitialize = I('function() { this.setValue("");}')
))


output$spatial_plot <- renderPlot({
  
  if(input$geneSymbolSpatial == ""){
    p <- ggplot() +
      theme_void() +
      geom_text(aes(0,0,label= "")) +
      xlab(NULL)
  }else{
  f <- function(srat, strain){
    SpatialFeaturePlot(srat,
                       features = input$geneSymbolSpatial, 
                       pt.size.factor = input$sizeSpatialSpot, 
                       alpha = input$alphaSpatialSpot, 
                       image.alpha = input$alphaSpatialImage) + ggtitle(strain) + theme(plot.title = element_text(size = 16))
  }
  
  p_list <- mapply(srat = all_spatial, 
              strain = names(all_spatial), 
              FUN = f, 
              SIMPLIFY = FALSE, 
              USE.NAMES = TRUE)
  
  p <- wrap_plots(p_list[[1]], p_list[[2]], p_list[[3]], p_list[[4]], ncol = 1)
  
  }
  
  p
  
})


#### SINGLE CELL - SPATIAL INTEGRATION TAB ################################

# update selectize input for displaying genes in reference UMAP

updateSelectizeInput(session  = session,
                     inputId  = "geneInRef",
                     server   = TRUE,
                     choices  = all_reference_symbols,
                     selected = character(0),
                     options = list(
                       placeholder = 'eg. Gfap',
                       onInitialize = I('function() { this.setValue("");}')
                     ))

# when integration tab is selected set default Seurat assay to predictions
observeEvent(input$integration, {
  all_spatial <- sapply(all_spatial, FUN = function(x){DefaultAssay(x) <- "predictions"
  x}, 
  simplify = FALSE, 
  USE.NAMES = TRUE)
})

# when gene search tab is selected set default Seurat assay to "SCT


observeEvent(input$gene_search, {
  all_spatial <- sapply(all_spatial, FUN = function(x){DefaultAssay(x) <- "SCT"
  x}, 
  simplify = FALSE, 
  USE.NAMES = TRUE)
})

output$refPlot <- renderPlot({
  
  if(input$geneInRef == ""){
    p <- DimPlot(allen_cortex, group.by = "subclass", label = TRUE)
  }else{
    p <- FeaturePlot(allen_cortex, features = input$geneInRef)
  }
  
  p

})

output$cellInSpatial <- renderPlot({
  SpatialFeaturePlot(all_spatial[[input$strainToIntegrate]], features = input$cellToDisplay , pt.size.factor = 1.6, ncol = 2, crop = TRUE)
})




} 
  
 shinyApp(ui = ui, server = server)
 
 