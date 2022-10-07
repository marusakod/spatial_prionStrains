

library(Seurat)
library(tidyverse)
library(shiny)
library(patchwork)
library(strex)
library(DT)

all_spatial <- readRDS("all_spatial_seurat_normalized.rds")
allen_cortex <-  readRDS("allen_cortex.rds")
Brain_regions_DE_results <- readRDS("Brain_regions_DE_results.rds")
all_SVFs <- readRDS("all_SVFs.rds")
spatial_merged <- readRDS("spatial_merged.rds")
Brain_regions_DE_withinSlice_res <- readRDS("Brain_regions_DE_withinSlice_res.rds")
all_ORA_results_between_strains_dfs <- readRDS("all_ORA_results_between_strains_dfs.rds")
all_ORA_results_within_strains_dfs <- readRDS("all_ORA_results_within_strains_dfs.rds")

all_ORA_results_between_strains_dfs_merged <- bind_rows(all_ORA_results_between_strains_dfs)
all_ORA_results_within_strains_dfs_merged <- bind_rows(all_ORA_results_within_strains_dfs)


all_spatial_gene_symbols <- rownames(GetAssayData(all_spatial$RML6, assay = "Spatial", slot ="counts"))
all_reference_symbols <- rownames(allen_cortex)
simages <- Images(spatial_merged)
names(simages) <- names(all_spatial)


# all unique brain regions in DE results

Brain_regions <- sapply(Brain_regions_DE_results, FUN = names, simplify = TRUE)
Brain_regions <- sapply(Brain_regions, str_after_nth, pattern = "_", n = 2, simplify = TRUE)
Brain_regions <- Reduce(unique, Brain_regions)

# add comparison and brain regions to DE results

complete_DE_results <- function(DE_list){
  names <- names(DE_list)
  comparisons <- str_before_nth(names, pattern = "_", n =2 )
  regions <- str_after_nth(names, pattern = "_", n = 2)
  
  for(i in 1:length(DE_list)){
    DE_list[[i]]$comparison <- comparisons[i]
    DE_list[[i]]$region <- regions[i]
    DE_list[[i]]$gene <- rownames(DE_list[[i]])
  }
  
  DE_list
}

Brain_regions_DE_results_completed <- sapply(Brain_regions_DE_results, 
                                             FUN = complete_DE_results, 
                                             simplify = FALSE, 
                                             USE.NAMES = TRUE)

Brain_regions_DE_results_completed_merged <- unlist(Brain_regions_DE_results_completed, recursive = FALSE)
Brain_regions_DE_results_completed_merged <- bind_rows(Brain_regions_DE_results_completed_merged) %>% as_tibble()
Brain_regions_DE_results_completed_merged$comparison <- gsub("_", " vs ", Brain_regions_DE_results_completed_merged$comparison)


# get unique brain regions for each strain (used as conditional select input for within slide DEA results tab)
strainSpecific_brainregions <- names(Brain_regions_DE_withinSlice_res)
names(strainSpecific_brainregions) <- gsub("[.].*", "", strainSpecific_brainregions)
strainSpecific_brainregions <- gsub(".*[.]", "", strainSpecific_brainregions)
strainSpecific_brainregions_list <- sapply(names(all_spatial), FUN = function(x){ 
                                     y <- strainSpecific_brainregions[names(strainSpecific_brainregions) == x]
                                     y <- unname(y)
                                     y
                                     }
                                     , simplify = FALSE)



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
           )), 
  tabPanel("Brain Regions DE", 
           value = "regions_DE", 
           sidebarLayout(
             sidebarPanel = sidebarPanel(width = 2,
                                         selectInput("selectComparison", 
                                                     "Select Comparison", 
                                                     choices = unique(Brain_regions_DE_results_completed_merged$comparison), 
                                                     selected = "NBH vs RML6", 
                                                     multiple = FALSE, 
                                                     width = "200px"), 
                                         selectInput("selectRegion", 
                                                     "Select brain region", 
                                                     choices = Brain_regions, 
                                                     selected = "Thalamic_Nuclei", 
                                                     multiple = FALSE, 
                                                     width = "200px")), 
             mainPanel = mainPanel(width = 10, 
                                   plotOutput("comparison1_plot"), 
                                   plotOutput("comparison2_plot"), 
                                   dataTableOutput("DE_table"), 
                                   dataTableOutput("GO_results_between"))
           )),
  tabPanel("Brain Regions within slice DE", 
           value = "within_slice_regions_DE", 
           sidebarLayout(
             sidebarPanel = sidebarPanel(with = 2, 
                                         radioButtons("selectStrain", 
                                                      "select Strain", 
                                                      choices = names(all_spatial), 
                                                      selected = "NBH", 
                                                      width = "200px"), 
                                         uiOutput("selectRegionOneSlice")), 
             mainPanel = mainPanel(dataTableOutput("withinSlice_DEgenes"), 
                                   dataTableOutput("GO_results_within")), 
                                   
           )
           
  ), 
  tabPanel("Module scores", 
           value = "module_scores", 
           sidebarLayout(
             sidebarPanel = sidebarPanel(width = 2, 
                             radioButtons("selectGO", 
                                          "Select_Gene_Ontology_term", 
                                          choices = c("Neuron Apoptotic Process" = "neuron_apoptotic_process1", 
                                                      "Neuroinflammatory response" = "neuroinflammatory_response2"), 
                                          width = "200px")), 
             mainPanel = mainPanel(plotOutput("moduleScores_plot", width = "1500px", height = "800px"))
           )
           ), 
  tabPanel("Spatially variable features", 
           value = "SVFs", 
           sidebarLayout(
             sidebarPanel = sidebarPanel(width = 2, 
                                         radioButtons("SVF_strain", 
                                                      "Select strain", 
                                                      choices = names(all_spatial), 
                                                      selected = "NBH", 
                                                      width = "200px")), 
             mainPanel = mainPanel(dataTableOutput("SVF_table"), 
                                   plotOutput("SVF_plot", width = "800px", height = "700px"))
           )
           )
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


# DIFFERENTIAL EXPRESSION RESULTS FOR BRAIN REGIONS TAB

comparison1 <- reactive({
  str_before_first(input$selectComparison, pattern = " vs ")
})

comparison2 <- reactive({
  str_after_first(input$selectComparison, pattern = " vs ")
})


# plot with brain regions colored in each strain in comparison

spatial_merged_reactive <- reactive({
  spatial_merged@meta.data <- spatial_merged@meta.data %>% mutate(Selected_region = case_when(Brain.Regions_new == input$selectRegion ~ input$selectRegion, 
                                                                                             TRUE ~ "other"))
  spatial_merged
})


output$comparison1_plot <- renderPlot({
 
  SpatialDimPlot(spatial_merged_reactive(), images = simages[comparison1()], group.by = "Selected_region") + theme(legend.position = "none")
})

output$comparison2_plot <- renderPlot({
  SpatialDimPlot(spatial_merged_reactive(), images = simages[comparison2()], group.by = "Selected_region") + theme(legend.position = "none")
  
})


# filter DE results based on the region and comparison selected
output$DE_table  <- renderDataTable({
  
  Brain_regions_DE_results_completed_merged %>% filter(comparison == input$selectComparison, region == input$selectRegion)
  
})


output$GO_results_between <- renderDataTable({
  
  all_ORA_results_between_strains_dfs_merged %>% filter(comparison == gsub(" vs ", "_", input$selectComparison), brainRegion == input$selectRegion)
  
})




# BRAIN REGIONS WITHIN SLICE DE TAB

# populate select input with brain regions that exist in selected strain
output$selectRegionOneSlice <-renderUI({
  
  selectInput("RegionOneSlice", 
              "Select brain region", 
              choices = strainSpecific_brainregions_list[[input$selectStrain]], 
              multiple = FALSE, 
              width = "200px")
})

# get the correct table from DE results based on brain region and strain selected

output$withinSlice_DEgenes <- renderDataTable({
  Brain_regions_DE_withinSlice_res[[paste(input$selectStrain, input$RegionOneSlice, sep = ".")]]
})

# table with within slice DE results

output$GO_results_within <- renderDataTable({
  
  all_ORA_results_within_strains_dfs_merged %>% filter(comparison == input$selectStrain, brainRegion == input$RegionOneSlice)
  
})




# MODULE SCORE TAB 

output$moduleScores_plot <- renderPlot({
  
  SpatialFeaturePlot(spatial_merged, features = input$selectGO)
  
})



# SPATIALLY VARIABLE FEATURES TAB

output$SVF_table <- renderDataTable({
  
  all_SVFs %>% filter(Strain == input$SVF_strain)
  
})

# display top 10 spatially variable features

output$SVF_plot <- renderPlot({
  
  top_10_svf <- head(all_SVFs %>%
                       filter(Strain == input$SVF_strain) %>% 
                       dplyr::select("gene") %>% 
                       flatten_chr(),
                       n = 10)

  SpatialFeaturePlot(all_spatial[[input$SVF_strain]], features = top_10_svf, ncol = 3, alpha = c(0.1, 1))
  
})





} 
  
 shinyApp(ui = ui, server = server)
 
 