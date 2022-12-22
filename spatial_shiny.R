

library(Seurat)
library(tidyverse)
library(shiny)
library(patchwork)
library(strex)
library(DT)
library(RColorBrewer)

# all seurat objects 
all_spatial_seurat_normalized <- readRDS("all_spatial_seurat_normalized.rds")
seurat_merged_timepoints <- readRDS("seurat_merged_timepoints.rds")

# all diffeentially expressed genes
DEGs_conditions <- readRDS("all_timepoint_condition_region_DEGs_merged.rds")
DEGs_regions <- readRDS("Brain_regions_DE_withinSlice_res_merged.rds")

# all ORA results for DEGs

ORA_results_between_strains_dfs_merged <- readRDS("ORA_results_between_strains_dfs_merged.rds")
ORA_results_within_strains_dfs_merged <- readRDS("ORA_results_within_strains_dfs_merged.rds")
                                                                     

# all spatially variable features 
all_SVFs <- readRDS("all_SVFs.rds")
all_SVFs_merged <- bind_rows(all_SVFs)

# get all gene symbols for select imput in Gene Search tab
all_spatial_gene_symbols <- rownames(GetAssayData(all_spatial_seurat_normalizedNBH_30, assay = "Spatial", slot ="counts"))


simages <- mapply(x = seurat_merged_timepoints, name = names(seurat_merged_timepoints), FUN =  function(x, name){
  I <- Images(x)
  names(I) <- paste(c("NBH", "RML6"), name, sep = "_")
  I
}, 
SIMPLIFY  = FALSE, 
USE.NAMES = TRUE)



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
  
  tabPanel("Brain Regions DE",
           value = "regions_DE",
           sidebarLayout(
             sidebarPanel = sidebarPanel(width = 2,
                                         selectInput("select_BS_DEA_timepoint",
                                                     "Select timepoint",
                                                     choices = c("27", "30", "TS"),
                                                     selected = "27",
                                                     multiple = FALSE,
                                                     width = "200px"),
                                         uiOutput("selectBetweenSlice")
             ),
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
                                                      "Select Strain",
                                                      choices = c("NBH", "RML6"),
                                                      selected = "NBH",
                                                      width = "200px"),
                                         radioButtons("select_WS_DEA_timepoint",
                                                      "Select Timepoint",
                                                      choices = c("27", "30", "TS"),
                                                      selected = "27",
                                                      width = "200px"),
                                         uiOutput("selectRegionOneSlice"),
                                         uiOutput("selectComparisonRegion")),
             mainPanel = mainPanel(dataTableOutput("withinSlice_DEgenes"),
                                   dataTableOutput("GO_results_within")),
             
           )
           
  ),
  tabPanel("Module scores",
           value = "module_scores",
           sidebarLayout(
             sidebarPanel = sidebarPanel(width = 2,
                                         radioButtons("selectGO",
                                                      "Select Gene Ontology term",
                                                      choices = c("Neuron Apoptotic Process" = "neuron_apoptotic_process",
                                                                  "Neuroinflammatory response" = "neuroinflammatory_response"),
                                                      width = "200px")),
             mainPanel = mainPanel(plotOutput("moduleScores_plot", width = "1300px", height = "3000px"))
           )
  ),
  tabPanel("Spatially variable features",
           value = "SVFs",
           sidebarLayout(
             sidebarPanel = sidebarPanel(width = 2,
                                         radioButtons("SVF_condition",
                                                      "Condition",
                                                      choices = c("NBH", "RML6"),
                                                      selected = "NBH",
                                                      width = "200px"),
                                         radioButtons("SVF_timepoint",
                                                      "Time point",
                                                      choices = c("27", "30" ,"TS"),
                                                      selected = "27",
                                                      width = "200px")),
             mainPanel = mainPanel(dataTableOutput("SVF_table"),
                                   plotOutput("SVF_plot", width = "1300px", height = "3000px"))
           )
  )
)
)

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
  f <- function(srat, images){
    # get maximum expression of a gene in all strains
    
    all_gene_counts <- sapply(seurat_merged_timepoints, FUN = function(x){ 
      mx <- GetAssayData(x, assay = "SCT")
      mx2 <- GetAssayData(x, assay = "Spatial")
      
      if(input$geneSymbolSpatial %in% rownames(mx)){
        y <- mx[input$geneSymbolSpatial, ]
      }else{
        y <- mx2[input$geneSymbolSpatial, ]
      }
      y
      
      }, simplify = FALSE)
    
    all_gene_counts_merged <- unlist(all_gene_counts)
    max  <- max(all_gene_counts_merged)
    
    
    NBH_plot <- SpatialFeaturePlot(srat,
                       features = input$geneSymbolSpatial, 
                       pt.size.factor = input$sizeSpatialSpot, 
                       alpha = input$alphaSpatialSpot, 
                       image.alpha = input$alphaSpatialImage, 
                       images = images[1]) + 
                       
      ggtitle(names(images[1])) + 
      theme(plot.title = element_text(size = 18), 
            legend.title = element_text(size = 16), 
            legend.text = element_text(size = 14)) +
      scale_fill_gradientn(colours= rev(brewer.pal(name = "Spectral", n = 11)),
                           breaks=c(0,0.5,1),
                           limits=c(0,max)) +
      guides(fill=guide_colourbar(barwidth=20, barheight = 2))
    
    RML6_plot <- SpatialFeaturePlot(srat,
                                    features = input$geneSymbolSpatial, 
                                    pt.size.factor = input$sizeSpatialSpot, 
                                    alpha = input$alphaSpatialSpot, 
                                    image.alpha = input$alphaSpatialImage, 
                                    images = images[2]) + 
      
      ggtitle(names(images[2])) + 
      theme(plot.title = element_text(size = 18), 
            legend.title = element_text(size = 16), 
            legend.text = element_text(size = 14)) +
      scale_fill_gradientn(colours= rev(brewer.pal(name = "Spectral", n = 11)),
                           breaks=c(0,0.5,1),
                           limits=c(0,max)) +
      guides(fill=guide_colourbar(barwidth=20, barheight = 2))
    
    NBH_plot + RML6_plot
  } 
  
  p_list <- mapply(srat = seurat_merged_timepoints, 
              images = simages, 
              FUN = f, 
              SIMPLIFY = FALSE, 
              USE.NAMES = TRUE)
  
  p <- wrap_plots(p_list[[1]], p_list[[2]], p_list[[3]], ncol = 1)
  
  }
  
  p
  
})




# # DIFFERENTIAL EXPRESSION RESULTS FOR BRAIN REGIONS TAB


# populate select input with brain regions that can be compared between 2 slices

 output$selectBetweenSlice <-renderUI({

 tp <- DEGs_conditions %>% filter(timepoint == input$select_BS_DEA_timepoint)
 regions <- gsub("(.+)_RML6_vs_NBH", "\\1", unique(tp$comparison))

 selectInput("BetweenSliceRegion",
              "Select Region",
              choices = regions,
              multiple = FALSE,
              width = "200px")

})

comparison1 <- reactive({
  str_before_first(input$selectComparison, pattern = " vs ")
})

comparison2 <- reactive({
  str_after_first(input$selectComparison, pattern = " vs ")
})


# plot with brain regions colored in each strain in comparison

spatial_merged_reactive <- reactive({
   
   spatial_merged <- seurat_merged_timepoints[[input$select_BS_DEA_timepoint]]
   spatial_merged@meta.data <- spatial_merged@meta.data %>% dplyr::mutate(Selected_region = case_when(Brain.Regions == input$BetweenSliceRegion ~ input$BetweenSliceRegion,
                                                                                                      TRUE ~ "other"))
   spatial_merged
 })


simages_reactive <- reactive({
  
  simages[[input$select_BS_DEA_timepoint]]
  
})

 
output$comparison1_plot <- renderPlot({
  
   SpatialDimPlot(spatial_merged_reactive(), images = simages_reactive()[1], group.by = "Selected_region") + theme(legend.position = "none") +
    scale_fill_manual(values = c("grey", "#F35D6D"), breaks = c("other", input$BetweenSliceRegion))
})


output$comparison2_plot <- renderPlot({
  SpatialDimPlot(spatial_merged_reactive(), images = simages_reactive()[1], group.by = "Selected_region") + theme(legend.position = "none") +
    scale_fill_manual(values= c("grey", "#F35D6D"), breaks = c("other", input$BetweenSliceRegion))

})


# filter DE results based on the region and comparison selected
output$DE_table <- renderDataTable({

  DEGs_conditions %>% filter(timepoint == input$select_BS_DEA_timepoint,
                             str_detect(comparison, input$BetweenSliceRegion))

})


output$GO_results_between <- renderDataTable({

  ORA_results_between_strains_dfs_merged %>% filter(timepoint == input$select_BS_DEA_timepoint,
                                                    str_detect(comparison, input$BetweenSliceRegion))

})




# BRAIN REGIONS WITHIN SLICE DE TAB

# populate select input with brain regions that exist in selected strain

DEGs_regions_reactive <- reactive({
  DEGs_regions %>% filter(timepoint == input$select_WS_DEA_timepoint, condition == input$selectStrain)
})

before_vs  <- reactive({
  x <- str_before_first(unique(DEGs_regions_reactive()$comparison), pattern = " vs")
 unique(x)
})



output$selectRegionOneSlice <-renderUI({


  selectInput("RegionOneSlice",
              "Compare",
              choices = before_vs(),
              multiple = FALSE,
              width = "200px")

})

# display all possible comparisons for selected region

output$selectComparisonRegion <- renderUI({

  selected_region <- input$RegionOneSlice

  comparisonsW_selectedRegion <- unique(DEGs_regions_reactive()$comparison)[grepl(paste(selected_region, " vs", sep = ""),
                                                                                                unique(DEGs_regions_reactive()$comparison))]

  after_vs <- str_after_first(comparisonsW_selectedRegion, pattern = "vs ")
  after_vs <- unique(after_vs)

selectInput("RegionOneSlice2",
            "To",
            choices = after_vs,
            multiple = FALSE,
            width = "200px")

})


# get the correct table from DE results based on timepoint, strain and comparison selected

output$withinSlice_DEgenes <- renderDataTable({

  comp <- paste0(input$RegionOneSlice, " vs ", input$RegionOneSlice2)
  DEGs_regions_reactive() %>% filter(comparison == comp)

})

# table with within slice DE results

output$GO_results_within <- renderDataTable({
  
  comp <- paste0(gsub(" ", "_", input$RegionOneSlice), "_vs_", gsub(" ", "_", input$RegionOneSlice2))

  ORA_results_within_strains_dfs_merged %>% filter(timepoint == input$select_WS_DEA_timepoint, 
                                                   strain == input$selectStrain, 
                                                   comparison == comp)
  
  
})




# MODULE SCORE TAB

output$moduleScores_plot <- renderPlot({

  feature <- sapply(all_spatial_seurat_normalized, 
                    FUN = function(srat){
    srat@meta.data[, input$selectGO]
                    }, 
    simplify = FALSE)
  
  feature <- unlist(feature, use.names = FALSE)
  
  minF <- min(feature)
  maxF <- max(feature)
  

  if(input$selectGO == "neuroinflammatory_response"){
    title <- "Neuroinflammation score"
    breaks <- c(-0.05, 0.1)
  }else{
    title <- "Neuron Apoptotic Process score"
    breaks <- c(0, 0.05, 0.1)
  }


  p_list <- mapply(x = all_spatial_seurat_normalized,
                   name = names(all_spatial_seurat_normalized), 
                                FUN = function(x, name){

    SpatialFeaturePlot(x, features = input$selectGO) + labs(fill = title) + ggtitle(name) +
    scale_fill_gradientn(colours= rev(brewer.pal(name = "Spectral", n = 11)),
                         breaks= breaks,
                         limits= c(minF, maxF)) +
    guides(fill=guide_colourbar(barwidth=20, barheight = 2)) +
                                    theme(plot.title = element_text(size = 18), 
                                          legend.title = element_text(size = 16), 
                                          legend.text = element_text(size = 14))
                                  


  },
  SIMPLIFY  = FALSE)

  p <- wrap_plots(p_list[[3]], p_list[[5]], p_list[[1]], p_list[[2]], p_list[[4]], p_list[[6]], ncol = 2)

  p
})



# SPATIALLY VARIABLE FEATURES TAB

output$SVF_table <- renderDataTable({
  
  all_SVFs_merged %>% filter(condition == input$SVF_condition, timepoint == input$SVF_timepoint) 
  
})

# display top 10 spatially variable features

output$SVF_plot <- renderPlot({
  
  top_10_svf <- head(all_SVFs_merged %>%
                       filter(condition == input$SVF_condition, timepoint == input$SVF_timepoint) %>% 
                       dplyr::select("gene") %>% 
                       flatten_chr(),
                       n = 10)
  
  name <- paste(input$SVF_condition, input$SVF_timepoint, sep = "_" )

  SpatialFeaturePlot(all_spatial_seurat_normalized[[name]], features = top_10_svf, ncol = 3, alpha = c(0.1, 1))

  
})





} 
  
 shinyApp(ui = ui, server = server)
 
 