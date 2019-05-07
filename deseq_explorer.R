library(shiny)
library(plotly)
library(DESeq2)
library(tidyverse)
library(heatmaply)
library(shinythemes)

#Very basic and general shiny app to explore DESeq2 results. The app expects the dds, vst and ds objects to be in the ./data folder

get_pca <- function(counts, top_n = 500, flip = F){
  row_ix <- head(order(rowVars(counts), decreasing = T), n = top_n)
  pca_matrix <- counts[row_ix,]
  if (flip){
    pca_matrix <- t(pca_matrix)
  }
  pca <- prcomp(t(pca_matrix), center = T, scale. = T)
  pca
}

scree_plot <- function(pca){
  vars <- pca$sdev**2
  rel_vars <- vars / sum(vars)
  as.data.frame(rel_vars) %>% mutate(PC = paste0("PC", 1:length(rel_vars))) %>% ggplot(aes(x = reorder(PC, -rel_vars), y = rel_vars)) + geom_col() + labs(title = "Scree plot", subtitle = "Relative variance of each PC", x = "PC", y = "quotient of variance") + theme_bw()
}

get_top_genes <- function(ds, n_top = 100){
  res <- results(ds, tidy = T)
  res <- res[!is.na(res$padj),]
  top_ix <- head(order(res$padj, decreasing = F), n_top)
  res$row[top_ix]
}

dds <- readRDS("./data/ta_dds.rds")
vst <- readRDS("./data/ta_vst.rds")
colnames(vst) <- as.character(colData(dds)$sample)
ds <- readRDS("./data/ta_ds.rds")



colData(dds)$sample <- factor(colData(dds)$sample)
colData(dds)$cond <- factor(colData(dds)$cond)


pca <- as.data.frame(get_pca(assay(vst))$x)
pca$cond <- as.character(colData(dds)$cond)


ui <- fluidPage(theme = shinytheme("united"),
  tabsetPanel(
    tabPanel(
      "Gene browser",
      fluidRow(
        column(4, 
               textInput("selected_gene", "Input gene name (MGI nomenclature)", value = "Myog")),
        column(8,
               plotlyOutput("plotcounts")
               )
      )
    ),
    tabPanel(
      "Global overview",
      fluidRow(
        column(12,
               titlePanel("Sample overview"),
               plotlyOutput("pca"),
               plotlyOutput("heatmap")
               )
        
      )
    )
  )
)

server <- function(input, output){
  
  get_gene <- reactive({
    querygene <- input$selected_gene
    querygene
    })
  
  output$plotcounts <- renderPlotly({
    counts <- plotCounts(dds, gene = get_gene(), intgroup = "cond", returnData = T)
    p <- counts %>% ggplot(aes(x = cond, y = count)) + geom_point() + labs(title = "Gene expression")
    ggplotly(p)
  })
  
  
  
  output$pca <- renderPlotly({
    plot_ly(data = pca, x = ~PC1, y = ~PC2, z = ~PC3, color = ~cond, type = "scatter3d", title = "PCA")
  })
  
 
  output$heatmap <- renderPlotly({
    ix <- get_top_genes(ds)
    heatmaply(as.matrix(assay(vst))[ix, ])
  })
  
  
  
}

shinyApp(ui, server)
