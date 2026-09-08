library(shiny)
library(DT)
library(ggplot2)
library(dplyr)

# Load data (plain data.frame, no DESeq2 dependency needed)
res_table <- readRDS("data/res_table_full.rds")


# ============================================================
# UI
# ============================================================
ui <- fluidPage(
  titlePanel("RNA-seq Differential Expression Explorer: Old vs AD"),
  
  sidebarLayout(
    sidebarPanel(
      sliderInput("padj_cutoff", "Adjusted p-value cutoff:",
                  min = 0.001, max = 0.2, value = 0.1, step = 0.001),
      sliderInput("lfc_cutoff", "|log2FoldChange| cutoff:",
                  min = 0, max = 10, value = 2, step = 0.5)
    ),
    
    mainPanel(
      plotOutput("volcano_plot", height = "500px"),
      br(),
      DTOutput("results_table")
    )
  )
)

# ============================================================
# SERVER
# ============================================================
server <- function(input, output, session) {
  
  # Reactive: reclassify genes based on current slider values
  classified_data <- reactive({
    res_table %>%
      mutate(Expression = case_when(
        log2FoldChange >= input$lfc_cutoff & padj <= input$padj_cutoff ~ "Up-regulated",
        log2FoldChange <= -input$lfc_cutoff & padj <= input$padj_cutoff ~ "Down-regulated",
        TRUE ~ "Unchanged"
      ))
  })
  
  output$volcano_plot <- renderPlot({
    ggplot(classified_data(), aes(log2FoldChange, -log10(padj))) +
      geom_point(aes(col = Expression)) +
      scale_color_manual(values = c("Down-regulated" = "darkturquoise",
                                    "Unchanged" = "grey25",
                                    "Up-regulated" = "darkorange")) +
      geom_vline(xintercept = c(-input$lfc_cutoff, input$lfc_cutoff), 
                 colour = "grey", linetype = "dotted") +
      geom_hline(yintercept = -log10(input$padj_cutoff), 
                 colour = "grey", linetype = "dotted") +
      theme_bw() +
      theme(legend.position = "top")
  })
  
  output$results_table <- renderDT({
    classified_data() %>%
      filter(Expression != "Unchanged") %>%
      select(gene_id, log2FoldChange, padj, baseMean, Expression) %>%
      arrange(padj)
  }, options = list(pageLength = 10))
}

# ============================================================
shinyApp(ui = ui, server = server)