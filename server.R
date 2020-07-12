library(shiny)
library(grid)
library(gridExtra)
library(dplyr)
memory.size(1000000)
library(tidyr)
library(gginnards)
library(psych)
library(ggpmisc)
library(readxl)

avg_info <- read.csv('AverageExpression_DonorInformation.csv')
rownames(avg_info) <- avg_info$X
avg_info <- subset(avg_info, select = -c(X))

dataframe <-  read.csv('DataFrame_incl_Genes.csv')
rownames(dataframe) <- dataframe$X
dataframe <- subset(dataframe, select = -c(X))

table1 <- function(){
  dataframe$age[dataframe$age == '1 mo'] <- '1'
  dataframe$age <- as.numeric(dataframe$age)
  dataframe$bmi <- as.numeric(dataframe$bmi)
  
  data <- subset(dataframe, select = c('orig.ident', 'cell_type1', 'donor', 'diabetes_type', 'age', 
                                       'bmi', 'sex', 'nCount_RNA', 'nFeature_RNA', 'method'))
  colnames(data) <- c('Dataset', 'Cell type', 'Donor', 'Diabetes type', 'Age', 'BMI', 
                      'Gender', 'Number of genes', 'Number of UMIs', 'Technology')
  return(data)
}
table2 <- function(){
  amount <- table(dataframe$cell_type1, dataframe$orig.ident)
  dataset <- as.vector(table(dataframe$orig.ident))
  celltype <- as.vector(table(dataframe$cell_type1))
  celltype<- c(sum(as.vector(table(dataframe$cell_type1))), celltype)
  amount <- rbind(dataset, amount)
  rownames(amount)[1] <- 'Total Number'
  amount <- cbind(celltype, amount)
  colnames(amount)[1] <- 'Total Number'
  return(amount)
}
data <- table1()
data2 <- table2()

shinyServer(function(input, output, session){
  ##### Clustering Page #####
  #Embedding type:
  coordinates <- reactiveValues(x = 'tSNE_1', y = 'tSNE_2')
  observeEvent(input$plt_type, {
    if(input$plt_type == 't-SNE'){
      coordinates$x <- 'tSNE_1'
      coordinates$y <- 'tSNE_2'
    }
    else{
      coordinates$x <- 'UMAP_1'
      coordinates$y <- 'UMAP_2'
    }
  })

  #Datasets checkboxes:
  observeEvent(input$check_all,{
    if(length(input$check_sets) == 7 & input$check_all == TRUE){
      updatePrettyCheckboxGroup(session, 'check_sets', selected = levels(unique(dataframe$orig.ident)))
    }
    else if(length(input$check_sets) == 7 & input$check_all == FALSE){
      updatePrettyCheckboxGroup(session, 'check_sets', selected = levels(unique(dataframe$orig.ident))[[1]])
    }
    else if(length(input$check_sets) != 7 & input$check_all == TRUE){
      updatePrettyCheckboxGroup(session, 'check_sets', selected = levels(unique(dataframe$orig.ident)))
    }
  }) 
  observeEvent(input$check_sets, {
    if(length(input$check_sets) == 7){
      updatePrettyCheckbox(session, 'check_all', value = TRUE)
    }
    else if(length(input$check_sets) < 1){
      print('Error: Too little datasets chosen')
    }
    else{
      updatePrettyCheckbox(session, 'check_all', value = FALSE)
    }
  })

  #Creating ggplots:
  show <- reactiveValues()
  observeEvent(c(input$check_sets, input$col_by, input$plt_type),{
    unselected <- as.list(setdiff(levels(unique(dataframe$orig.ident)), input$check_sets))
    dataframe$selected <- 1
    dataframe$selected[dataframe$orig.ident %in% unselected] <- 0
    show$plot1 <- ggplot(dataframe[dataframe$selected == 1,], aes_string(coordinates$x, coordinates$y, label = 'cell_type1')) +
      geom_point(size = 1.2, alpha = 0.2, aes_string(color = input$col_by)) +
      theme_bw(base_size = 14) + xlab(coordinates$x)+ylab(coordinates$y) +
      guides(color = guide_legend(ncol = 4, override.aes = list(alpha = 1,size = 5))) +
      theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),
          axis.ticks = element_blank(), panel.border = element_blank(), legend.text = element_text(size=12),
          legend.title = element_blank(),legend.background = element_blank(),
          panel.grid = element_blank(), legend.position = 'bottom')
  })
  observeEvent(c(input$check_sets, input$plt_type),{
    unselected <- as.list(setdiff(levels(unique(dataframe$orig.ident)), input$check_sets))
    dataframe$selected <- 1
    dataframe$selected[dataframe$orig.ident %in% unselected] <- 0
    show$plot2 <- ggplot(dataframe[dataframe$selected == 1,], aes_string(coordinates$x, coordinates$y)) +
      geom_point(aes(color = get(input$genes)), size = 1) +
      scale_colour_gradientn(colours=c(rev(terrain.colors(100))))+
      theme_bw(base_size = 14)+ xlab(coordinates$x) + ylab(coordinates$y)+
      theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),
      panel.border=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
      theme(legend.title=element_blank(), legend.position = 'bottom')
    })
  
  output$hover_info <- renderPrint({
    if(! is.null(input$plot1_hover)){
      point <- nearPoints(dataframe, input$plot1_hover, xvar = coordinates$x, yvar = coordinates$y, maxpoints = 1, threshold = 5)
      if(nrow(point) > 0){
        cat('Hover Information:\n')
        cat('Dataset: ', point[[1]], '\n')
        cat('Cell type: ', point[[4]], '\n')
        cat('Donor:', point[[5]], '\n')
        cat('Gender:', point[[6]], '\n')
        cat('Age:', point[[7]], '\n')
        cat('BMI:', point[[8]], '\n')
        cat('Diabetes type:', point[[9]], '\n')
        cat('Technology:', point[[10]])
      }
    }
  })
  output$hover2_info <- renderPrint({
    if(! is.null(input$plot2_hover)){
      point <- nearPoints(dataframe, input$plot2_hover, xvar = coordinates$x, yvar = coordinates$y, maxpoints = 1, threshold = 5)
      if(nrow(point) > 0){
        cat('Hover Information:\n')
        cat(input$genes, "expression: ", point[[input$genes]])
      }
    }
  })

  #Embeddings plot 1:
  output$plot1 <- renderPlot ({
    show$plot1
    }) 
  #Embeddeings plot 2:
  output$plot2 <- renderPlot ({
    if(input$genes != ""){
      show$plot2
    }
  })
  #Download button:
  output$download <- downloadHandler(
    filename = function(){
      paste(input$plt_type,'_embedded_', input$col_by,'_', input$genes, '.pdf', sep = '')},
    content = function(file){
      pdf(file)
      arrangeGrob(print(show$plot1),
                  print(show$plot2), ncol = 2)
      dev.off()
      })
  ##### DotPlot Page #####
  observeEvent(c(input$col_by_dot, input$genes_dot), {
    if(length(input$genes_dot) > 0){
      yaxis <- ifelse(input$col_by_dot == 'cell_type1', 'Cell Type', ifelse(input$col_by_dot == 'orig.ident', 'Datasets', ifelse(input$col_by_dot == 'diabetes_type', 'Diabetes Type', 'Technology')))
      df <- dataframe[, c(input$col_by_dot, input$genes_dot)] %>% gather(key = Gene, value = LogNorm, 2:length(dataframe[, c(input$col_by_dot, input$genes_dot)]))
      show$dotplot <- ggplot(df, aes_string(x='Gene', y = input$col_by_dot)) +
        geom_point(aes(color = LogNorm), size = 9) +
        scale_color_gradient(low="blue", high="red") +
        ylab(yaxis) + xlab('Selected Genes') + 
        theme_minimal() +
        theme(text = element_text(size = 20))
      }
    })
  output$dotplot <- renderPlot({
      show$dotplot
  })
  output$download_dot <- downloadHandler(
    filename = function(){
      paste('Dotplot', '.pdf', sep = '')},
    content = function(file){
      pdf(file)
      print(show$dotplot)
      dev.off()
    })
  ##### ScatterPlot Page #####
  observeEvent(c(input$col_by_sc, input$genes_sc1, input$genes_sc2),{
    if(input$col_by_sc != "" & input$genes_sc1 != "" & input$genes_sc2 != ""){
      show$scatter <- ggplot(dataframe, aes_string(x = input$genes_sc1, y = input$genes_sc2)) +
        geom_point(size = 1.5, alpha = 0.3, aes_string(color = input$col_by_sc)) +
        ylab(input$genes_sc2) + xlab(input$genes_sc1) + 
        theme_minimal() +
        theme(panel.border=element_blank(), legend.text=element_text(size=18),legend.title=element_blank(),legend.background=element_blank()) +
        guides(color=guide_legend(ncol=1, override.aes = list(alpha = 1,size=3)))
    }
  })
  output$scatter <- renderPlot({
    show$scatter
  })
  output$download_sc <- downloadHandler(
    filename = function(){
      paste('ScatterPlot_of_',input$genes_sc1, '_and_', input$genes_sc2, '.pdf', sep = '')},
    content = function(file){
      pdf(file)
      print(show$scatter)
      dev.off()
    })
  
  
  ##### Data Page #####
  output$meta_table <- DT::renderDataTable(
    DT::datatable(data, extensions = 'Scroller', filter = 'top',
                  options = list(lengthMenu = list(c(20, 50, 100), c('20', '50', '100')),
                                 pageLength = 20))
  )
  output$cells_table <- DT::renderDataTable(
    DT::datatable(data2, options = list(paging = FALSE,
                                         searching = FALSE))
  )
  ##### Regression Page #####
  observeEvent(c(input$col_by_corr, input$genes_corr, input$feature_corr),{
    if(input$col_by_corr != "" & input$genes_corr != "" & input$feature_corr != ""){
      my.formula <- y ~ x
      show$corr <- ggplot(avg_info, aes_string(x = input$feature_corr, y = input$genes_corr)) +
                 geom_point(aes_string(color = input$col_by_corr), size = 3) +
                 ylab(paste(input$genes_corr, 'expression', sep = ' ')) + xlab(input$feature_corr) +
                 theme_minimal() +
                 geom_smooth(method = 'lm', se = FALSE, col = 'black', formula = my.formula) +
                 theme(legend.text = element_text(size = 13),
                       legend.title = element_blank(),
                       axis.text = element_text(size = 13),
                       axis.title=element_text(size=15)) +
                 stat_poly_eq(formula = my.formula, aes(label = paste(..rr.label.., ..p.value.label.., sep = '~~~')), parse = TRUE)
    }
  })
  
  output$corr <- renderPlot({
    show$corr
  })
  output$download_corr <- downloadHandler(
    filename = function(){
      paste('Correlation_of_',input$genes_corr, '_expression_and_', input$feature_corr, '.pdf', sep = '')},
    content = function(file){
      pdf(file)
      print(show$corr)
      dev.off()
    })
})