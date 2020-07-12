library(shiny);library(RColorBrewer); library(shinythemes);library(shinydashboardPlus);library(ggplot2)

library(shiny)
library(shinyWidgets)
library(shinydashboard)
library(shinycssloaders)
library(plotly)

avg_info <- read.csv('AverageExpression_DonorInformation.csv')
rownames(avg_info) <- avg_info$X
avg_info <- subset(avg_info, select = -c(X))

dataframe <-  read.csv('DataFrame_incl_Genes.csv')
rownames(dataframe) <- dataframe$X
dataframe <- subset(dataframe, select = -c(X))

gene.options <- colnames(dataframe)[16:5020]
options(spinner.color="#e95420", spinner.color.background="#ffffff", spinner.size=1)


shinyUI(
    navbarPage('Visualization of Seurat Integration - scRNAseq Pancreatic Data',
               theme = shinytheme('united'),
        tabPanel('Clustering',
            fluidPage(
                fluidRow(
                 h2('Cluster Overview'),
                 p('Data integration performed on 37,615 pancreatic cells across 64 donors and 7 studies.'),
                 hr(style = 'width: 40%; margin-left: 0;'),
                 
                 column(2, prettyRadioButtons('plt_type', h4('Reduction type:'), choices = list("t-SNE", "UMAP"),
                              inline = TRUE, status = 'danger'),)),
                 column(5, offset = 1, withSpinner(plotOutput('plot1', hover = hoverOpts('plot1_hover')), type = 3)),
                 column(5, withSpinner(plotOutput('plot2', hover = hoverOpts('plot2_hover')), type = 3))
                ),
               fluidRow(
                   column(2, verbatimTextOutput("hover_info")),
                   column(2, offset = 7, verbatimTextOutput("hover2_info")),
                   column(1, offset = 9, downloadButton('download', 'Download Plots'))
               ),
               fluidRow(
                   hr(),
                   column(width = 2, offset = 1, h5("Incorporating datasets:"),
                             prettyCheckbox(inputId = 'check_all' , "All",  value = TRUE, status = 'danger'),
                             prettyCheckboxGroup(inputId = "check_sets","",
                                                 choices = levels(unique(dataframe$orig.ident)),
                                                 selected = levels(unique(dataframe$orig.ident)),
                                                 status = 'danger', inline = FALSE)),
                   column(2, offset = 1, selectInput('col_by', h5('Colour corresponds to:'),
                                         choices = list('Cell Type' = 'cell_type1', 'Datasets' = 'orig.ident', 'Technology' = 'method', 'Diabetes Type' = 'diabetes_type', 'Gender' = 'sex', 'Number of genes' = 'nFeature_RNA', 'Number of UMIs' = 'nCount_RNA' ))),
                   column(2, offset = 1, selectizeInput('genes', h5('Gene Expression:'), width = '100%',
                                            multiple = FALSE, choices = gene.options, selected = 'INS', 
                                            options = list(maxOptions = 5004)))
                )
        ),
        tabPanel('Data',
            fluidPage(
                fluidRow(
                    h2('Data Overview'),
                    p('The integrating resulted in a comprehensive pancreatic atlas that contains 37,615 cells across 64 donors and 7 studies.'),
                    hr(style = 'width: 50%; margin-left: 0;'),
                    column(10, offset = 1, h4("Feature information:"), DT::dataTableOutput('meta_table'))
                 ),
                 fluidRow(
                     hr(),
                     column(width = 10, offset = 1, h4("Distributions of number of cells per cell type and dataset:"),
                            DT::dataTableOutput('cells_table')),
                     column(width = 1, p(""))
                 )
            )
        ),
        navbarMenu('Biological Analysis',
                   tabPanel('DotPlot',
                            fluidPage(
                                fluidRow(
                                    h2('Gene expression plots - DotPlot'),
                                    p('Choose one or multiple genes to visualize their expression against cell type, datasets, diabetes type or technology.'),
                                    hr(style = 'width: 50%; margin-left: 0;'),
                                column(7, offset = 2, withSpinner(plotOutput('dotplot'), type = 3))
                            ),
                            fluidRow(
                                column(1, p("")),
                                column(1, offset = 9, downloadButton('download_dot', 'Download Plot'))
                            ),
                            fluidRow(
                                hr(),
                                column(4, offset = 1, selectInput('col_by_dot', h5('Colour corresponds to:'),
                                                                  choices = list('Cell Type' = 'cell_type1', 'Datasets' = 'orig.ident', 'Technology' = 'method', 'Diabetes Type' = 'diabetes_type'))),
                                column(4, offset = 0, selectizeInput('genes_dot', h5('Expression of gene(s):'), width = '100%',
                                                                     multiple = TRUE, choices = gene.options, 
                                                                     options = list(maxOptions = 5004),
                                                                     selected = c('INS', 'GCG', 'SST'))),
                                column(width = 1, p(""))
                            )
                            )),
                   tabPanel('Scatter Plot',
                            fluidPage(
                                fluidRow(
                                    h2('Gene expression plots - ScatterPlot'),
                                    p('Select two genes to visualize their expression agianst eachother. The cells can be coloured by cell type, dataset, diabetes type or technology.'),
                                    hr(style = 'width: 60%; margin-left: 0;'),
                                    column(7, offset = 2, withSpinner(plotOutput('scatter'), type = 3))
                                ),
                                fluidRow(
                                    column(1, p("")),
                                    column(1, offset = 9, downloadButton('download_sc', 'Download Plot'))
                                ),
                                fluidRow(
                                    hr(),
                                    column(4, offset = 1, selectInput('col_by_sc', h5('Colour corresponds to:'),
                                                                      choices = list('Cell Type' = 'cell_type1', 'Datasets' = 'orig.ident', 'Technology' = 'method', 'Diabetes Type' = 'diabetes_type'))),
                                    column(2, offset = 0, selectizeInput('genes_sc1', h5('Choose first gene expression:'), width = '100%',
                                                                         multiple = FALSE, choices = gene.options, selected = 'INS',
                                                                         options = list(maxOptions = 5004))),
                                    column(2, offset = 0, selectizeInput('genes_sc2', h5('Choose second gene expression'), width = '100%',
                                                                         multiple = FALSE, choices = gene.options, selected = 'GCG',
                                                                         options = list(maxOptions = 5004))),
                                    column(width = 1, p(""))
                                )
                            )),
                   tabPanel('Correlation',
                            fluidPage(
                                fluidRow(
                                    h2('Correlation between genes and features - Linear regression'),
                                    p('Select a feature (BMI or age) and a gene to calculate the correlation by using linear regression. The average donor expression can by coloured by dataset, technology, or diabetes type.'),
                                    hr(style = 'width: 75%; margin-left: 0;'),
                                    column(5, offset = 3, withSpinner(plotOutput('corr'), type = 3))
                                ),
                                fluidRow(
                                    column(1, p("")),
                                    column(1, offset = 9, downloadButton('download_corr', 'Download Plot'))
                                ),
                                fluidRow(
                                    hr(),
                                    column(4, offset = 1, selectInput('col_by_corr', h5('Colour corresponds to:'),
                                                                      choices = list('Datasets' = 'orig.ident', 'Technology' = 'method', 'Diabetes Type' = 'diabetes_type'))),
                                    column(2, offset = 0, selectizeInput('genes_corr', h5('Choose a gene expression:'), width = '100%',
                                                                         multiple = FALSE, choices = gene.options, selected = 'INS',
                                                                         options = list(maxOptions = 5004))),
                                    column(2, offset = 0, selectizeInput('feature_corr', h5('Choose the feature:'), width = '100%',
                                                                         multiple = FALSE, choices = list('BMI' = 'bmi', 'Age' = 'age'), selected = 'BMI',
                                                                         options = list(maxOptions = 4))),
                                    column(width = 1, p("")),
                                    column(width = 1, p("")),
                                    column(width = 1, p(""))
                                )
                            ))
        )
    )
)
