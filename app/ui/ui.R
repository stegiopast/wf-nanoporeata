#### ui.R

# UI
metadata_table_box <- function(input, output, session) {
  box(
    id = "input_data",
    title = "Input Data",
    status = "primary",
    solidHeader = T,
    collapsible = T,
    width = 12,
    collapsed = F,
    fluidRow(
      column(12, DTOutput("metadata.out") %>% withSpinner(color = "#0dc5c1"))
    )
  )
}
design_matrix_box <- function(input, output, session, coll = T) {
  box(
    id = "design_matrix",
    title = "Design matrix",
    status = "primary",
    solidHeader = TRUE,
    collapsible = TRUE,
    width = 12,
    collapsed = coll,
    fluidRow(
      column(4, uiOutput("design_column_in")),
      column(4, uiOutput("feature_column_A")),
      column(4, uiOutput("feature_column_B"))
    ),
    fluidRow(
      column(4, offset = 4, uiOutput("color_feature_A")), 
      column(4, uiOutput("color_feature_B"))
    )
  )
}

metadata_table_box2 <- function(input, output, session) {
  box(
    id = "input_data2",
    title = "Input Data",
    status = "primary",
    solidHeader = T,
    collapsible = T,
    width = 12,
    height = 12,
    collapsed = F,
    fluidRow(
      column(12, DTOutput("metadata.out2") %>% withSpinner(color = "#0dc5c1"))
    )
  )
}

# Figure Descriptions 

fig_des <- list("read_length" = "The distribution of read lengths derived from generated fastq files is plotted per sample (left) and per condition (right). All reads of length over the 99 % quantile of all lengths are removed from this visualizations. This plot updates automatically when newly generated data is processed. ",
             "inner_var_samp" = "The change in gene composition (CGC) plot shows the sum of absolute differences of relative gene (%gi) counts per individual gene (g) in iteration (i) to the relative gene counts in iteration (i-1) for each respective sample. \nCGC = sum(|(%g(i-1))-(%g(i))|). ",
             "inner_var_cond" = "The change in gene composition (CGC) plot shows the mean sum of absolute differences of relative gene (%gi) counts per individual gene (g) in iteration (i) to the relative gene counts in iteration (i-1) per sample in each condition (c). ",
             "gene_counts_var_samp" = "The number of identified genes (> 0 reads counted) is plotted after each iteration. Each line corresponds to one sample. When no more additional genes are detected, the lines reach a plateau.",
             "gene_counts_var_cond" = "The sum of the number of identified genes in each condition (> 0 reads counted) is plotted after each iteration. When no more additional genes are detected, the lines reach a plateau. ",
             "process_time" = "The bar plot shows the time in seconds all preprocessing steps needed per iteration. One iteration processes at maximum 30 files from all samples. This plot updates automatically when a new process has finished. ",
             "gene_counts" = "The raw and normalized gene counts from FeatureCounts [Liao et al. 2014] are visualized for selected genes per condition using boxplots. The median-of-ratios normalization method from DESeq2 [Love et al. 2014] was used for normalization. ",
             "gene_body_coverage" = "The raw and normalized gene counts from FeatureCounts [Liao et al. 2014] are visualized for selected genes per condition using boxplots. The median-of-ratios normalization method from DESeq2 [Love et al. 2014] was used for normalization. ",
             "volcano_genes" = "The differential gene expression analysis results from DESeq2 are shown by plotting the log2FoldChange between the conditions of interest against the –log10 adjusted p-value per gene observed from the Wald-Test integrated in DESeq2 [Love, Huber, and Anders et al. 2014]. The top 10 significant genes are labeled by their symbol. ",
             "pca_genes" = "From the Principal Component Analysis of the top 500 genes with the highest variance across all samples, the two Principal Components – PC1 and PC2 - that explain the most variance in the dataset are plotted against each other. Each dot corresponds to one sample and is coloured by the respective condition. ",
             "sample2sample_genes" = "The euclidean distance between the gene expression patterns of all samples to each other is plotted using heatmap. Normalized gene counts were used for distance computation.  ",
             "heatmap_genes" = "The expression of the top 20 differentially expressed genes is plotted using Heatmap by coloring the number of reads per gene. The legend on the right side describes which colors correspond to highly and lowly expressed genes.  ",
             "volcano_transcripts" = "The differential transcript expression analysis results from DESeq2 are visualized by plotting the log2FoldChange between the conditions of interest against the –log10 adjusted p-value per transcript observed from the Wald-Test integrated in DESeq2 [Love, Huber, and Anders et al. 2014]. The top 10 significant transcripts are labeled by their symbol. ",
             "pca_transcripts" = "From the Principal Component Analysis of the top 500 transcripts with the highest variance across all samples, the two Principal Components – PC1 and PC2 - that explain the most variance in the dataset are plotted against each other. Each dot corresponds to one sample and is coloured by the respective condition. ",
             "sample2sample_transcrips" = "The euclidean distance between the transcript expression patterns of all samples to each other is plotted using heatmap. Normalized transcript counts were used for distance computation.  ",
             "heatmap_transcripts" = "The expression of the top 20 differentially expressed transcripts is plotted using Heatmap by coloring the number of reads per transcript. The legend on the right side describes which colors correspond to highly and lowly expressed transcripts.  ",
             "volcano_deu" = "The differential transcript usage analysis results from DEXSeq [Anders, Reyes, and Huber 2012] are visualized by plotting the log2FoldChange between the conditions of interest against the –log10 adjusted p-value per transcript. The top 10 significantly differentially used transcripts are labeled by their symbol. ",
             "single_gene_deu" ="The proportional expression of each transcript from the selected gene is shown via boxplots for each condition. Differentially used transcripts between the conditions derived from DRIMSeq [Nowicka and Robinson 2016] are highlighted by a red significance line above the respective boxplots.  ")  


# _______________________________________________________________________________
# UI ####

ui <- dashboardPage(
  dashboardHeader(title = "NanopoReaTA"),
  # ____________________________________________________________________________
  ## SIDEBAR ####
  dashboardSidebar(
    width = 270,
    sidebarMenu(
      id = "tabs",
      ### 1. Start page ####
      #menuItem("Welcome", tabName = "start_preprocessing_page", icon = icon("hand-spock")),
      ### 2. Preprocessing ####
      menuItem("Preprocessing", tabName = "pre_settings", icon = icon("cog", lib = "glyphicon")),
      ### 3. Run Overview ####
      menuItem("Run Overview", tabName = "run_overview", icon = icon("compass")),
      ### 4. Gene-wise Analysis ####
      menuItem("Gene-wise Analysis", tabName = "qualitycontrol", icon = icon("balance-scale-left")),
      ### 5. Differential Expression Analysis ####
      menuItem("Differential Expression Analysis", tabName = "dea_main", icon = icon("chart-bar"))
    ),
    column(12,align = "center",offset = 0,uiOutput("stop_preprocessing")),
    column(12,align = "center",offset = 0,uiOutput("resume_preprocessing"))
  ),
  ## DASHBOARD ####
  dashboardBody(
    ### THEME
    shinyDashboardThemes(
      theme = "grey_dark"
    ),
    ### OVERFLOW
    tags$head(tags$style(
      HTML(".wrapper {height: auto !important; position:relative; overflow-x:hidden; overflow-y:hidden}")
    )),

    ### Fonts size for labels
    tags$head(
      tags$style(HTML(
        "label { font-size:120%;}"
      ))
    ),
    ### Fonts size for box headers
    tags$style(type = "text/css", ".box-header h3.box-title {
           font-size: 24px;
        }
        "),
    ### Fonts size for selectInput
    # different sizes of mapping stuff comes from here
    tags$style(type = "text/css", ".selectize-input { font-size: 15px} .selectize-dropdown { font-size: 15px}"),
    ### Fonts size for textInput
    tags$style(type = "text/css", "[type = 'text'] { font-size: 15px}"),

    # Change bg color in selectInput, when dropdowns are selected
    tags$style(type = "text/css", ".selectize-control.single .selectize-input.input-active {background: #46505a}"),
    # Change color in selectInput dropdowns
    tags$style(type = "text/css", ".selectize-dropdown-content .option {background-color: #46505a; color: white}"),
    tags$style(type = "text/css", ".selectize-dropdown-content .option:hover {background-color: #28323c; color: white}"),

    # BG fileInput
    tags$style(type = "text/css", ".input-group .form-control {background: #44505a}"),
    tags$style(type = "text/css", ".progress-bar {background-color: grey}"),
    tags$style(type = "text/css", ".spinner {height:30px; width: 30px}"),

    # Change color for Radio Button selection
    # https://moderncss.dev/pure-css-custom-styled-radio-buttons/
    tags$style(type = "text/css", 'input[type="radio"] {
              -webkit-appearance: none;
              -moz-appearance: none;
              appearance: none;
              background-color: rgb(200, 200, 200);
              margin: 0;
              font: inherit;
              color: currentColor;
              width: 0.79em;
              height: 0.79em;
              border: 0.1em solid currentColor;
              border-radius: 50%;
              transform: translateY(0.3em);
              display: grid;
              place-content: center;
    } input[type="radio"]::before {
              content: "";
              width: 0.5em;
              height: 0.5em;
              border-radius: 50%;
              transform: scale(0);
              transition: 120ms transform ease-in-out;
              box-shadow: inset 1em 1em var(--form-control-color);
              /*background-color: #117733;*/
              /*background-color: #689f38;*/
              background-color: #46505a;
              /*background-color: #6699CC;*/

    } input[type="radio"]:checked::before {
              transform: scale(1);
    }
               '),

    # Change color for sliderInput
    # add new ".js-irs-*" when adding new sliderInput!!!!
    tags$head(tags$style(HTML(".js-irs-0 .irs-single, .js-irs-0 .irs-bar-edge, .js-irs-0 .irs-bar {
                                                  background: #6699CC;
                                                  border-top: 1px solid #6699CC ;
                                                  border-bottom: 1px solid #6699CC ;}
                               .js-irs-1 .irs-single, .js-irs-1 .irs-bar-edge, .js-irs-1 .irs-bar {
                                                  background: #6699CC;
                                                  border-top: 1px solid #6699CC ;
                                                  border-bottom: 1px solid #6699CC ;}

                            /* changes the colour of the number tags */
                           .irs-from, .irs-to, .irs-single { background: #6699CC }"))),
    
    ## Input files
    tags$head(
      tags$style(HTML("
            code {
                display:block;
                padding:9 px;
                margin:0 0 10px;
                margin-top:10px;
                font-size:13px;
                line-height:20px;
                word-break:break-all;
                word-wrap:break-word;
                white-space:pre-wrap;
                background-color:#F5F5F5;
                border:1px solid rgba(0,0,0,0.15);
                border-radius:4px; 
                font-family:monospace;
            }"))),
    
    ### 1. Start page ####
    tabItems(
      ## 2. Preprocessing ####
      tabItem(
        tabName = "pre_settings",
        fluidPage(
          tags$head(
            tags$style(HTML("hr {border-top: 0.5px solid #5c5e5c;}"))
          ),
          tabBox(
            width = 12,
            # The id lets us use input$tabset1 on the server to find the current tab
            id = "tabset1",

            tabPanel(
              title = "Sample settings", value = "panel2", fluidRow(
                column(12, uiOutput("metadata.tables.out")),
              ),
              fluidRow(
                column(12, uiOutput("design.matrix.out"))
              )
            ),
          )
        )
       ),
      ### 3. Run Overview ####
      tabItem(
        tabName = "run_overview",
        fluidPage(
          fluidRow(
            column(12,

                withSpinner(DT::dataTableOutput("table_seq_run"), color = "#0dc5c1")
              ),
            column(
            12, 
            tabBox(
              width = 12,
              title = "",
              id = "main_quality_control",
              tabPanel(
                title = "Read length distribution", value = "read_length_dist_panel",
                fluidRow(
                  column(
                    6,
     
                    div(style = "display: inline-block; vertical-align:top", downloadButton("samplewise_read_length.down")),
                    div(style = "display: inline-block; width: 80px; vertical-align: top", selectInput("samplewise_read_length.type", NULL, choices = c("png", "pdf"))),
                    plotOutput("samplewise_read_length.out") %>% withSpinner(color = "#0dc5c1")
                  ),
                  column(
                    6,
                    div(style = "display: inline-block; vertical-align:top", downloadButton("groupwise_read_length.down")),
                    div(style = "display: inline-block; width: 80px; vertical-align: top", selectInput("groupwise_read_length.type", NULL, choices = c("png", "pdf"))),
                    plotOutput("groupwise_read_length.out") %>% withSpinner(color = "#0dc5c1")
                  ),
                  column(12, br(),
                  p(em(fig_des$read_length))       
                  )
                )
              ),
              tabPanel(
                title = "Gene expression variability", value = "gene_variablity_panel",
                fluidRow(
                  column(
                    12,
                    box(column(
                      6,
                      div(style = "display: inline-block; vertical-align:top", downloadButton("total_genes_sample.down")),
                      div(style = "display: inline-block; width: 80px; vertical-align: top", selectInput("total_genes_sample.type", NULL, choices = c("png", "pdf"))),
                      plotOutput("total_genes_sample.out")%>% withSpinner(color = "#0dc5c1")
                    ),
                    column(
                      6,
                      div(style = "display: inline-block; vertical-align:top", downloadButton("inner_var_sample.down")),
                      div(style = "display: inline-block; width: 80px; vertical-align: top", selectInput("inner_var_sample.type", NULL, choices = c("png", "pdf"))),
                      plotOutput("inner_var_sample.out")%>% withSpinner(color = "#0dc5c1")
                    ),
                    column(6, br(), p(em(fig_des$gene_counts_var_samp))),
                    column(6, br(), p(em(fig_des$inner_var_samp))),
                    title = "Sample-wise", solidHeader = T, status = "primary", width = 12, collapsible = F
                    )
                  )
                ),
                fluidRow(
                  column(
                    12,
                    box(column(
                      6,
                      div(style = "display: inline-block; vertical-align:top", downloadButton("total_genes_condition.down")),
                      div(style = "display: inline-block; width: 80px; vertical-align: top", selectInput("total_genes_condition.type", NULL, choices = c("png", "pdf"))),
                      plotOutput("total_genes_condition.out")%>% withSpinner(color = "#0dc5c1")
                    ),
                    column(
                      6,
                      div(style = "display: inline-block; vertical-align:top", downloadButton("inner_var_condition.down")),
                      div(style = "display: inline-block; width: 80px; vertical-align: top", selectInput("inner_var_condition.type", NULL, choices = c("png", "pdf"))),
                      plotOutput("inner_var_condition.out")%>% withSpinner(color = "#0dc5c1")
                    ),
                    column(6, br(), p(em(fig_des$gene_counts_var_cond))),
                    column(6, br(), p(em(fig_des$inner_var_cond))),
                    title = "Condition-wise", solidHeader = T, status = "primary", width = 12, collpsible = F
                    )
                  )
                ),
                column(12, br(), p(em(fig_des$inner_var)))
              )
            )
          ))
        )
      ),
      ### 4. Gene-wise analysis ####
      tabItem(
        tabName = "qualitycontrol",
        fluidPage(
          tabBox(
            width = 12,
            # The id lets us use input$tabset1 on the server to find the current tab
            id = "tabset3",
            tabPanel(
              title = "Gene Counts", value = "tea.tab",
              fluidRow(
                column(6, box(
                  title = "Select genes of interest",
                  width = 12,
                  status = "primary",
                  solidHeader = T,
                  collapsible = T,
                  collapsed = F,
                  column(
                    12,
                   
                    withSpinner(DT::dataTableOutput("table_of_genes_df"), color = "#0dc5c1")
                  )
                )), 
                column(6, box(
                  title = "Selected genes",
                  width = 12,
                  status = "primary",
                  solidHeader = T,
                  collapsible = T,
                  collapsed = F,
                  column(
                    12,
                    withSpinner(DT::dataTableOutput("selectedGenes"), color = "#0dc5c1")
                  ),

                ))
              ),
              fluidRow(
                useShinyjs(), column(12, uiOutput("submit_gene_selection.out"), align = "center", style = "margin-bottom: 10px; margin-top: -10px")
              ),
              fluidRow(column(12, uiOutput("counts_plots_box")))
            ),
            tabPanel(
              title = "Gene Body Coverage",
              value = "genebodycoverage.tab",
              fluidRow(
                column(6, box(
                  title = "Select one gene of interest",
                  width = 12,
                  status = "primary",
                  solidHeader = T,
                  collapsible = T,
                  collapsed = F,
                  column(
                    12,
                    
                    withSpinner(DT::dataTableOutput("table_of_genes_df_gC"), color = "#0dc5c1")
                  )
                )), 
                column(6, box(
                  title = "Selected gene",
                  width = 12,
                  status = "primary",
                  solidHeader = T,
                  collapsible = T,
                  collapsed = F,
                  column(
                    12,
                    withSpinner(DT::dataTableOutput("selectedGene_gC"), color = "#0dc5c1")
                  ),

                ))
              ),
              fluidRow(
                useShinyjs(), column(12, uiOutput("submit_gene_selection_gC.out"), align = "center", style = "margin-bottom: 10px; margin-top: -10px")
              ),
              fluidRow(column(12, uiOutput("geneBodyCoveragePlot")),
                       column(12, br(), p(em(fig_des$gene_body_coverage))))
            )
          )
        )
      ),
      ### 5. DEA ####
      tabItem(
         tabName = "dea_main",
         fluidPage(
           fluidRow(
           box(
                 title = "Differential Expression Analysis",
                 width = 12,
                 status = "primary",
                 solidHeader = T,
                 collapsible = F,
                 collapsed = F,
                column(12, h3(p(em(strong("Significance threshold")))), align = "left", style = "margin-bottom: 10px; margin-top: -10px"),
                column(12, textInput("pvalue","adjusted p-value", value = 0.05)),
                column(4, uiOutput("submit_dge"), align = "center", style = "margin-bottom: 10px;", style = "margin-top: -10px;"),
                column(4, uiOutput("submit_dte"), align = "center", style = "margin-bottom: 10px;", style = "margin-top: -10px;"),
              column(4, uiOutput("submit_dt_preprocess"), align = "center", style = "margin-bottom: 10px;", style = "margin-top: -10px;")
           ),
          column(12,
           tabBox(
             width = 12, id = "dea.box",
       #### DEA results ####
             tabPanel(
               title = "Gene expression", value = "dea.tab",
               tabItem(
                tabName = "dea",
                uiOutput("dea_results.out")
               )
             ), # eof tabPanel
             tabPanel(
               title = "Transcript expression", value = "dte.tab",
               tabItem(
                 tabName = "dte",
                 uiOutput("dte_results.out")
               )
             ),
             #### DEU results ####
             tabPanel(
               title = "Transcript usage", value = "deu.tab",
               tabItem(
                 tabName = "deu",
                 fluidRow(
                   column(12,
                      fluidRow(column(12,box( 
                                        title = "Differential transcript usage",
                                        width = 12,
                                        tabBox(
                                            width = 12,
                                            tabPanel(
                                            title = "DEXSeq analysis",
                                            width = 12,
                                            status = "primary",
                                            solidHeader = T,
                                            collapsible = T,
                                            collapsed = F,
                                            DT::dataTableOutput("dex_tab") %>% withSpinner(color = "#0dc5c1")
                                            ),
                                            tabPanel(
                                            title = "DRIMSeq analysis",
                                            width = 12,
                                            status = "primary",
                                            solidHeader = T,
                                            collapsible = T,
                                            collapsed = F,
                                            DT::dataTableOutput("drim_tab") %>% withSpinner(color = "#0dc5c1")
                                            )
                                          )))),
                      tabBox(
                        width = 12, title = "Transcript analysis options",
                        # The id lets us use input$tabset1 on the server to find the current tab
                        id = "tabset.deu",
                        tabPanel(
                          title = "General", value = "dtu_dex.tab",
                          fluidPage(
                            fluidRow(column(12, 
                              div(style = "display: inline-block; vertical-align: top", downloadButton("down_vol_dex")),
                              div(style = "display: inline-block; vertical-align: top; width: 80px", selectInput("down_vol_dex.type", NULL, choices = c("png", "pdf"))),
                              plotOutput("volcano_dex_output") %>% withSpinner(color = "#0dc5c1")
                                    ),
                              column(12, br(), p(em(fig_des$volcano_deu))))
                          )
                        ),
                        tabPanel(
                            title = "Gene specific", value = "dtu.tab",
                            fluidPage(
                              fluidRow(
                                column(6, box( 
                                  title = "Select one gene of interest",
                                  width = 12,
                                  status = "primary",
                                  solidHeader = T,
                                  collapsible = T,
                                  collapsed = F,
                                  column(
                                    12,
                                    withSpinner(DT::dataTableOutput("dtu_tab"), color = "#0dc5c1")
                                  )
                                )), 
                                column(6, box(
                                  title = "Selected gene",
                                  width = 12,
                                  status = "primary",
                                  solidHeader = T,
                                  collapsible = T,
                                  collapsed = F,
                                  column(
                                    12,
                                    withSpinner(DT::dataTableOutput("selectedGene_dtu"), color = "#0dc5c1")
                                  ),

                                ))
                              ),
                              fluidRow(
                                useShinyjs(), column(12, uiOutput("submit_dtu"), align = "center", style = "margin-bottom: 10px; margin-top: -10px")
                              ),
                              fluidRow(
                                column(12, uiOutput("dtu_boxplot_UI"))
                              )
              
                            ) # fluidPage closed
                        ) #TabPanel closed
                        ) # tabBox DEU
                      ) # column
                      ) #fluidRow
                    ) #tabItem
                  ) # tabPanel Close
           ))) # tabBox
         )
      ) # eof tabItem DEA
    ) # tabItems eof
  ) # dashboardBody
) # ui
