# MilkyWay
# version: 0.3.9b
# Author: Hee Jong Kim, William Barshop

#library(org.Hs.eg.db)
library(httr)
library(shiny)
library(shinyTree)
library(DT)
library(shinydashboard)
library(RColorBrewer)
library(plotly)
library(ggrepel)
library(ggfortify)
library(GalaxyConnector)
library(ComplexHeatmap)
library(circlize)
library(reshape2)
library(plyr)
library(NbClust)
library(splitstackshape)
library(Hmisc)
library(DiagrammeR)
library(scales)
library(rhdf5)

library(shinylorikeet)
library(sequenceViewer)
library(featureViewer)

library(RCurl)
library(rhandsontable)
library(rPython)

library(openxlsx)
library(rsvg)
library(magrittr)
library(DiagrammeRsvg)
library(data.table)
library(UpSetR)
library(cowplot)

source("./func/gx_get2.R")                    # Galaxy dataset download function - wget version
source("./func/volcano_preprocess.R")         # MSStats files preprocessing for volcano plot
source("./func/selected_volcano_preprocess.R")# Selected Protein's volcano plot data processing
source("./func/complex_heatmap.R")            # Processing and Ploting complext heatmap
source("./func/arrange_vars.R")               # arrange df vars by position
source("./func/exp_diagram.R")                # Experiment Diagram Drawing Preprocessor
source("./func/feature_parser.R")             # Protein sequence coverage - "Feature Viewer" Parser
source("./func/excel_report_generator.R")     # Excel report generator for LFQ
source("./func/excel_report_generator_qual.R")     # Excel report generator for Qual
source("./func/custom_shinytree_get_selected.R") # Custom output of get_selected function from ShinyTree

########### Shiny UI
Logged = FALSE;
ui <- dashboardPage(dashboardHeader(title = "MilkyWay"), 
                    dashboardSidebar(width = 250, 
                                     uiOutput("sidebarpanel")), 
                    dashboardBody(tags$head(tags$script(src = "enrichr.js")),
                                  uiOutput("mainbodypanel")))
########### /Shiny UI

########### Shiny SERVER
server <- function(input, output, session) {
  ##### LOGIN coupled UI randering part
  USER <- reactiveValues(Logged = Logged, api_key = "")
  
  observe({ 
    if (USER$Logged == FALSE) {
      if (!is.null(input$Login)) {
        if (input$Login > 0) {
          Username <- isolate(input$userName)
          Password <- isolate(input$passwd)
          if (length(Username) > 0 & length(Password) > 0) {
            galaxy_con <- GET("http://192.168.2.102/api/authenticate/baseauth", authenticate(Username, Password))
            if (status_code(galaxy_con) == 200) {
              USER$Logged <- TRUE
              USER$api_key <- content(galaxy_con)$api_key
            }
            else{
              cat(content(galaxy_con))
            }
          }
        } 
      }
    }    
  })
  # Sidebar rendering only if logged in
  output$sidebarpanel <- renderUI({
    if(USER$Logged == TRUE){
      div(
        sidebarMenu( id='side1', # placeholder id <- doing nothing
                     menuItem("Dataset Browser",
                              tabName = "galaxy_history_browser",
                              icon = icon("th")),
                     conditionalPanel(
                       condition = "output.analysis_selector == 'lfq'",
                       sidebarMenu(id='lfq_side', # actual conditionalpanel purpose id
                                   menuItem("LFQ Analysis Viewer",
                                            tabName = "lfq_analysis_viewer",
                                            icon = icon("bullseye")))
                       
                     ),
                     conditionalPanel(
                       condition = "output.analysis_selector == 'qual'",
                       sidebarMenu(id="qual_side", 
                                   menuItem("Qualitative Analysis Viewer",
                                            tabName = "qual_analysis_viewer",
                                            icon = icon("university")))
                     ),
                     menuItem("Default Job Submitter",
                              tabName = "galaxy_job_submitter",
                              icon = icon("dashboard")
                              ),
                     menuItem("Paired DIA+DDA Job Submitter",
                              tabName = "galaxy_job_submitter_dia_dda",
                              icon = icon("dashboard")
                              ),
                     menuItem("TMT Job Submitter",
                              tabName = "galaxy_job_submitter_tmt",
                              icon = icon("dashboard")
                              )
                     
                     
        ),
        # Qual Sidebar
        conditionalPanel(
          condition = "input.qual_side == 'qual_analysis_viewer'",
          conditionalPanel(
            condition = "input.qual_analyzer_tabset == 'Experiment Overview'",
            column(12,
                   align="center",
                   helpText("Experiment Overview Parameters")),
            numericInput("overview_qvalue_cutoff_qual", "Percolator q Value Cutoff:",
                         min = 0,
                         max = 1,
                         value = 0.01,
                         step= 0.01,
                         width = '800px')
          ),conditionalPanel(
            condition = "input.qual_analyzer_tabset == 'Protein'",
            column(12,
                   align="center",
                   helpText("Peptide Level Parameters")),
            numericInput("peptide_qvalue_cutoff_qual", "Percolator q Value Cutoff:",
                         min = 0,
                         max = 1,
                         value = 0.01,
                         step= 0.01,
                         width = '800px')
          ),
          conditionalPanel(
            condition = "input.qual_analyzer_tabset == 'PSM-ID'",
            column(12,
                   align="center",
                   helpText("PSM-ID Parameters")),
            numericInput("psm_id_qvalue_cutoff_qual", "Percolator q Value Cutoff:",
                         min = 0,
                         max = 1,
                         value = 0.01,
                         step= 0.01,
                         width = '800px'),
            column(12,
                   align="center",
                   actionButton("ms2_loading_qual", "Load MS2 Data")),
            column(12,
                   align="center",
                   helpText("In order to use Lorikeet viewer, you need to download MS2 dataset first")),
            column(12,
                   align="center",
                   tags$style("#download_psm_table_qual { color: #444; }"),
                   downloadButton("download_psm_table_qual", "Download PSM table")),
            column(12,
                   align="center",
                   helpText("Download PSM table based on selected columns and percolator q-value cutoff"))
          )
        ),
        # LFQ Sidebar
        conditionalPanel(
          condition = "input.lfq_side == 'lfq_analysis_viewer'",
          conditionalPanel(
            condition = "input.analyzer_tabset == 'Experiment Overview'",
            column(12,
                   align="center",
                   helpText("Experiment Overview Parameters")),
            numericInput("overview_qvalue_cutoff_lfq", "Percolator q Value Cutoff:",
                         min = 0,
                         max = 1,
                         value = 0.01,
                         step= 0.01,
                         width = '800px')
          ),
          conditionalPanel(
            condition = "input.analyzer_tabset == 'Analysis Metrics'",
            column(12,
                   align="center",
                   helpText("Analysis Metrics Parameters")),
            numericInput("unnorm_peptide_dist_violin_plot_mprophet_qvalue", "mProphet q Value Cutoff:",
                         min = 0,
                         max = 1,
                         value = 0.01,
                         step= 0.01,
                         width = '800px'),
            column(12,
                   align="center",
                   helpText("For Peptide Level Violin Plot"))
          ),
          conditionalPanel(
            condition = "input.analyzer_tabset == 'Protein'",
            column(12,
                   align="center",
                   helpText("Peptide Level Parameters")),
            numericInput("peptide_qvalue_cutoff", "Percolator q Value Cutoff:",
                         min = 0,
                         max = 1,
                         value = 0.01,
                         step= 0.01,
                         width = '800px')
          ),
          conditionalPanel(
            condition = "input.analyzer_tabset == 'Volcano Plot'",
            column(12,
                   align="center",
                   helpText("Volcano Plot Parameters")),
            uiOutput("exp_label"),
            column(12,
                   align="center",
                   helpText("Volcano Plot and Table")),
            numericInput("volcano_plot_pvalue_cutoff", 
                         "p Value Cutoff",
                         min = 0,
                         max = 1,
                         value = 0.01,
                         step= 0.01,
                         width = '800px'),
            numericInput("log2FC_cutoff",
                         "log2FC Cutoff",
                         min = 0,
                         max = 50,
                         value = 1,
                         step= 1),
            column(12,
                   align="center",
                   helpText("Condition Plot Parameters")),
            selectInput("error_bar_value",
                        "Condition Plot Value",
                        choices = c("Confidence Interval"="CI","Standard Deviation"="SD")),
            column(12,
                   align="center",
                   helpText("Enrichr Parameters")),
            selectInput("enrichr_grouping",
                        label = "Fold Change for Enrichr",
                        choices = c("Both", "Positive", "Negative")),
            column(12,
                   align="center",
                   actionButton("enrichr_action", "Enrichr", width='80%')),
            column(12,
                   align="center",
                   helpText("pvalue, fold change, its directionality, and selected condition will be applied"))
          ),
          conditionalPanel(
            condition = "input.analyzer_tabset == 'Heatmap'",
            column(12,
                   align="center",
                   helpText("Heatmap Parameters")),
            numericInput("heatmap_pvalue_cutoff", "p Value Cutoff:",
                         min = 0,
                         max = 1,
                         value = 0.01,
                         step= 0.01,
                         width = '800px'),
            selectInput("unidirectional",
                        label = "Directionality",
                        choices = c(TRUE, FALSE)),
            selectInput("quantCentering",
                        label = "Center Intensity Values",
                        choices = c("CONTROL", "ALL", "NONE")),
            numericInput("manual_k",
                         "Number of clusters (k-means)",
                         min = -1, 
                         max = 50, 
                         value = 3),
            column(12,
                   align="center",
                   helpText("-1 will automatically calculate optimal number of clusters via NbClust"))
          ),
          conditionalPanel(
            condition = "input.analyzer_tabset == 'PSM-ID'",
            column(12,
                   align="center",
                   helpText("PSM-ID Parameters")),
            numericInput("psm_id_qvalue_cutoff", "Percolator q Value Cutoff:",
                         min = 0,
                         max = 1,
                         value = 0.01,
                         step= 0.01,
                         width = '800px'),
            column(12,
                   align="center",
                   actionButton("ms2_loading", "Load MS2 Data")),
            column(12,
                   align="center",
                   helpText("In order to use Lorikeet viewer, you need to download MS2 dataset first")),
            column(12,
                   align="center",
                   tags$style("#download_psm_table { color: #444; }"),
                   downloadButton("download_psm_table", "Download PSM table")),
            column(12,
                   align="center",
                   helpText("Download PSM table based on selected columns and percolator q-value cutoff"))
          ),
          conditionalPanel(
            condition = "input.analyzer_tabset == 'Peptide-Quant'",
            column(12,
                   align="center",
                   helpText("Peptide-Quant Parameters")),
            numericInput("psm_quant_qvalue_cutoff", "Percolator q Value Cutoff:",
                         min = 0,
                         max = 1,
                         value = 0.01,
                         step= 0.01,
                         width = '800px'),
            numericInput("mprophet_qvalue_cutoff", "mProphet q Value Cutoff for Boundaries and bar plots:",
                         min = 0,
                         max = 1,
                         value = 0.01,
                         step= 0.01,
                         width = '800px'),
            selectInput("grid_row_num",
                        label = "Choose the number of grid rows",
                        choices = c(1,2,3,4,5,6,7,8,9,10),
                        selected = 3),
            selectInput("chrom_intensity_transform",
                        label = "Choose the type of intensity transformation",
                        choices = c("Raw", "Log2"),
                        selected = "Raw")
          )
        )
      )
    }
  })
  # Main page renderring only if logged in. Otherwise, login panel
  output$mainbodypanel <- renderUI({
    if(USER$Logged == TRUE){
      tabItems(

        ####################################################################### Galaxy Upload tool server side code
        tabItem(tabName = "galaxy_job_submitter",
      fluidRow(
          box(
          #title="Galaxy Job Design and Upload tool",
          width=12,
          id="galaxy_job_submitter_box_lfq",


        conditionalPanel(
        condition="!output.fastafileReceived",
        h2("1. Choose an Experiment Name and provide required annotations:"),
        textInput("historyName","Experiment Name:"),
        wellPanel(h4("Enter full PI names below:"),
              splitLayout(        inputPanel(textInput("pifirstName","PI First Name (Full):")),#,
                        #helpText("e.g. \"James\"")),
                        inputPanel(textInput("pilastName","PI Last Name (Full):"))#,
                        #helpText("e.g. \"Wohlschlegel\""))
              )
        ),
        textInput("sampleContactName","Collaboration Contact Name:"),
        helpText("e.g. \"Hee Jong Kim\" or \"Buck Strickland\" - This is usually the person who generated the biological material.")
        ),
        conditionalPanel(
        condition = "input.historyName.length > 0 && input.pilastName.length > 0 && input.pifirstName.length > 0 && input.sampleContactName.length > 0 && !output.fastafileReceived",
        #h3("2. Upload files:"),
        h2("2. Upload protein FASTA database:"),
        helpText("Target sequences only!"),
        fileInput('fastafile','Select FASTA file',accept=c(".fasta",".FASTA"))
        ),
        conditionalPanel(
          condition = "output.fastafileReceived && !output.skylinefileReceived && !input.QAcheck",
          h2("3. Upload empty Skyline file:"),
          helpText("File should be set up for the desired analysis, modifications, and acquisition parameters."),
          fileInput('skylinefile','Select Skyline file', accept=c(".sky")),
          checkboxInput('QAcheck','Identification Analysis Only (No Intensity Analysis)',value=FALSE)

        ),
        conditionalPanel(
          condition = "(input.QAcheck || output.skylinefileReceived) && !output.diawindowfileReceived && !input.DDAcheck",
          h2("4. If DIA, upload DIAUmpire window definitions, otherwise check DDA box:"),
          helpText("Window definitions must be saved without quotes or commas, etc, in a tsv or csv file with two columns (start and end m/z values)"),
          fileInput('diawindowfile','Select DIA Window file', accept=c(".csv",".tsv",".txt")),
          checkboxInput('DDAcheck','DDA File Upload: SKIP',value=FALSE)
        ),

        conditionalPanel(
          condition = "!output.datafilesReceived && (input.DDAcheck || output.diawindowfileReceived)",
          h2("5. Upload mass spec data:"),
          helpText("mzML files should be zlib compressed and centroided"),
          fileInput('files', 'Choose raw/mzML files', accept=c('.raw','.mzML','.mzml'),multiple=TRUE)
        ),
        conditionalPanel(condition="output.datafilesReceived",
                 wellPanel(
                 h2("6. Edit the table below, and click the save button below to send\nthe experimental design to Galaxy"),
                 #h3("Save"), 
                 actionButton("save", "Save table")
                 )
              ),
                      conditionalPanel(
          condition=FALSE,
          verbatimTextOutput('uploadedFASTA'),
          verbatimTextOutput('uploadedSkyline'),
          verbatimTextOutput('uploadedDIAwindowfile'),
          verbatimTextOutput('uploaded')
          )

        )#end of the box
      ),
      fluidRow(
         title="rhandson_box_lfq",
         id="handson_box_lfq",
         width="12",
         #height="600px",
         box(
          title="Experimental Design Table",
          id="interior_handson_box_lfq",
          width=12,
          
          rHandsontableOutput("hot")
         )
      )


        ),
    
    #DIA+DDA UPLOAD SECTION
    
    tabItem(tabName = "galaxy_job_submitter_dia_dda",
      fluidRow(
        box(
          #title="Galaxy Job Design and Upload tool",
          width=12,
          id="galaxy_job_submitter_box_dia_dda",
        conditionalPanel(
        condition="!output.fastafileReceivedDIADDA",
        h2("1. Choose an Experiment Name and provide required annotations:"),
        textInput("historyNameDIADDA","Experiment Name:",value=""),
        wellPanel(h4("Enter full PI names below:"),
              splitLayout(        inputPanel(textInput("pifirstNameDIADDA","PI First Name (Full):")),#,
                        #helpText("e.g. \"James\"")),
                        inputPanel(textInput("pilastNameDIADDA","PI Last Name (Full):"))#,
                        #helpText("e.g. \"Wohlschlegel\""))
              )
        ),
        textInput("sampleContactNameDIADDA","Collaboration Contact Name:",value=""),
        helpText("e.g. \"Hee Jong Kim\" or \"Buck Strickland\" - This is usually the person who generated the biological material.")
        ),
        conditionalPanel(
        condition = "input.historyNameDIADDA.length > 0 && input.pilastNameDIADDA.length > 0 && input.pifirstNameDIADDA.length > 0 && input.sampleContactNameDIADDA.length > 0 && !output.fastafileReceivedDIADDA",
        #h3("2. Upload files:"),
        h2("2. Upload protein FASTA database:"),
        helpText("Target sequences only!"),
        fileInput('fastafileDIADDA','Select FASTA file',accept=c(".fasta",".FASTA"))
        ),
        conditionalPanel(
          condition = "output.fastafileReceivedDIADDA && !output.skylinefileReceivedDIADDA",
          h2("3. Upload empty Skyline file:"),
          helpText("File should be set up for the desired analysis, modifications, and acquisition parameters."),
          fileInput('skylinefileDIADDA','Select Skyline file', accept=c(".sky"))
        ),

        conditionalPanel(
          condition = "!output.datafilesReceivedDIADDA",
          h2("4. Upload DDA (searchable) mass spec data:"),
          helpText("mzML files should be zlib compressed and centroided"),
          fileInput('filesDIADDA', 'Choose raw/mzML files', accept=c('.raw','.mzML','.mzml'),multiple=TRUE)
        ),
        conditionalPanel(condition="output.datafilesReceivedDIADDA",
                 wellPanel(
                 h2("5. Edit the DDA file table below, and click the save button below to send\nthe experimental design to Galaxy"),
                 actionButton("saveDIADDA", "Save DDA table")
                 )
              ),
        conditionalPanel(
          condition = "output.datafilesReceivedDIADDA",
          h2("6. Upload DIA (quantitative) mass spec data:"),
          helpText("mzML files should be zlib compressed and centroided"),
          fileInput('filesDIADDA_DIA', 'Choose raw/mzML files', accept=c('.raw','.mzML','.mzml'),multiple=TRUE)
        ),
        conditionalPanel(condition="output.datafilesReceivedDIADDA",
                 wellPanel(
                 h2("7. Edit the DIA file table below, and click the save button below to send\nthe experimental design to Galaxy"),
                 actionButton("saveDIADDA_DIA", "Save DIA table")
                 )
              ),


                      conditionalPanel(
          condition=FALSE,
          verbatimTextOutput('uploadedFASTADIADDA'),
          verbatimTextOutput('uploadedSkylineDIADDA'),
          verbatimTextOutput('uploadedDIAwindowfileDIADDA'),
          verbatimTextOutput('uploadedDIADDA'),
          verbatimTextOutput('uploadedDIADDA_DIA')
          )

        )#end of the box
      ),
      fluidRow(
         title="DDA Files",
         id="handson_box_dia_dda",
         width="12",
         #height="600px",
         box(
          title="DDA Experimental Design Table",
          id="interior_handson_box_dia_dda",
          width=12,
          
          rHandsontableOutput("hotDIADDA")
           )
        ),#endofHandsonFluidRow
      fluidRow(
         title="DIA Files",
         id="handson_box_dia_dda_DIA",
         width="12",
         #height="600px",
         box(
          title="DIA Experimental Design Table",
          id="interior_handson_box_dia_dda_DIA",
          width=12,
          
          rHandsontableOutput("hotDIADDA_DIA")
           )
        )#endofHandsonFluidRow
    ),
    
    
    #TMT UPLOAD SECTION
    
    tabItem(tabName = "galaxy_job_submitter_tmt",
      fluidRow(
        box(
          #title="Galaxy Job Design and Upload tool",
          width=12,
          id="galaxy_job_submitter_box_tmt"
          ), #below this goes another rhansdsontable
      fluidRow(
         title="rhandson_box",
         id="handson_box_tmt",
         width="12",
         #height="600px",
         box(
          title="Experimental Design Table",
          id="interior_handson_box_tmt",
          width=12,
          
          rHandsontableOutput("hot_tmt")
           )
        )#endofHandsonFluidRow
          
      )
    ),





        ####################################################################### Galaxy Upload tool server side code











        ####################################################################### /Galaxy Upload tool server side code v.0.2.0
        tabItem(tabName = "galaxy_history_browser",
                fluidRow(
                  box(title = "Collaborator",
                      width = 4,
                      status = "info",
                      dataTableOutput("pi_table")),
                  box(title = "Experiments",
                      width = 8,
                      status = "info",
                      dataTableOutput("history_table"))
                ),
                fluidRow(
                  box(title = "History Progress",
                      width = 8,
                      status = "info",
                      dataTableOutput("selected_history_total_dataset")),
                  box(title = "Loadable Data Set",
                      width = 4,
                      status = "info",
                      dataTableOutput("selected_history_rdataset"))
                ),
                fluidRow(
                  column(width = 8,
                         ""),
                  column(width = 4,
                         actionButton("loading", "Load the Analysis", 
                                      icon("paper-plane"),
                                      style="color: #fff; background-color: #337ab7; border-color: #2e6da4; padding:4px; font-size:150%", 
                                      width='100%'))
                ),
                
                "Galaxy History Browser"),
        # qual_analysis_viewer content
        tabItem(tabName = "qual_analysis_viewer",
                fluidRow(
                  # hide error msg from the front end
                  tags$style(type="text/css",
                             ".shiny-output-error { visibility: hidden; }",
                             ".shiny-output-error:before { visibility: hidden; }"),
                  tabBox(
                    title = "Qualitative Analysis Viewer",
                    width = 12,
                    id = "qual_analyzer_tabset",
                    tabPanel("Experiment Overview",
                             box(title = "Experimental Design Diagram",
                                 width = 6,
                                 status = "info",
                                 grVizOutput('exp_diagram_qual', 
                                             width = "100%", 
                                             height = "600px")),
                             box(title = "Delta Mass per Run",
                                 width = 6,
                                 status = "warning",
                                 plotOutput('ppm_hist_qual', 
                                            width = "100%", 
                                            height = "300px")),
                             box(title = "PSM Count per Run",
                                 width = 6,
                                 status = "warning",
                                 plotOutput('psm_count_plot_qual',
                                            width = "100%",
                                            height = "300px")),
                             box(title = "Experimental Design Table",
                                 width = 12,
                                 status = "info",
                                 dataTableOutput('experimental_design_table_qual')),
                             tags$hr(),
                             "Experimental Data Overview"
                    ),
                    tabPanel("Workflow Parameters",
                             box(title = "Workflow Tool Tree",
                                 width = 6,
                                 status = "info",
                                 shinyTree("ParameterTree_qual")),
                             box(title = "Selected Entry's Parameters",
                                 width = 6,
                                 status = "info",
                                 htmlOutput("SelectedParameterText_qual", container = tags$p)),
                             tags$hr(),
                             "Workflow Parameters"
                    ),
                    tabPanel("Analysis Metrics",
                             box(title = "FIDO ROC Curve",
                                 width = 4,
                                 status = "warning",
                                 plotOutput("fido_roc_curve_qual")),
                             box(title = "PSM ROC Curve per Run",
                                 width = 4,
                                 status = "warning",
                                 plotOutput("psm_roc_curve_qual",
                                            height = "600px")),
                             box(title = "Protein ROC Curve per Run",
                                 width = 4,
                                 status = "warning",
                                 plotOutput("protein_roc_curve_qual",
                                            height = "600px")),
                             box(title = "PCA Plots and Tables",
                                 width = 12,
                                 status = "info",
                                 tabBox(
                                   title = "PCA",
                                   id = "pca_tabset1",
                                   width = 12,
                                   tabPanel("SpC based", plotlyOutput("spc_pca_plot_qual"), dataTableOutput("spc_pca_table_qual")),
                                   tabPanel("NSAF based", plotlyOutput("nsaf_pca_plot_qual"), dataTableOutput("nsaf_pca_table_qual")))),
                             "Analysis Metrics for Quality Control"
                    ),
                    tabPanel("Protein",
                             box(title = "Protein Table with SAINTexpress",
                                 solidHeader = TRUE,
                                 width = 12,
                                 status = "info",
                                 dataTableOutput("protein_saint_table_qual")),
                             box(title = "Protein Quantitation Table",
                                 width = 12,
                                 status = "info",
                                 tabBox(
                                   title = "Protein Quantitation Measurements",
                                   id = "protein_quant_tabset1",
                                   width = 12,
                                   tabPanel("SpC", dataTableOutput("spc_table_qual")),
                                   tabPanel("NSAF", dataTableOutput("nsaf_table_qual")))),
                             box(title = "Selected Protein Sequence",
                                 width = 5,
                                 status = "info",
                                 sequenceViewerOutput("protein_seq_qual")),
                             box(title = "Selected Protein Sequence Coverage and Feature Map",
                                 solidHeader = TRUE,
                                 width = 7,
                                 status = "info",
                                 featureViewerOutput("protein_feature_qual"),
                                 style = 'display:block;width:100%;overflow-y: scroll;'),
                             box(title = "Unique Peptide Table Column Selection",
                                 width=12,
                                 status = "warning",
                                 checkboxInput("unique_peptide_column_show_qual", "Show Column Selections"),
                                 conditionalPanel(condition = "input.unique_peptide_column_show_qual == true",
                                                  uiOutput("unique_peptide_columns_ui_qual"))),
                             box(title = "Unique Peptide Table from the Selected Protein",
                                 width = 12,
                                 sataus = "info",
                                 dataTableOutput("unique_peptide_table_qual")),
                             box(title = "Peptide SpC Quantitation Table",
                                 width = 12,
                                 status = "info",
                                 dataTableOutput("pep_spc_table_qual")),
                             "Protein Level Data Inspection"
                    ),
                    tabPanel("PSM-ID",
                             box(title = "PSM Table Column Selection",
                                 width=12,
                                 status = "warning",
                                 checkboxInput("psm_id_column_show_qual", "Show Column Selections"),
                                 conditionalPanel(condition = "input.psm_id_column_show_qual == true",
                                                  uiOutput("psm_id_columns_ui_qual"))),
                             box(title = "PSM Table",
                                 solidHeader = TRUE,
                                 width=12,
                                 status = "primary",
                                 dataTableOutput("psm_id_table_qual")),
                             box(title = "Selected Peptide-Spectrum-Match Table",
                                 solidHeader = TRUE,
                                 width=4,
                                 status = "primary",
                                 dataTableOutput("selected_psm_id_table_qual")),
                             box(title = "Identified MS2 Spectra Annotation",
                                 solidHeader = TRUE,
                                 width=8,
                                 status = "primary",
                                 shinylorikeetOutput("ms2lorikeet_qual", 
                                                     height = "700px")),
                             "PSM Identification"
                    ),
                    tabPanel("Excel Report",
                             box(title = "Report File Name",
                                 status = "primary",
                                 width = 6,
                                 textInput("excelFilename_qual", 
                                           NULL, 
                                           value="Report"),
                                 column(12,
                                        align="center",
                                        helpText(".xlsx extension will be added automatically"))
                             ),
                             box(title = "UpSet Plots for peptide and protein",
                                 status = "primary",
                                 width = 6,
                                 numericInput("excelPeptideQvalue_qual", 
                                              "Peptide q Value Cutoff:",
                                              min = 0,
                                              max = 1,
                                              value = 0.01,
                                              step= 0.01),
                                 numericInput("excelProteinQvalue_qual", 
                                              "Protein q Value Cutoff:",
                                              min = 0,
                                              max = 10,
                                              value = 0.01,
                                              step= 0.01)
                             ),
                             actionButton("generating_excel_report_qual", "Build the Excel Report"),
                             conditionalPanel(condition = "output.reportbuilt_qual",
                                              br(),
                                              downloadButton("excel_report_download_qual", "Download"))
                             
                    )
                  )
                )),
        # lfq_analysis_viewer content
        tabItem(tabName = "lfq_analysis_viewer",
                fluidRow(
                  # hide error msg from the front end
                  tags$style(type="text/css",
                             ".shiny-output-error { visibility: hidden; }",
                             ".shiny-output-error:before { visibility: hidden; }"),
                  tabBox(
                    title = "LFQ Analysis Viewer",
                    width = 12,
                    # The id lets us use input$tabset1 on the server to find the current tab
                    id = "analyzer_tabset",
                    
                    tabPanel("Experiment Overview",
                             # box(title = "Parameters",
                             #     width = 12,
                             #     status = "info",
                             #     valueBoxOutput("exp_overview_tab_percolator_qvalue_lfq")),
                             box(title = "Experimental Design Diagram",
                                 width = 6,
                                 status = "info",
                                 grVizOutput('exp_diagram_lfq', 
                                             width = "100%", 
                                             height = "600px")),
                             box(title = "Delta Mass per Run",
                                 width = 6,
                                 status = "warning",
                                 plotOutput('ppm_hist_lfq', 
                                            width = "100%", 
                                            height = "300px")),
                             box(title = "PSM Count per Run",
                                 width = 6,
                                 status = "warning",
                                 plotOutput('psm_count_plot_lfq', 
                                            width = "100%", 
                                            height = "300px")),
                             box(title = "Experimental Design Table",
                                 width = 12,
                                 status = "info",
                                 dataTableOutput('experimental_design_table_lfq')),
                             tags$hr(),
                             "Experimental Data Overview"
                    ),
                    tabPanel("Workflow Parameters",
                             box(title = "Workflow Tool Tree",
                                 width = 6,
                                 status = "info",
                                 shinyTree("ParameterTree_lfq")),
                             box(title = "Selected Entry's Parameters",
                                 width = 6,
                                 status = "info",
                                 htmlOutput("SelectedParameterText_lfq", container = tags$p)),
                             tags$hr(),
                             "Workflow Parameters"
                    ),
                    tabPanel("Analysis Metrics",
                             box(title = "mProphet Target vs Decoy Histogram",
                                 width = 6,
                                 status = "warning",
                                 plotOutput('mprophet_histogram', 
                                            width = "100%")),
                             box(title = "FIDO ROC Curve",
                                 width = 6,
                                 status = "warning",
                                 plotOutput("fido_roc_curve_lfq")),
                             box(title = "PSM ROC Curve per Run",
                                 width = 6,
                                 status = "warning",
                                 plotOutput("psm_roc_curve_lfq", 
                                            height = "600px")),
                             box(title = "Protein ROC Curve per Run",
                                 width = 6,
                                 status = "warning",
                                 plotOutput("protein_roc_curve_lfq", 
                                            height = "600px")),
                             box(title = "Intensity Violin Plot per Run",
                                 width = 12,
                                 status = "warning",
                                 tabBox(
                                   title = "Violin Plots",
                                   id = "violin_tabset1",
                                   width = 12,
                                   tabPanel("Protein Level Distributions", plotOutput("protein_intensity_violin_plot_lfq",height = "600px")),
                                   tabPanel("Unnormalized Peptide Level Distributions", plotOutput("unnorm_peptide_dist_violin_plot_lfq"),height="600px")
                                 )),
                             box(title = "PCA Plots and Tables",
                                 width = 12,
                                 status = "info",
                                 tabBox(
                                   title = "PCA",
                                   id = "pca_tabset1", 
                                   width = 12,
                                   tabPanel("SpC based", plotlyOutput("spc_pca_plot_lfq"), dataTableOutput("spc_pca_table_lfq")),
                                   tabPanel("NSAF based", plotlyOutput("nsaf_pca_plot_lfq"), dataTableOutput("nsaf_pca_table_lfq")),
                                   tabPanel("MSStats Quant based", plotlyOutput("msstats_quant_pca_plot_lfq"), dataTableOutput("msstats_quant_pca_table_lfq")))),
                             "Analysis Metrics for Quality Control"
                    ),
                    tabPanel("Protein",
                             #plotlyOutput("plotly_test"),
                             #plotlyOutput("quant_type_ratio_plot"),
                             plotlyOutput("quantratio_barplot", height = "130px"),
                             box(title = "Protein Table with SAINTexpress",
                                 solidHeader = TRUE,
                                 width = 12,
                                 status = "info",
                                 dataTableOutput("protein_saint_table")),
                             box(title = "Protein Quantitation Table",
                                 width = 12,
                                 status = "info",
                                 tabBox(
                                   title = "Protein Quantitation Measurements",
                                   id = "protein_quant_tabset1", 
                                   width = 12,
                                   tabPanel("SpC", dataTableOutput("spc_table")),
                                   tabPanel("NSAF", dataTableOutput("nsaf_table")))),
                             box(title = "Selected Protein Sequence",
                                 width = 5,
                                 status = "info", 
                                 sequenceViewerOutput("protein_seq")),
                             box(title = "Selected Protein Sequence Coverage and Feature Map",
                                 solidHeader = TRUE,
                                 width = 7,
                                 status = "info", 
                                 featureViewerOutput("protein_feature"),
                                 style = 'display:block;width:100%;overflow-y: scroll;'),
                             box(title = "Unique Peptide Table Column Selection",
                                 width=12,
                                 status = "warning",
                                 checkboxInput("unique_peptide_column_show", "Show Column Selections"),
                                 conditionalPanel(condition = "input.unique_peptide_column_show == true",
                                                  uiOutput("unique_peptide_columns_ui"))),
                             box(title = "Unique Peptide Table from the Selected Protein",
                                 width = 12,
                                 sataus = "info",
                                 dataTableOutput("unique_peptide_table")),
                             box(title = "Peptide SpC Quantitation Table",
                                 width = 12,
                                 status = "info",
                                 dataTableOutput("pep_spc_table")),
                             "Protein Level Data Inspection"
                    ),
                    tabPanel("Volcano Plot",
                             box(title = "Volcano Plot",
                                 solidHeader = TRUE,
                                 width = 6,
                                 status = "primary",
                                 plotlyOutput('volcano_plot')),
                             box(title = "Condition Plot from Volcano Plot's Dot Selection",
                                 width = 6,
                                 status = "warning",
                                 plotOutput("condition_plot_from_volcano")),
                             box(title = "Protein Table",
                                 solidHeader = TRUE,
                                 width = 12,
                                 status = "primary",
                                 dataTableOutput('comparison_table')),
                             box(title = "Selected Protein's Volcano Plot",
                                 width = 6,
                                 status = "primary",
                                 plotOutput("selected_protein_volcano_plot")),
                             box(title = "Selected Protein's Condition Plot",
                                 width = 6,
                                 status = "warning",
                                 plotOutput("condition_plot_from_table")),
                             box(title = "Outlier Inspector",
                                 solidHeader = TRUE,
                                 width = 12,
                                 status = "danger",
                                 dataTableOutput('comparison_outlier_table')),
                             box(title = "Outlier Quant Barplot",
                                 width = 12,
                                 plotOutput("selected_outlier_protein_quant_barplot")),
                             "Volcano Plot Coupled with filterable DataTable"
                    ),
                    tabPanel("Heatmap",
                             box(title = "Heatmap for Probability and Intensity",
                                 width = 12,
                                 status = "primary",
                                 plotOutput('complexHeatmapPlot', 
                                            height='800px')),
                             "Complex Heatmap: Enriched Probability and Intensity"
                    ),
                    tabPanel("PSM-ID",
                             box(title = "PSM Table Column Selection",
                                 width=12,
                                 status = "warning",
                                 checkboxInput("psm_id_column_show", "Show Column Selections"),
                                 conditionalPanel(condition = "input.psm_id_column_show == true",
                                                  uiOutput("psm_id_columns_ui"))),
                             box(title = "PSM Table",
                                 solidHeader = TRUE,
                                 width=12,
                                 status = "primary",
                                 dataTableOutput("psm_id_table")),
                             box(title = "Selected Peptide-Spectrum-Match Table",
                                 solidHeader = TRUE,
                                 width=4,
                                 status = "primary",
                                 dataTableOutput("selected_psm_id_table")),
                             box(title = "Identified MS2 Spectra Annotation",
                                 solidHeader = TRUE,
                                 width=8,
                                 status = "primary",
                                 shinylorikeetOutput("ms2lorikeet", 
                                                     height = "700px")),
                             "PSM Identification"
                    ),
                    tabPanel("Peptide-Quant",
                             box(title = "Peptide Table Column Selection",
                                 width = 12,
                                 status = "warning",
                                 checkboxInput("psm_quant_column_show", "Show Column Selections"),
                                 conditionalPanel(condition = "input.psm_quant_column_show == true",
                                                  uiOutput("psm_quant_columns_ui"))),
                             box(title = "Peptide Table",
                                 solidHeader = TRUE,
                                 width=12,
                                 status = "primary",
                                 dataTableOutput("psm_quant_table")),
                             box(title = "Selected Peptide Quantitation Bar Plot",
                                 width = 12,
                                 status = "warning",
                                 plotOutput("selected_peptide_quant_barplot")),
                             box(title = "Selected Chromatogram Plot",
                                 solidHeader = TRUE,
                                 width = 12,
                                 status = "primary",
                                 plotlyOutput('selected_chromatogram_plot', 
                                              height='800px')),
                             "Peptide Quantification"
                    ),
                    tabPanel("Excel Report",
                             box(title = "Report File Name",
                                 status = "primary",
                                 width = 3,
                                 textInput("excelFilename", 
                                           NULL, 
                                           value="Report"),
                                 column(12,
                                        align="center",
                                        helpText(".xlsx extension will be added automatically"))
                             ),
                             box(title = "Volcano Plot and Conditional Format",
                                 status = "primary",
                                 width = 3,
                                 numericInput("excelPvalue", 
                                              "p Value Cutoff:",
                                              min = 0,
                                              max = 1,
                                              value = 0.01,
                                              step= 0.01),
                                 numericInput("excelLogFC", 
                                              "log2 Fold Change Cutoff:",
                                              min = 0,
                                              max = 10,
                                              value = 1,
                                              step= 0.01)
                             ),
                             box(title = "UpSet Plots for peptide and protein",
                                 status = "primary",
                                 width = 3,
                                 numericInput("excelPeptideQvalue", 
                                              "Peptide q Value Cutoff:",
                                              min = 0,
                                              max = 1,
                                              value = 0.01,
                                              step= 0.01),
                                 numericInput("excelProteinQvalue", 
                                              "Protein q Value Cutoff:",
                                              min = 0,
                                              max = 10,
                                              value = 0.01,
                                              step= 0.01)
                             ),
                             box(title = "Heatmap",
                                 status = "primary",
                                 width = 3,
                                 numericInput("excelPvalueHeatMap", 
                                              "p Value Cutoff:",
                                              min = 0,
                                              max = 1,
                                              value = 0.05,
                                              step= 0.01),
                                 selectInput("excelUnidirectional",
                                             label = "Directionality",
                                             choices = c(TRUE, FALSE)),
                                 selectInput("excelQuantCentering",
                                             label = "Center Intensity Values",
                                             choices = c("CONTROL", "ALL", "NONE")),
                                 numericInput("excelManual_k",
                                              "Number of clusters (k-means)",
                                              min = -1, 
                                              max = 50, 
                                              value = 3),
                                 column(12,
                                        align="center",
                                        helpText("-1 will automatically calculate optimal number of clusters via NbClust"))
                             ),
                             actionButton("generating_excel_report", "Build the Excel Report"),
                             conditionalPanel(condition = "output.reportbuilt",
                                              br(),
                                              downloadButton("excel_report_download", "Download"))
                             
                    )
                  )
                )
        )
      )
    }else{ # show login panel
      wellPanel(textInput("userName", "Galaxy Username (Email Address Form)"),
                passwordInput("passwd", "Password"),
                br(), actionButton("Login", "Log in"),
                br(), helpText("Version 0.3.8 - Updated: 2017/08/21"))
    }
  })
  ##### /LOGIN page associated server side code
  
  ##### Galaxy Data Browser
  # Preprocessing for History Browser
  history_values <- reactiveValues(dcast_history = "", pi_df = "")
  observe({
    if (USER$api_key != ""){
      gx_init(USER$api_key,
              GALAXY_URL='http://192.168.2.102/',
              HISTORY_ID = "")
      history <- gx_list_histories()
      history <- cbind(history_index = rownames(history), history)
      if(sum(lapply(history$tags,length)>0)>0){
        tag_info <- melt(setNames(history$tags, seq_along(history$tags)))
        colnames(tag_info) <- c("value", "history_index")
        tag_info <- with(tag_info, 
                         cbind(history_index, 
                               colsplit(tag_info$value, 
                                        pattern = ":", 
                                        names = c('role', 'fullname'))))
        
        history <- merge(x=history, 
                         y=tag_info, 
                         by="history_index", 
                         all = TRUE)
        history$tags <- NULL
        dcast_history <- dcast(history, name + annotation + id ~ role, value.var = "fullname", na.rm = TRUE)
        #history <- subset(history, role=="cc")
        history_values$dcast_history <- dcast_history
        
        pi_df <- data.frame(unique(tag_info[tag_info$role == 'pi',]$fullname))
        colnames(pi_df) <- "Collaborator"
        history_values$pi_df <- pi_df
      }else{
        history_values$dcast_history <- data.frame(name=character(),
                                                   annotation=character(),
                                                   id=character(),
                                                   cc=character(),
                                                   pi=character())
        history_values$pi_df <- data.frame(Collaborator=character())
      }
      
    }
  })
  
  # Collaborator table
  output$pi_table <- renderDataTable({
    history_values$pi_df
  }, selection = 'single', options = list(
    initComplete = JS(
      "function(settings, json) {",
      "$(this.api().table().header()).css({'background-color': '#74b3dc', 'color': '#000000'});",
      "}"))
  )
  # Experiments table
  output$history_table <- renderDataTable({
    history_values$dcast_history[which(history_values$dcast_history$pi == as.character(history_values$pi_df[input$pi_table_row_last_clicked,])), c("name","annotation", "cc")]
  }, selection = 'single', options = list(
    initComplete = JS(
      "function(settings, json) {",
      "$(this.api().table().header()).css({'background-color': '#74b3dc', 'color': '#000000'});",
      "}"))
  )
  # Selected Experiment's history progress table
  output$selected_history_total_dataset <- renderDataTable({
    selected_history_id <- history_values$dcast_history[which(history_values$dcast_history$pi == as.character(history_values$pi_df[input$pi_table_row_last_clicked,])),][input$history_table_row_last_clicked,"id"]
    if(length(selected_history_id)!=0){
      gx_switch_history(selected_history_id)
      datasets <- gx_list_history_datasets()
      datatable(datasets[order(-as.numeric(row.names(datasets))), c("name", "extension", "id", "visible", "state", "hid")]) %>% formatStyle(
        'state',
        target = 'row',
        backgroundColor = styleEqual(c("paused","error","ok","running","queued"),c("#D3EAF6","#F9BFBD","#A4F0A3","#FFFFC3","#EBEBEB"))
      )
    }
  })
  
  # Parse only Rdata containing dataset
  RData_datasets <- reactive({
    selected_history_id <- history_values$dcast_history[which(history_values$dcast_history$pi == as.character(history_values$pi_df[input$pi_table_row_last_clicked,])),][input$history_table_row_last_clicked,"id"]
    if(length(selected_history_id)!=0){
      gx_switch_history(selected_history_id)
      datasets <- gx_list_history_datasets()
      ok_datasets <- datasets[datasets$state == "ok" & datasets$visible == TRUE & datasets$deleted == FALSE, ]
      RData_datasets <- ok_datasets[grepl("Rdata Output", ok_datasets$name), ]
      RData_datasets <- RData_datasets[order(-RData_datasets$hid),c("name", "extension", 'hid', 'id')]
    }
  })
  # Loadable Rdata set table
  output$selected_history_rdataset <- renderDataTable({
    RData_datasets()
  }, selection = 'single', options = list(
    initComplete = JS(
      "function(settings, json) {",
      "$(this.api().table().header()).css({'background-color': '#74b3dc', 'color': '#000000'});",
      "}"))
  )
  ##### /Galaxy Data Browser
  
  ##### Pre-processing
  # Event based RData retrival
  e = new.env()
  # RData download and load
  observeEvent(input$loading, {
    rdata_file <- RData_datasets()[row.names(RData_datasets()[input$selected_history_rdataset_row_last_clicked,]),]$hid
    withProgress(message = "Downloading and Loading the Data", value = 1, {
      load(gx_get2(rdata_file), envir=e)
      incProgress(5/10)
    })
    withProgress(message = "Preprocessing", value = 0, {
      incProgress(1/23)
      # ONLY LFQ PRE-PROCESSING
      if(e[["analysis_type"]] == "lfq"){
        # PSM preprocessing for mProphet histogram
        e[["targets"]] <- e[["msstats_skyline_input"]][e[["msstats_skyline_input"]]$`Protein Name` != 'Decoys', c('Protein Name', 'annotation_Score')]
        e[["targets"]]$TD <- "Target"
        incProgress(2/23)
        e[["decoys"]] <- e[["msstats_skyline_input"]][e[["msstats_skyline_input"]]$`Protein Name` == 'Decoys', c('Protein Name', 'annotation_Score')]
        e[["decoys"]]$TD <- "Decoy"
        incProgress(4/23)
        e[["mprophet_hist_df"]] <- merge(e[["targets"]], e[["decoys"]], by = c("Protein Name", "annotation_Score", "TD"), all = TRUE)
        incProgress(5/23)
        # Post-processing for chromatogram and msstats_skyline_input
        names(e[["msstats_skyline_input"]]) <- gsub('\\s+', '', names(e[["msstats_skyline_input"]]))
        e[["msstats_skyline_input"]] <- e[["msstats_skyline_input"]][!(e[["msstats_skyline_input"]]$ProteinName=="Decoys"),]
        incProgress(7/23)
        e[["merged_msstats_skyline_chromatograms"]] <- merge(e[["chromatograms"]], e[["msstats_skyline_input"]], by=c("FileName","PeptideModifiedSequence","PrecursorCharge", "ProductCharge", "FragmentIon", "IsotopeLabelType"))
        names(e[["boundaries_csv"]]) <- gsub('\\s+', '', names(e[["boundaries_csv"]]))
        e[["boundaries_csv"]] <- e[['boundaries_csv']][e[["boundaries_csv"]]$PrecursorIsDecoy=="False",]
        e[["boundaries_csv"]]$MinStartTime <- as.numeric(as.character(e[["boundaries_csv"]]$MinStartTime))
        e[["boundaries_csv"]]$MaxEndTime <- as.numeric(as.character(e[["boundaries_csv"]]$MaxEndTime))
        # unique peptide list for chromatogram
        e[["peptides_for_chromatograms"]] <- as.character(unique(e[["chromatograms"]]$PeptideModifiedSequence))
        incProgress(12/23)
        # Processing for Volcano Plot table
        logFC_colname_comparison_csv <- colnames(e[["comparison_csv"]])[grepl("FC", colnames(e[["comparison_csv"]]))]
        e[["MSstats_comparison_df"]] = merge(x=e[["comparison_csv"]],y=e[["fasta_df"]],by.x="Protein", by.y="id", all.x=TRUE)
        e[["MSstats_comparison_df"]]$Protein <- as.character(e[["MSstats_comparison_df"]]$Protein)
        e[["MSstats_comparison_df"]]$description <- as.character(e[["MSstats_comparison_df"]]$description)
      }
      ## COMMON PROCESSING FOR BOTH LFQ AND QUAL
      # Post-processing for Protein tab
      incProgress(20/23)
      
      if(exists("saint_table", envir = e)){
        e[["protein_tab_fasta_df"]] <- e[["fasta_df"]]
        rownames(e[["protein_tab_fasta_df"]]) <- e[["protein_tab_fasta_df"]]$id
        e[["protein_tab_fasta_df"]]$id <- NULL
        e[["protein_tab_fasta_df"]]$Sequence <- NULL
        rownames(e[["saint_table"]]) <- e[["saint_table"]]$Prey
        e[["saint_table"]]$Prey <- NULL
        colnames(e[["saint_table"]]) <- paste("SAINT", colnames(e[["saint_table"]]), sep=".")
        e[["protein_saint_table"]] <- merge(x=e[["protein_tab_fasta_df"]], y=e[["saint_table"]], by='row.names', all.y = TRUE)
        colnames(e[["protein_saint_table"]])[colnames(e[["protein_saint_table"]])=="Row.names"] <- "Protein"
        e[["protein_saint_table"]]$Protein <- as.character(e[["protein_saint_table"]]$Protein)
        e[["protein_saint_table"]]$description <- as.character(e[["protein_saint_table"]]$description)
      }else{
        e[["protein_saint_table"]] <- merge(x=e[["fasta_df"]], y=e[["fido_roc"]], by.x="id", by.y="protein group", all.y = TRUE)
        e[["protein_saint_table"]]$Sequence <- NULL
        colnames(e[["protein_saint_table"]])[colnames(e[["protein_saint_table"]])=="id"] <- "Protein"
        e[["protein_saint_table"]]$Protein <- as.character(e[["protein_saint_table"]]$Protein)
        e[["protein_saint_table"]]$description <- as.character(e[["protein_saint_table"]]$description)
      }
      
      e[["exp_design_for_pViz"]] <- e[["qual_experiment_design_df"]][,c("Original File Name", "Biological Condition")]
      colnames(e[["exp_design_for_pViz"]])[colnames(e[["exp_design_for_pViz"]]) == 'Original File Name'] <- "category"
      incProgress(23/23)
    })
  })
  
  # Assign variables from e[[]] environment for easy access
  ## COMMON LOADING
  job_history <- eventReactive(input$loading, {e[["job_history"]]})
  fasta_df <- eventReactive(input$loading, {e[["fasta_df"]]})
  protein_saint_table <- eventReactive(input$loading, {e[["protein_saint_table"]]})
  qual_experiment_design_df <- eventReactive(input$loading, {e[["qual_experiment_design"]]})
  exp_design_for_pViz <- eventReactive(input$loading, {e[["exp_design_for_pViz"]]})
  nsaf_df <- eventReactive(input$loading, {e[["nsaf_table"]]})
  spc_df <- eventReactive(input$loading, {e[["spc_table"]]})
  peptide_nsaf_df <- eventReactive(input$loading, {e[["peptide_nsaf_table"]]})
  peptide_spc_df <- eventReactive(input$loading, {e[["peptide_spc_table"]]})
  psm_df <- eventReactive(input$loading, {e[["psm_table"]]})
  fido_roc <- eventReactive(input$loading, {e[["fido_roc"]]})
  ms2_scans_df <- eventReactive(input$loading, {e[["ms2_scans"]]})
  
  ## Reactive value from analysis_type to apply conditionalpanel
  analysis_type <- eventReactive(input$loading, {e[["analysis_type"]]})
  output$analysis_selector <- eventReactive(input$loading, {e[["analysis_type"]]})
  outputOptions(output, 'analysis_selector', suspendWhenHidden = FALSE)
  
  ## ONLY LFQ LOADING
  quant_experiment_design_df <- eventReactive(input$loading, {e[["quant_experiment_design"]]})
  MSstats_comparison_df <- eventReactive(input$loading, if(analysis_type() == "lfq"){e[["MSstats_comparison_df"]]})
  MSstats_quant_df <- eventReactive(input$loading, if(analysis_type() == "lfq"){e[["quantification_csv"]]})
  MSstats_condition_plot_df <- eventReactive(input$loading, if(analysis_type() == "lfq"){e[["condition_plot_csv"]]})
  chromatogram_df <-eventReactive(input$loading, if(analysis_type() == "lfq"){e[["chromatograms"]]})
  msstats_skyline_input_df <-eventReactive(input$loading, if(analysis_type() == "lfq"){e[["msstats_skyline_input"]]})
  mprophet_hist_df <- eventReactive(input$loading, if(analysis_type() == "lfq"){e[["mprophet_hist_df"]]})
  merged_msstats_skyline_chromatograms_df <-eventReactive(input$loading, if(analysis_type() == "lfq"){e[["merged_msstats_skyline_chromatograms"]]})
  boundaries_df <- eventReactive(input$loading, if(analysis_type() == "lfq"){e[["boundaries_csv"]]})
  peptides_for_chromatograms <- eventReactive(input$loading, if(analysis_type() == "lfq"){e[["peptides_for_chromatograms"]]})
  
  
  # Download MS2 HDF File by triggered loading button
  ## For LFQ
  ms2_hdf <- eventReactive(input$ms2_loading, {
    if(!is.null(e[["ms2_dataset_id"]])){
      selected_history_dataset_for_ms2_hdf <- gx_list_history_datasets()
      hdf_hid <- selected_history_dataset_for_ms2_hdf[selected_history_dataset_for_ms2_hdf$id == e[["ms2_dataset_id"]],]$hid
      withProgress(message = "Downloading MS2 peaks...", value = 5, {
        # Download MS2 HDF file
        ms2_hdf <- gx_get2(hdf_hid)
      })
    }
  })
  ## For Qual
  ms2_hdf_qual <- eventReactive(input$ms2_loading_qual, {
    if(!is.null(e[["ms2_dataset_id"]])){
      selected_history_dataset_for_ms2_hdf <- gx_list_history_datasets()
      hdf_hid <- selected_history_dataset_for_ms2_hdf[selected_history_dataset_for_ms2_hdf$id == e[["ms2_dataset_id"]],]$hid
      withProgress(message = "Downloading MS2 peaks...", value = 5, {
        # Download MS2 HDF file
        ms2_hdf <- gx_get2(hdf_hid)
      })
    }
  })
  ##### /Pre-processing
  
  ##### Experiment Overview
  ### Experiment Design Diagram
  exp_diagram_process <- reactive({
    withProgress(message = "Drawing Experiment Diagram", value = 0, {
      graph <- exp_diagram(qual_experiment_design_df())
    })
  })
  # For LFQ
  output$exp_diagram_lfq <- renderGrViz({
    grViz({
      exp_diagram_process()$dot_code
    })
  })
  # For Qual
  output$exp_diagram_qual <- renderGrViz({
    grViz({
      exp_diagram_process()$dot_code
    })
  })
  ### /Experiment Design Diagram
  
  ### Delta Mass Plot
  ppm_hist_process <- reactive({
    withProgress(message = "Drawing Delta Mass Plot", value = 0, {
      exp_design <- qual_experiment_design_df()
      if(analysis_type() == "lfq"){q_value <- input$overview_qvalue_cutoff_lfq}
      if(analysis_type() == "qual"){q_value <- input$overview_qvalue_cutoff_qual}
      # Sort the table
      exp_design <- exp_design[order(-as.numeric(exp_design$`Test or Control`), exp_design$BioReplicate, exp_design$`Original File Name`),]
      exp_design$`Original File Name` <- paste(exp_design$`Original File Name`, '.mzML', sep='')
      exp_design$unique_factor <- paste('Crux', exp_design$`Crux File Integer`, sep=":")
      psm_subset <- psm_df()[psm_df()$`percolator q-value`<=q_value,c("file", "ppm")]
      psm_subset$unique_factor <- mapvalues(psm_subset$file, from = exp_design$`Original File Name`, to = exp_design$unique_factor)
      # Sort the factor order of unique_factor to have a consistent order
      psm_subset$unique_factor <- factor(psm_subset$unique_factor, as.character(exp_design$unique_factor))
      psm_subset
    })
  })
  # For LFQ
  output$ppm_hist_lfq <- renderPlot({
    # Formulating delta mass histogram
    psm_ppm_hist <- ggplot(ppm_hist_process(), aes(x=ppm)) + 
      geom_histogram(aes(fill=unique_factor), binwidth=.1) +
      facet_wrap(~unique_factor, ncol=4) +
      theme(legend.position = "none")
    psm_ppm_hist
  })
  # For Qual
  output$ppm_hist_qual <- renderPlot({
    # Formulating delta mass histogram
    psm_ppm_hist <- ggplot(ppm_hist_process(), aes(x=ppm)) + 
      geom_histogram(aes(fill=unique_factor), binwidth=.1) +
      facet_wrap(~unique_factor, ncol=4) +
      theme(legend.position = "none")
    psm_ppm_hist
  })
  ### /Delta Mass Plot
  
  ### PSM count Barplot
  psm_count_plot_process <- reactive({
    withProgress(message = "Drawing PSM Count Plot", value = 0, {
      exp_design <- qual_experiment_design_df()
      if(analysis_type() == "lfq"){q_value <- input$overview_qvalue_cutoff_lfq}
      if(analysis_type() == "qual"){q_value <- input$overview_qvalue_cutoff_qual}
      # Sort the table
      exp_design <- exp_design[order(-as.numeric(exp_design$`Test or Control`), exp_design$BioReplicate, exp_design$`Original File Name`),]
      exp_design$`Original File Name` <- paste(exp_design$`Original File Name`, '.mzML', sep='')
      exp_design$unique_factor <- paste('Crux', exp_design$`Crux File Integer`, sep=":")
      psm_subset <- psm_df()[psm_df()$`percolator q-value`<=q_value,c("file", "percolator q-value")]
      psm_subset$unique_factor <- mapvalues(psm_subset$file, from = exp_design$`Original File Name`, to = exp_design$unique_factor)
      # Sort the factor order of unique_factor to have a consistent order :: reverse the order since the bar plot is rotated
      psm_subset$unique_factor <- factor(psm_subset$unique_factor, as.character(exp_design$unique_factor))
    })
    psm_subset
  })
  # For LFQ
  output$psm_count_plot_lfq <- renderPlot({
    # Formulating psm count barplot
    psm_freq_hist <- ggplot(psm_count_plot_process(), aes(factor(unique_factor), fill=unique_factor)) + 
      geom_bar() + 
      xlab("MS Run Crux ID") +
      guides(fill=FALSE) 
    psm_freq_hist
  })
  # For Qual
  output$psm_count_plot_qual <- renderPlot({
    # Formulating psm count barplot
    psm_freq_hist <- ggplot(psm_count_plot_process(), aes(factor(unique_factor), fill=unique_factor)) + 
      geom_bar() + 
      xlab("MS Run Crux ID") +
      guides(fill=FALSE) 
    psm_freq_hist
  })
  ### /PSM count Barplot
  
  ### Experiment Design Table
  # For LFQ
  output$experimental_design_table_lfq <- renderDataTable(
    arrange.vars(qual_experiment_design_df(), c("Test or Control"=1,
                                           "Biological Condition"=2,
                                           "BioReplicate"=3,
                                           "Fractionation Group Name"=4,
                                           "Fractionation Group ID String"=5,
                                           "Original File Name"=6)),
    selection = 'single',
    rownames = FALSE,
    options = list(pageLength = 50,
                   order = list(list(0, 'desc'), list(2, 'asc'), list(5, 'asc')))
  )
  # For Qual
  output$experimental_design_table_qual <- renderDataTable(
    arrange.vars(qual_experiment_design_df(), c("Test or Control"=1,
                                           "Biological Condition"=2,
                                           "BioReplicate"=3,
                                           "Fractionation Group Name"=4,
                                           "Fractionation Group ID String"=5,
                                           "Original File Name"=6)),
    selection = 'single',
    rownames = FALSE,
    options = list(pageLength = 50,
                   order = list(list(0, 'desc'), list(2, 'asc'), list(5, 'asc')))
  )
  ### /Experiment Design Table
  
  ### Parameter Tree
  ## For LFQ
  # Parameter Tree
  output$ParameterTree_lfq <- renderTree({
    job_history()
  })
  # Selected Parameter Value
  output$SelectedParameterText_lfq <- renderText({
    selected_name <- get_selected(input$ParameterTree_lfq)
    selected_list <- job_history()
    for (i in 1:length(selected_name)){
      selected_list <- selected_list[[selected_name[i]]]
    }
    HTML(selected_list)
  })
  ## For Qual
  # Parameter Tree
  output$ParameterTree_qual <- renderTree({
    job_history()
  })
  # Selected Parameter Value
  output$SelectedParameterText_qual <- renderText({
    selected_name <- get_selected(input$ParameterTree_qual)
    selected_list <- job_history()
    for (i in 1:length(selected_name)){
      selected_list <- selected_list[[selected_name[i]]]
    }
    HTML(selected_list)
  })
  ### /Parameter Tree
  ##### /Experiment Overview
  
  ##### Analysis Metrics
  # mProphet Target vs Decoy Histogram
  output$mprophet_histogram <- renderPlot({
    ggplot(data=mprophet_hist_df(), aes(x=annotation_Score, fill=TD)) +
      geom_histogram(binwidth=.3, position="dodge")
  })
  ### FIDO ROC curve
  fido_roc_process <- reactive({
    subset_fido <- fido_roc()[fido_roc()$`q-value` <= 0.05,]
  })
  # For LFQ
  output$fido_roc_curve_lfq <- renderPlot({
    ggplot(fido_roc_process(), aes(x=`q-value`)) + 
      geom_step(aes(len=length(`q-value`), y=..y.. * len), stat="ecdf") + 
      scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
      labs(x = "Fido q value", y = "Protein Count")
  })
  # For Qual
  output$fido_roc_curve_qual <- renderPlot({
    ggplot(fido_roc_process(), aes(x=`q-value`)) + 
      geom_step(aes(len=length(`q-value`), y=..y.. * len), stat="ecdf") + 
      scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
      labs(x = "Fido q value", y = "Protein Count")
  })
  ### /FIDO ROC curve
  
  ### PSM ROC curve
  psm_roc_curve_process <- reactive({
    # http://stackoverflow.com/questions/18379933/plotting-cumulative-counts-in-ggplot2
    subset_psm <- psm_df()[psm_df()$`percolator q-value` <= 0.05, c("file", "percolator q-value")]
    subset_psm <- ddply(subset_psm, .(file), transform, len=length(`percolator q-value`))
  })
  # For LFQ
  output$psm_roc_curve_lfq <- renderPlot({
    ggplot(psm_roc_curve_process(), aes(x=`percolator.q.value`, color=file)) + 
      geom_step(aes(len=len, y=..y.. * len), stat="ecdf") + 
      scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
      labs(x = "Percolator q value", y = "PSM Count") + 
      theme(legend.position = "bottom", legend.direction="vertical")
  })
  # For Qual
  output$psm_roc_curve_qual <- renderPlot({
    ggplot(psm_roc_curve_process(), aes(x=`percolator.q.value`, color=file)) + 
      geom_step(aes(len=len, y=..y.. * len), stat="ecdf") + 
      scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
      labs(x = "Percolator q value", y = "PSM Count") + 
      theme(legend.position = "bottom", legend.direction="vertical")
  })
  ### PSM ROC curve
  
  ### Protein ROC curve
  protein_roc_curve_process <- reactive({
    splited_psm_protien_df <- psm_df()[psm_df()$`percolator q-value` <= 0.05, c("file","protein id","protein q-values")]
    splited_psm_protien_df <- cSplit(splited_psm_protien_df, c("protein id","protein q-values"), ",", "long")
    splited_psm_protien_df <- splited_psm_protien_df[splited_psm_protien_df$`protein q-values` <= 0.05,]
    splited_psm_protien_df <- splited_psm_protien_df[order(splited_psm_protien_df$file, splited_psm_protien_df$`protein id`, splited_psm_protien_df$`protein q-values`), ]
    splited_psm_protien_df <- splited_psm_protien_df[!duplicated(splited_psm_protien_df[,c("file", "protein id")]),]
    splited_psm_protien_df <- ddply(splited_psm_protien_df, .(file), transform, len=length(`protein q-values`))
  })
  # For LFQ
  output$protein_roc_curve_lfq <- renderPlot({
    ggplot(protein_roc_curve_process(), aes(x=`protein.q.values`, color=file)) +
      geom_step(aes(len=len, y=..y.. * len), stat="ecdf") +
      scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
      labs(x = "protein q-values", y = "Protein Count") +
      theme(legend.position = "bottom", legend.direction="vertical")
  })
  # For Qual
  output$protein_roc_curve_qual <- renderPlot({
    ggplot(protein_roc_curve_process(), aes(x=`protein.q.values`, color=file)) +
      geom_step(aes(len=len, y=..y.. * len), stat="ecdf") +
      scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
      labs(x = "protein q-values", y = "Protein Count") +
      theme(legend.position = "bottom", legend.direction="vertical")
  })
  ### /Protein ROC curve
  
  ### protein_intensity_violin_plot_lfq
  output$protein_intensity_violin_plot_lfq <- renderPlot({
    # Preprocessing for the plot
    df <- MSstats_quant_df()
    df[,1] <- NULL
    melt_df <- melt(df, id.vars="Protein")
    # Y-axis name parsing
    logFC_colname <- colnames(MSstats_comparison_df())[grepl("FC", colnames(MSstats_comparison_df()))]
    # Violin + boxplot
    ggplot(melt_df, aes(x=factor(variable), y = value)) +
      geom_violin(aes(fill=factor(variable))) +
      geom_boxplot(width=0.2, alpha = 0.3) +
      labs(x = "Runs", y = paste0("Protein ",substr(logFC_colname,1,nchar(logFC_colname)-2),"(intensity)")) +
      theme(legend.title=element_blank())
  })
  ### /protein_intensity_violin_plot_lfq
  ### Unnormalized Peptide Area Distributions violin_plot_lfq
  output$unnorm_peptide_dist_violin_plot_lfq <- renderPlot({
    df <- msstats_skyline_input_df()[which(msstats_skyline_input_df()$annotation_QValue <= input$unnorm_peptide_dist_violin_plot_mprophet_qvalue),c("FileName", "Area", "annotation_QValue")]
    exp_df <- quant_experiment_design_df()
    # Add filename with extension
    exp_df$FilenameWithExt <- paste(exp_df$`Original File Name`, 
                                    tools::file_ext(levels(df$FileName)[1]), 
                                    sep=".")
    # Give Factor Order by Biological condition 
    df$FileName <- factor(df$FileName,
                          levels = exp_df[order(exp_df$`Biological Condition`, exp_df$FilenameWithExt), "FilenameWithExt"],
                          ordered = TRUE)
    # Add biological condtion with Biorep # & Crux #
    exp_df$cond_rep <- paste(exp_df$`Biological Condition`,exp_df$BioReplicate, exp_df$`Crux File Integer`, sep="_")
    # Merge df and exp_df
    df <- merge(x=df, y=exp_df[,c("FilenameWithExt","cond_rep")], by.x="FileName", by.y="FilenameWithExt", all=TRUE)
    # Violin + boxplot
    ggplot(df, aes(x=factor(cond_rep), y = log(Area))) +
      geom_violin(aes(fill=factor(cond_rep))) +
      geom_boxplot(width=0.2, alpha = 0.3) +
      labs(x = "Runs", y = paste0("Peptide log(Area)")) +
      theme(legend.title=element_blank())
  })
  ### /Unnormalized Peptide Area Distributions violin_plot_lfq
  
  ### SpC PCA
  spc_pca <- reactive({
    # PCA SpC per Run
    spc_df <- spc_df()
    rownames(spc_df) <- spc_df[,1]
    spc_df[1] <- NULL
    t_spc_df <- as.data.frame(t(spc_df))
    t_spc_df_for_label <- cbind(files = rownames(t_spc_df), t_spc_df)
    t_spc_df_prcomp <- prcomp(t_spc_df)
    list(prcomp=t_spc_df_prcomp, label=t_spc_df_for_label)
  })
  ## For LFQ
  # SpC PCA plot
  output$spc_pca_plot_lfq <- renderPlotly({
    ggplotly(autoplot(spc_pca()$prcomp, data=spc_pca()$label, colour="files"))
  })
  # SpC PCA table
  output$spc_pca_table_lfq <- renderDataTable({
    spc_pca()$prcomp$rotation[order(abs(spc_pca()$prcomp$rotation[,"PC1"]),decreasing=TRUE),]
  })
  ## For Qual
  # SpC PCA plot
  output$spc_pca_plot_qual <- renderPlotly({
    ggplotly(autoplot(spc_pca()$prcomp, data=spc_pca()$label, colour="files"))
  })
  # SpC PCA table
  output$spc_pca_table_qual <- renderDataTable({
    spc_pca()$prcomp$rotation[order(abs(spc_pca()$prcomp$rotation[,"PC1"]),decreasing=TRUE),]
  })
  ### /SpC PCA
  
  ### NSAF PCA
  nsaf_pca <- reactive({
    # PCA NSAF per Run
    nsaf_df <- nsaf_df()
    rownames(nsaf_df) <- nsaf_df[,1]
    nsaf_df[1] <- NULL
    t_nsaf_df <- as.data.frame(t(nsaf_df))
    t_nsaf_df_for_label <- cbind(files = rownames(t_nsaf_df), t_nsaf_df)
    t_nsaf_df_prcomp <- prcomp(t_nsaf_df)
    list(prcomp=t_nsaf_df_prcomp, label=t_nsaf_df_for_label)
  })
  # For LFQ
  # NSAF PCA plot
  output$nsaf_pca_plot_lfq <- renderPlotly({
    ggplotly(autoplot(nsaf_pca()$prcomp, data=nsaf_pca()$label, colour="files"))
  })
  # NSAF PCA table
  output$nsaf_pca_table_lfq <- renderDataTable({
    nsaf_pca()$prcomp$rotation[order(abs(nsaf_pca()$prcomp$rotation[,"PC1"]),decreasing=TRUE),]
  })
  # For Qual
  # NSAF PCA plot
  output$nsaf_pca_plot_qual <- renderPlotly({
    ggplotly(autoplot(nsaf_pca()$prcomp, data=nsaf_pca()$label, colour="files"))
  })
  # NSAF PCA table
  output$nsaf_pca_table_qual <- renderDataTable({
    nsaf_pca()$prcomp$rotation[order(abs(nsaf_pca()$prcomp$rotation[,"PC1"]),decreasing=TRUE),]
  })
  ### /NSAF PCA
  
  ### MSStats PCA
  msstats_qunat_pca <- reactive({
    # PCA MSstats quant per condition
    quant_df <- MSstats_quant_df()
    if(colnames(quant_df[1]) == ""){
      quant_df[1] <- NULL
    }
    rownames(quant_df) <- quant_df[,1]
    quant_df[1] <- NULL
    quant_df <- quant_df[complete.cases(quant_df),]
    t_quant_df <- as.data.frame(t(quant_df))
    t_quant_df_for_label <- cbind(conditions = rownames(t_quant_df), t_quant_df)
    t_quant_df_prcomp <- prcomp(t_quant_df)
    list(prcomp=t_quant_df_prcomp, label=t_quant_df_for_label)
  })
  ## For LFQ
  # MSStats Quant PCA plot
  output$msstats_quant_pca_plot_lfq <- renderPlotly({
    ggplotly(autoplot(msstats_qunat_pca()$prcomp, data=msstats_qunat_pca()$label, colour="conditions"))
  })
  # MSStats Quant PCA table
  output$msstats_quant_pca_table_lfq <- renderDataTable({
    msstats_qunat_pca()$prcomp$rotation[order(abs(msstats_qunat_pca()$prcomp$rotation[,"PC1"]),decreasing=TRUE),]
  })
  ### /MSStats PCA
  ##### /Analysis Metrics
  
  ##### Protein
  ### SpC vs MSstats Protein list comparison barplot
  ## Handling the Quant Type Ratio processing (including associated protein lists)
  quant_type_ratio_values <- reactiveValues(documents = NULL)
  observe({
    if(analysis_type() == "lfq"){
      comparison_csv <- MSstats_comparison_df()
      spc_table <- spc_df()
      # pivot comparison_csv with ImputationPercentage column <- something that contains NO NA value
      cast_comparison <- dcast(comparison_csv, Protein ~ Label, value.var = "ImputationPercentage")
      merged_spc_comp <- merge(x=cast_comparison, y=spc_table, by.x="Protein", by.y="ProteinID", all=TRUE)
      
      # extract NA rows or !NA rows
      spc_only_df <- merged_spc_comp[rowSums(is.na(merged_spc_comp))>0,]
      msstats_df <- merged_spc_comp[rowSums(is.na(merged_spc_comp))==0,]
      
      # convert protein column to vector to be used for filtering the protein table
      quant_type_ratio_values$spc_only_protein_list <- as.vector(spc_only_df[['Protein']])
      quant_type_ratio_values$msstats_protein_list <- as.vector(msstats_df[['Protein']])
      
      # number of proteins
      quant_type_ratio_values$spc_only_count <- length(quant_type_ratio_values$spc_only_protein_list)
      quant_type_ratio_values$msstats_count <- length(quant_type_ratio_values$msstats_protein_list)
    }
    
  })
  ## SpC vs MSstats Quant Type Ratio horizontal barplot (interactive plotly)
  output$quantratio_barplot <- renderPlotly({
    y <- c('Quant Type Ratio')
    x1 <- round((quant_type_ratio_values$spc_only_count/(quant_type_ratio_values$spc_only_count+quant_type_ratio_values$msstats_count))*100)
    x2 <- round((quant_type_ratio_values$msstats_count/(quant_type_ratio_values$spc_only_count+quant_type_ratio_values$msstats_count))*100)
    
    data <- data.frame(y, x1, x2)

    top_labels <- c('SpC Only', 'MSstats')
    
    p <- plot_ly(data, x = ~x1, y = ~y, type = 'bar', orientation = 'h', height = 120,
                 marker = list(color = 'rgba(38, 24, 74, 0.8)',
                               line = list(color = 'rgb(248, 248, 249)', width = 1))) %>%
      add_trace(x = ~x2, marker = list(color = 'rgba(71, 58, 131, 0.8)')) %>%
      layout(xaxis = list(title = "",
                          showgrid = FALSE,
                          showline = FALSE,
                          showticklabels = FALSE,
                          zeroline = FALSE,
                          domain = c(0.15, 1)),
             yaxis = list(title = "",
                          showgrid = FALSE,
                          showline = FALSE,
                          showticklabels = FALSE,
                          zeroline = FALSE),
             barmode = 'stack',
             paper_bgcolor = 'rgb(248, 248, 255)', plot_bgcolor = 'rgb(248, 248, 255)',
             #margin = list(l = 120, r = 10, t = 140, b = 80),
             showlegend = FALSE) %>%
      # labeling the y-axis
      add_annotations(xref = 'paper', yref = 'y', x = 0.14, y = y,
                      xanchor = 'right',
                      text = y,
                      font = list(family = 'Arial', size = 12,
                                  color = 'rgb(67, 67, 67)'),
                      showarrow = FALSE, align = 'right') %>%
      # labeling the percentages of each bar (x_axis)
      add_annotations(xref = 'x', yref = 'y',
                      x = x1 / 2, y = y,
                      text = paste(data[,"x1"], '%'),
                      font = list(family = 'Arial', size = 12,
                                  color = 'rgb(248, 248, 255)'),
                      showarrow = FALSE) %>%
      add_annotations(xref = 'x', yref = 'y',
                      x = x1 + x2 / 2, y = y,
                      text = paste(data[,"x2"], '%'),
                      font = list(family = 'Arial', size = 12,
                                  color = 'rgb(248, 248, 255)'),
                      showarrow = FALSE) %>%
      # labeling the first Likert scale (on the top)
      add_annotations(xref = 'x', yref = 'paper',
                      x = c(x1 / 2, x1 + x2 / 2),
                      y = 1.15,
                      text = top_labels,
                      font = list(family = 'Arial', size = 12,
                                  color = 'rgb(67, 67, 67)'),
                      showarrow = FALSE)
    p
  })
  
  ### /SpC vs MSstats Protein list comparison barplot
  ### SAINT + FASTA description table
  # For LFQ
  # Selected Table & Protein based on horizontal barplot filter
  horizontal_barplot_based_values <- reactiveValues(documents = NULL)
  observe({
    if(is.null(event_data("plotly_click")$curveNumber)){
      horizontal_barplot_based_values$table <- protein_saint_table()
    } else if(as.numeric(event_data("plotly_click")$curveNumber) == 0){
      horizontal_barplot_based_values$table <- protein_saint_table()[protein_saint_table()$Protein %in% quant_type_ratio_values$spc_only_protein_list,]
    }else if(event_data("plotly_click")$curveNumber == 1){
      horizontal_barplot_based_values$table <- protein_saint_table()[protein_saint_table()$Protein %in% quant_type_ratio_values$msstats_protein_list,]
    }
    horizontal_barplot_based_values$selected_protein <- as.character(horizontal_barplot_based_values$table[input$protein_saint_table_rows_selected,]$Protein)
  })
  output$protein_saint_table <- renderDataTable({
    horizontal_barplot_based_values$table
  }, filter = 'top', 
  rownames = FALSE, 
  extensions = 'Buttons',
  options = list(colReorder = TRUE,
                 scrollX = TRUE,
                 lengthMenu = list(c(10, 15, 20, 25, -1), c('10', '15', '20', '25', 'All')),
                 dom = 'lBrtip',
                 buttons = c('copy', 'csv', 'excel', 'pdf', 'print')) 
  )
  
  # For Qual
  output$protein_saint_table_qual <- renderDataTable({
    protein_saint_table()
  }, filter = 'top', 
  rownames = FALSE, 
  extensions = 'Buttons',
  options = list(colReorder = TRUE,
                 scrollX = TRUE,
                 lengthMenu = list(c(10, 15, 20, 25, -1), c('10', '15', '20', '25', 'All')),
                 dom = 'lBrtip',
                 buttons = c('copy', 'csv', 'excel', 'pdf', 'print')) 
  )
  
  ### SpC table
  # For LFQ
  output$spc_table <- renderDataTable({
    selected_protein <- horizontal_barplot_based_values$selected_protein
    if (is.null(input$protein_saint_table_rows_selected))
      return(spc_df())
    spc_df()[spc_df()$ProteinID %in% selected_protein, ]
  }, selection = 'single', 
  filter = 'top', 
  rownames = FALSE, 
  extensions = 'Buttons',
  options = list(colReorder = TRUE,
                 scrollX = TRUE,
                 lengthMenu = list(c(10, 15, 20, 25, -1), c('10', '15', '20', '25', 'All')),
                 dom = 'lBrtip',
                 buttons = c('copy', 'csv', 'excel', 'pdf', 'print')) 
  )
  # For Qual
  output$spc_table_qual <- renderDataTable({
    selected_protein <- as.character(protein_saint_table()[input$protein_saint_table_qual_rows_selected,]$Protein)
    if (is.null(input$protein_saint_table_qual_rows_selected))
      return(spc_df())
    spc_df()[spc_df()$ProteinID %in% selected_protein, ]
  }, selection = 'single', 
  filter = 'top', 
  rownames = FALSE, 
  extensions = 'Buttons',
  options = list(colReorder = TRUE,
                 scrollX = TRUE,
                 lengthMenu = list(c(10, 15, 20, 25, -1), c('10', '15', '20', '25', 'All')),
                 dom = 'lBrtip',
                 buttons = c('copy', 'csv', 'excel', 'pdf', 'print')) 
  )
  
  ### NSAF table
  # For LFQ
  output$nsaf_table <- renderDataTable({
    selected_protein <- horizontal_barplot_based_values$selected_protein
    if (is.null(input$protein_saint_table_rows_selected))
      return(nsaf_df())
    nsaf_df()[nsaf_df()$ProteinID %in% selected_protein, ]
  }, selection = 'single', 
  filter = 'top', 
  rownames = FALSE, 
  extensions = 'Buttons',
  options = list(colReorder = TRUE,
                 scrollX = TRUE,
                 lengthMenu = list(c(10, 15, 20, 25, -1), c('10', '15', '20', '25', 'All')),
                 dom = 'lBrtip',
                 buttons = c('copy', 'csv', 'excel', 'pdf', 'print')) 
  )
  # For Qual
  output$nsaf_table_qual <- renderDataTable({
    selected_protein <- as.character(protein_saint_table()[input$protein_saint_table_qual_rows_selected,]$Protein)
    if (is.null(input$protein_saint_table_qual_rows_selected))
      return(nsaf_df())
    nsaf_df()[nsaf_df()$ProteinID %in% selected_protein, ]
  }, selection = 'single', 
  filter = 'top', 
  rownames = FALSE, 
  extensions = 'Buttons',
  options = list(colReorder = TRUE,
                 scrollX = TRUE,
                 lengthMenu = list(c(10, 15, 20, 25, -1), c('10', '15', '20', '25', 'All')),
                 dom = 'lBrtip',
                 buttons = c('copy', 'csv', 'excel', 'pdf', 'print')) 
  )
  
  ### sequence viewer
  # For LFQ
  output$protein_seq <- renderSequenceViewer({
    sequenceViewer(fasta_df()[fasta_df()$id == horizontal_barplot_based_values$selected_protein,]$Sequence, 
                   horizontal_barplot_based_values$selected_protein)
  })
  # For Qual
  output$protein_seq_qual <- renderSequenceViewer({
    sequenceViewer(fasta_df()[fasta_df()$id == as.character(protein_saint_table()[input$protein_saint_table_qual_row_last_clicked,]$Protein),]$Sequence, 
                   as.character(protein_saint_table()[input$protein_saint_table_qual_row_last_clicked,]$Protein))
  })
  ### Feature Viewer
  ## For LFQ
  # feature viewer dataframe parsing
  proteinFeature_df <- reactive(
    feature_parser(protein_name = horizontal_barplot_based_values$selected_protein, 
                   fasta_df = fasta_df(), 
                   psm_table = psm_df()[psm_df()$`percolator q-value`<=input$peptide_qvalue_cutoff, , drop=FALSE], 
                   exp_design_for_pViz = exp_design_for_pViz())
  )
  # feature viewer - sequence coverage
  output$protein_feature <- renderFeatureViewer({
    featureViewer(feature_df = proteinFeature_df(),
                  sequence = as.character(fasta_df()[fasta_df()$id == horizontal_barplot_based_values$selected_protein,]$Sequence))
  })
  ## For Qual
  # feature viewer dataframe parsing
  proteinFeature_df_qual <- reactive(
    feature_parser(protein_name = as.character(protein_saint_table()[input$protein_saint_table_qual_row_last_clicked,]$Protein), 
                   fasta_df = fasta_df(), 
                   psm_table = psm_df()[psm_df()$`percolator q-value`<=input$peptide_qvalue_cutoff_qual, , drop=FALSE], 
                   exp_design_for_pViz = exp_design_for_pViz())
  )
  # feature viewer - sequence coverage
  output$protein_feature_qual <- renderFeatureViewer({
    featureViewer(feature_df = proteinFeature_df_qual(),
                  sequence = as.character(fasta_df()[fasta_df()$id == as.character(protein_saint_table()[input$protein_saint_table_qual_row_last_clicked,]$Protein),]$Sequence))
  })
  
  ### Unique Peptide Table
  ## For LFQ
  # processing selected protein and selected peptide table
  unique_peptide_df <- reactive({
    unique_peptide_table <- feature_parser(protein_name = horizontal_barplot_based_values$selected_protein, 
                                           fasta_df = fasta_df(), 
                                           psm_table = psm_df()[psm_df()$`percolator q-value`<=input$peptide_qvalue_cutoff, , drop=FALSE], 
                                           exp_design_for_pViz = exp_design_for_pViz())
    if (is.null(input$protein_feature_selected))
      return(unique_peptide_table)
    unique_peptide_table$x <- as.numeric(unique_peptide_table$x)
    unique_peptide_table$y <- as.numeric(unique_peptide_table$y)
    setnames(unique_peptide_table,
             old = c("category", "description", "x", "y"),
             new = c("Condition", "Peptide Sequence", "Start Position", "End Position"))
    selected_peptide_n_terminal_position <- unique_peptide_table[unique_peptide_table$`Peptide Sequence`==input$protein_feature_selected$description,]$`Start Position`
    selected_peptide_c_terminal_position <- unique_peptide_table[unique_peptide_table$`Peptide Sequence`==input$protein_feature_selected$description,]$`End Position`
    unique_peptide_table[!(selected_peptide_n_terminal_position > unique_peptide_table$`End Position` | selected_peptide_c_terminal_position < unique_peptide_table$`Start Position`),]
  })
  
  # peptide table column
  output$unique_peptide_columns_ui <- renderUI({
    checkboxGroupInput("unique_peptide_columns",
                       "Select Columns",
                       choices=names(psm_df()),
                       selected=pre_selected_columns, # Using psm common pre_selected_columns
                       inline = T)
  })
  # unique Peptide table
  output$unique_peptide_table <- renderDataTable({
    # Processing unique peptide per run
    selected_peptide_list <- unique(as.vector(unique_peptide_df()$`Peptide Sequence`))
    selected_protein <- horizontal_barplot_based_values$selected_protein
    filtered_psm_table <- psm_df()[psm_df()$`percolator q-value`<=input$peptide_qvalue_cutoff & grepl(escapeRegex(selected_protein), psm_df()$`protein id`), , drop=FALSE]
    unique_psm_table <- filtered_psm_table[order(filtered_psm_table$file, filtered_psm_table$sequence, filtered_psm_table$`percolator q-value`), ] #sort by id and reverse of abs(value)
    unique_psm_table <- unique_psm_table[!duplicated(unique_psm_table[,c("file", "sequence")]),]
    # Column selection
    unique_peptide_visible_columns <- pre_selected_columns
    if (!is.null(input$unique_peptide_columns)){
      unique_peptide_visible_columns <- input$unique_peptide_columns
    }
    if (is.null(selected_peptide_list))
      return(unique_psm_table[,unique_peptide_visible_columns])
    unique_psm_table[unique_psm_table$sequence %in% selected_peptide_list, unique_peptide_visible_columns, drop=FALSE]
  }, selection = 'single', 
  filter = 'top', 
  rownames = FALSE, 
  extensions = 'Buttons',
  options = list(colReorder = TRUE,
                 scrollX = TRUE,
                 lengthMenu = list(c(10, 15, 20, 25, -1), c('10', '15', '20', '25', 'All')),
                 dom = 'lBrtip',
                 buttons = c('copy', 'csv', 'excel', 'pdf', 'print')) 
  )
  
  ## For Qual
  # processing selected protein and selected peptide table
  unique_peptide_df_qual <- reactive({
    unique_peptide_table <- feature_parser(protein_name = as.character(protein_saint_table()[input$protein_saint_table_qual_row_last_clicked,]$Protein), 
                                           fasta_df = fasta_df(), 
                                           psm_table = psm_df()[psm_df()$`percolator q-value`<=input$peptide_qvalue_cutoff_qual, , drop=FALSE], 
                                           exp_design_for_pViz = exp_design_for_pViz())
    if (is.null(input$protein_feature_qual_selected))
      return(unique_peptide_table)
    unique_peptide_table$x <- as.numeric(unique_peptide_table$x)
    unique_peptide_table$y <- as.numeric(unique_peptide_table$y)
    setnames(unique_peptide_table,
             old = c("category", "description", "x", "y"),
             new = c("Condition", "Peptide Sequence", "Start Position", "End Position"))
    selected_peptide_n_terminal_position <- unique_peptide_table[unique_peptide_table$`Peptide Sequence`==input$protein_feature_qual_selected$description,]$`Start Position`
    selected_peptide_c_terminal_position <- unique_peptide_table[unique_peptide_table$`Peptide Sequence`==input$protein_feature_qual_selected$description,]$`End Position`
    unique_peptide_table[!(selected_peptide_n_terminal_position > unique_peptide_table$`End Position` | selected_peptide_c_terminal_position < unique_peptide_table$`Start Position`),]
  })
  
  # peptide table column
  output$unique_peptide_columns_ui_qual <- renderUI({
    checkboxGroupInput("unique_peptide_columns_qual",
                       "Select Columns",
                       choices=names(psm_df()),
                       selected=pre_selected_columns, # Using psm common pre_selected_columns
                       inline = T)
  })
  # unique Peptide table
  output$unique_peptide_table_qual <- renderDataTable({
    # Processing unique peptide per run
    selected_peptide_list <- unique(as.vector(unique_peptide_df_qual()$`Peptide Sequence`))
    selected_protein <- protein_saint_table()[input$protein_saint_table_qual_row_last_clicked,]$Protein
    filtered_psm_table <- psm_df()[psm_df()$`percolator q-value`<=input$peptide_qvalue_cutoff_qual & grepl(escapeRegex(selected_protein), psm_df()$`protein id`), , drop=FALSE]
    unique_psm_table <- filtered_psm_table[order(filtered_psm_table$file, filtered_psm_table$sequence, filtered_psm_table$`percolator q-value`), ] #sort by id and reverse of abs(value)
    unique_psm_table <- unique_psm_table[!duplicated(unique_psm_table[,c("file", "sequence")]),]
    # Column selection
    unique_peptide_visible_columns <- pre_selected_columns
    if (!is.null(input$unique_peptide_columns_qual)){
      unique_peptide_visible_columns <- input$unique_peptide_columns_qual
    }
    if (is.null(selected_peptide_list))
      return(unique_psm_table[,unique_peptide_visible_columns])
    unique_psm_table[unique_psm_table$sequence %in% selected_peptide_list, unique_peptide_visible_columns, drop=FALSE]
  }, selection = 'single', 
  filter = 'top', 
  rownames = FALSE, 
  extensions = 'Buttons',
  options = list(colReorder = TRUE,
                 scrollX = TRUE,
                 lengthMenu = list(c(10, 15, 20, 25, -1), c('10', '15', '20', '25', 'All')),
                 dom = 'lBrtip',
                 buttons = c('copy', 'csv', 'excel', 'pdf', 'print')) 
  )
  
  ### Peptide SpC
  ## For LFQ
  # Peptide SpC table
  output$pep_spc_table <- renderDataTable({
    # only works for single peptide selection
    selected_peptide_list <- input$protein_feature_selected$description
    psm_selected_peptide <- psm_df()[psm_df()$`percolator q-value`<=input$peptide_qvalue_cutoff & psm_df()$sequence %in% selected_peptide_list, c("file", "sequence"), drop=FALSE]
    filename_list <- as.data.frame(unique(psm_df()$file))
    colnames(filename_list) <- c("file")
    peptide_spc <- merge(x = ddply(psm_selected_peptide, .(file,sequence), nrow), y = filename_list, by="file", all = TRUE)
    dcast(peptide_spc, sequence ~ file)
  }, selection = 'single', 
  filter = 'top', 
  rownames = FALSE, 
  extensions = 'Buttons',
  options = list(colReorder = TRUE,
                 scrollX = TRUE,
                 lengthMenu = list(c(10, 15, 20, 25, -1), c('10', '15', '20', '25', 'All')),
                 dom = 'lBrtip',
                 buttons = c('copy', 'csv', 'excel', 'pdf', 'print')) 
  )
  ## For Qual
  output$pep_spc_table_qual <- renderDataTable({
    # only works for single peptide selection
    selected_peptide_list <- input$protein_feature_qual_selected$description
    psm_selected_peptide <- psm_df()[psm_df()$`percolator q-value`<=input$peptide_qvalue_cutoff_qual & psm_df()$sequence %in% selected_peptide_list, c("file", "sequence"), drop=FALSE]
    filename_list <- as.data.frame(unique(psm_df()$file))
    colnames(filename_list) <- c("file")
    peptide_spc <- merge(x = ddply(psm_selected_peptide, .(file,sequence), nrow), y = filename_list, by="file", all = TRUE)
    dcast(peptide_spc, sequence ~ file)
  }, selection = 'single', 
  filter = 'top', 
  rownames = FALSE, 
  extensions = 'Buttons',
  options = list(colReorder = TRUE,
                 scrollX = TRUE,
                 lengthMenu = list(c(10, 15, 20, 25, -1), c('10', '15', '20', '25', 'All')),
                 dom = 'lBrtip',
                 buttons = c('copy', 'csv', 'excel', 'pdf', 'print')) 
  )
  ##### /Protein
  
  ##### Volcano Plot
  # Select Experiment UI
  output$exp_label <- renderUI({
    selectInput("exp_label_selected",
                label = "Condition:",
                choices = as.character(unique(MSstats_comparison_df()$Label)),
                width = '800px')
  })
  
  # Main Volcano Plot with Plotly interactivity
  output$volcano_plot <- renderPlotly({
    volcanoData <- reactive(
      volcano_preprocess(MSstats_comparison_df(),
                         input$exp_label_selected,
                         input$volcano_plot_pvalue_cutoff,
                         input$log2FC_cutoff)
    )
    # manual coloring group
    volcano_group.colors <- c(NotSignificant = "#999999", 
                              Significant = "#d6d6d6", 
                              FoldChange = "#7f7f7f", 
                              SignificantANDNegativeFoldChange = "#0072B2", 
                              SignificantANDPositiveFoldChange = "#D55E00")
    p <- ggplot(volcanoData(),
                aes_string(x=colnames(volcanoData())[grepl("FC", colnames(volcanoData()))], 
                           y="-log2(adj.pvalue)",
                           key="Protein",
                           color="group",
                           label="adj.pvalue",
                           label2="logFC_str")) + 
      geom_point(size = 0.5, position = position_jitter(w = 0.6,h = 0.6)) +
      scale_color_manual(values=volcano_group.colors) +
      ggtitle(input$exp_label_selected) +
      theme(legend.text=element_text(size=8), legend.title=element_text(size=10)) +
      theme(plot.title = element_text(hjust=0.5, lineheight=.8, face="bold", size=12))
    ggplotly(p, tooltip = c("key", "adj.pvalue", "logFC_str")) %>% layout(dragmode = "select")
  })
  
  # Condition Bar Plot from Volcano Plot's dot selection
  output$condition_plot_from_volcano <- renderPlot({
    if (input$error_bar_value == 'CI'){
      error_value_colname <- colnames(MSstats_condition_plot_df())[grepl("CI", colnames(MSstats_condition_plot_df()))]
    }
    else if (input$error_bar_value == 'SD'){
      error_value_colname <- colnames(MSstats_condition_plot_df())[grepl("SD", colnames(MSstats_condition_plot_df()))]
    }
    
    ggplot(data=MSstats_condition_plot_df()[MSstats_condition_plot_df()$Protein %in% event_data("plotly_click")$key,],
           aes(x=Condition, y=Mean)) +
      geom_errorbar(width=.1,
                    aes_string(ymin=paste("Mean",error_value_colname, sep = '-'),
                               ymax=paste("Mean",error_value_colname, sep = '+')),
                    colour="red") +
      geom_point(size=4,
                 shape=21,
                 fill="white") +
      ggtitle(event_data("plotly_click")$key) +
      theme(plot.title = element_text(hjust=0.5, lineheight=.8, face="bold", size=12))
  })
  
  # Vocano Plot's Dataframe
  volcano_comparison_df <- reactive({
    logFC_colname <- colnames(MSstats_comparison_df())[grepl("FC", colnames(MSstats_comparison_df()))]
    MSstats_comparison_df()[MSstats_comparison_df()$adj.pvalue <= input$volcano_plot_pvalue_cutoff & abs(MSstats_comparison_df()[logFC_colname]) >= input$log2FC_cutoff & !is.infinite(MSstats_comparison_df()[,logFC_colname]) & !is.na(MSstats_comparison_df()[,"adj.pvalue"]), 
                            c("Protein", "Label", logFC_colname, "pvalue", "adj.pvalue", "SE", "Tvalue", "DF", "description", "Gene ID")]
  })
  
  # Volcano Plot's DataTable
  output$comparison_table <- renderDataTable({
    if (is.null(event_data("plotly_click")$key))
      return(volcano_comparison_df())
    volcano_comparison_df()[volcano_comparison_df()$Protein %in% event_data("plotly_click")$key,] ###############################
  }, selection = 'single', 
  filter = 'top', 
  rownames = FALSE,
  extensions = 'Buttons',
  options = list(lengthMenu = list(c(10, 15, 20, 25, -1), c('10', '15', '20', '25', 'All')),
                 dom = 'lBrtip',
                 buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                 autoWidth = TRUE,
                 columnDefs = list(list(width = '200px', targets = 1)))
  )
  
  # Selected Protein's Volcano Plot - Static
  output$selected_protein_volcano_plot <- renderPlot({
    logFC_colname <- colnames(MSstats_comparison_df())[grepl("FC", colnames(MSstats_comparison_df()))]
    filtered_MSstats_comparison_df <- volcano_comparison_df()
    if (is.null(event_data("plotly_click")$key)){
      volcano_plot_from_table_title <- as.character(filtered_MSstats_comparison_df[input$comparison_table_row_last_clicked,]$Protein)
      selectedVolcanoData <- reactive(
        selected_volcano_preprocess(MSstats_comparison_df(),
                                    MSstats_comparison_df()[MSstats_comparison_df()$Protein == volcano_plot_from_table_title,],
                                    input$volcano_plot_pvalue_cutoff,
                                    input$log2FC_cutoff)
      )
    }else{
      selectedVolcanoData <- reactive(
        selected_volcano_preprocess(MSstats_comparison_df(),
                                    MSstats_comparison_df()[MSstats_comparison_df()$Protein == event_data("plotly_click")$key,],
                                    input$volcano_plot_pvalue_cutoff,
                                    input$log2FC_cutoff)
      )
      volcano_plot_from_table_title <- event_data("plotly_click")$key
    }
    
    # manual coloring group
    volcano_group.colors <- c(NotSignificant = "#999999", 
                              Significant = "#d6d6d6", 
                              FoldChange = "#7f7f7f", 
                              SignificantANDNegativeFoldChange = "#0072B2", 
                              SignificantANDPositiveFoldChange = "#D55E00")
    # Inf value handling: set biggest finite value + 1
    logFCvector <- MSstats_comparison_df()[,logFC_colname]
    logFCvector <- logFCvector[!is.na(logFCvector)]
    finiteFCMax <- max(logFCvector[logFCvector!=max(logFCvector)]) + 1
    finiteFCMin <- min(logFCvector[logFCvector!=min(logFCvector)]) - 1
    # adj.pvalue = 0 to the lowest adj.pvalue 
    adjPvalueVector <- MSstats_comparison_df()[, "adj.pvalue"]
    adjPvalueVector <- adjPvalueVector[!is.na(adjPvalueVector)]
    finiteMinAdjPvalue <- min(adjPvalueVector[adjPvalueVector!=min(adjPvalueVector)]) / 2
    
    ggplot(selectedVolcanoData(),
           aes_string(x=colnames(selectedVolcanoData())[grepl("FC", colnames(selectedVolcanoData()))], 
                      y="-log2(adj.pvalue)",
                      key="Protein",
                      color="group")) + 
      geom_point() +
      scale_color_manual(values=volcano_group.colors) + 
      coord_cartesian(ylim=c(0, -log2(finiteMinAdjPvalue)), xlim=c(finiteFCMin, finiteFCMax)) +
      geom_text_repel(aes(label=Label)) +
      ggtitle(volcano_plot_from_table_title) +
      theme(plot.title = element_text(hjust=0.5, lineheight=.8, face="bold", size=12))
  })
  
  # Condition Bar Plot from table row selection - static
  output$condition_plot_from_table <- renderPlot({
    filtered_MSstats_comparison_df <- volcano_comparison_df()
    if (input$error_bar_value == 'CI'){
      error_value_colname <- colnames(MSstats_condition_plot_df())[grepl("CI", colnames(MSstats_condition_plot_df()))]
    }
    else if (input$error_bar_value == 'SD'){
      error_value_colname <- colnames(MSstats_condition_plot_df())[grepl("SD", colnames(MSstats_condition_plot_df()))]
    }
    selected_condition_plot_data <- NULL
    if(is.null(event_data("plotly_click")$key)){
      selected_condition_plot_data <- NULL
      condition_plot_from_table_title <- as.character(filtered_MSstats_comparison_df[input$comparison_table_row_last_clicked,]$Protein)
      selected_condition_plot_data <- MSstats_condition_plot_df()[MSstats_condition_plot_df()$Protein %in% condition_plot_from_table_title,]
      
    }else{
      selected_condition_plot_data <- MSstats_condition_plot_df()[MSstats_condition_plot_df()$Protein %in% event_data("plotly_click")$key,]
      condition_plot_from_table_title <- event_data("plotly_click")$key
    }
    ggplot(data=selected_condition_plot_data,
           aes(x=Condition, y=Mean)) +
      geom_errorbar(width=.1,
                    aes_string(ymin=paste("Mean",error_value_colname, sep = '-'),
                               ymax=paste("Mean",error_value_colname, sep = '+')),
                    colour="red") +
      geom_point(size=4,
                 shape=21,
                 fill="white") +
      ggtitle(condition_plot_from_table_title) +
      theme(plot.title = element_text(hjust=0.5, lineheight=.8, face="bold", size=12))
  })
  
  # Combined Dataframe for outlier Inspector
  spc_saint_msstats <- reactive({
    ##### Comparison of MSStats comparison_csv and SpC and SAINT
    # combine SAINT and SpC 
    saint_df <- protein_saint_table()
    saint_df$description <- NULL
    names(saint_df) <- gsub('SAINT.', '', names(saint_df))
    melt_saint <- melt(saint_df, id="Protein")
    colnames(melt_saint)[colnames(melt_saint)=="variable"] <- "Condition"
    colnames(melt_saint)[colnames(melt_saint)=="value"] <- "SAINT.AvgP"
    colnames(melt_saint)[colnames(melt_saint)=="Prey"] <- "ProteinID"
    melt_spc <- melt(spc_df(), id="ProteinID")
    colnames(melt_spc)[colnames(melt_spc)=="variable"] <- "Original File Name"
    colnames(melt_spc)[colnames(melt_spc)=="value"] <- "SpC"
    melt_spc <- merge(x=melt_spc, y=qual_experiment_design_df()[, c("Original File Name", "Biological Condition")], by="Original File Name", all=TRUE)
    melt_spc$`Original File Name` <- NULL
    colnames(melt_spc)[colnames(melt_spc)=="Biological Condition"] <- "Condition"
    melt_spc <- aggregate(melt_spc[,"SpC"], list(melt_spc$ProteinID, melt_spc$Condition), mean)
    colnames(melt_spc) <- c("Protein", "Condition", "AvgSpC")
    melt_saint_spc <- merge(x=melt_spc, y=melt_saint, by=c("Protein", "Condition"), all=TRUE)
    # combine SAINT.SpC and MSStats comparison_csv
    colnames(melt_saint_spc) <- c("Protein", "Label", "AvgSpC", "SAINT.AvgP")
    spc_saint_msstats <- merge(x=melt_saint_spc, y=MSstats_comparison_df(), by=c("Protein", "Label"), all.y=TRUE)
    spc_saint_msstats$Var.5 <- NULL
    spc_saint_msstats$description <- NULL
    # DSR SpC calculation
    test_control_map <- unique(qual_experiment_design_df()[,c("Biological Condition", "Test or Control")])
    colnames(test_control_map) <- c("Condition", "TestControl")
    TestCondition <- as.character(test_control_map[test_control_map$TestControl=="T", "Condition"])
    ControlCondition <- as.character(test_control_map[test_control_map$TestControl=="C", "Condition"])
    TestDataset <- melt_spc[melt_spc$Condition %in% TestCondition, ]
    ControlDataset <- melt_spc[melt_spc$Condition %in% ControlCondition, ]
    calc_dataset <- merge(x=TestDataset, y=ControlDataset, by="Protein", all=TRUE)
    calc_dataset$DSR_SpC <- (calc_dataset$AvgSpC.x - calc_dataset$AvgSpC.y - 0.5)/sqrt(calc_dataset$AvgSpC.x^2 + calc_dataset$AvgSpC.y^2)
    # combine SAINT.SpC.MSStats and DSR SpC
    spc_saint_msstats_dsrspc <- merge(x=spc_saint_msstats, y=calc_dataset[,c("Protein","Condition.x","DSR_SpC")], by.x=c("Protein", "Label"), by.y=c("Protein", "Condition.x"), all.x=TRUE)
    spc_saint_msstats_dsrspc$Protein <- as.character(spc_saint_msstats_dsrspc$Protein)
    # LogFC and adj.pvalue string columns
    logFC_colname <- colnames(spc_saint_msstats_dsrspc)[grepl("FC", colnames(spc_saint_msstats_dsrspc))]
    logFC_str <- "FoldChange_str"
    spc_saint_msstats_dsrspc[,logFC_str] <- as.character(spc_saint_msstats_dsrspc[,logFC_colname])
    spc_saint_msstats_dsrspc[is.na(spc_saint_msstats_dsrspc[,logFC_str]),logFC_str] <- "NA"
    spc_saint_msstats_dsrspc[,"adj.pvalue_str"] <- as.character(spc_saint_msstats_dsrspc[,"adj.pvalue"])
    spc_saint_msstats_dsrspc[is.na(spc_saint_msstats_dsrspc[,"adj.pvalue_str"]),"adj.pvalue_str"] <- "NA"
    spc_saint_msstats_dsrspc <- merge(x=spc_saint_msstats_dsrspc,y=fasta_df(),by.x="Protein", by.y="id", all.x=TRUE)
    spc_saint_msstats_dsrspc$description <- as.character(spc_saint_msstats_dsrspc$description)
    # Filter to include only either inf or NA
    logFC_colname <- colnames(MSstats_comparison_df())[grepl("FC", colnames(MSstats_comparison_df()))]
    spc_saint_msstats_dsrspc <- spc_saint_msstats_dsrspc[is.infinite(spc_saint_msstats_dsrspc[,logFC_colname]) | is.na(spc_saint_msstats_dsrspc[,"adj.pvalue"]),]
    spc_saint_msstats_dsrspc
  })
  
  # Outlier Inspector
  output$comparison_outlier_table <- renderDataTable({
    logFC_colname <- colnames(MSstats_comparison_df())[grepl("FC", colnames(MSstats_comparison_df()))]
    df <- spc_saint_msstats()[ ,c("Protein","Label","DSR_SpC","SAINT.AvgP","adj.pvalue_str","FoldChange_str","AvgSpC","issue","description","DF","adj.pvalue","pvalue",logFC_colname,"SE","Tvalue","MissingPercentage","ImputationPercentage")]
    if (is.null(event_data("plotly_click")$key))
      return(df)
    df[df$Protein %in% event_data("plotly_click")$key,] ###############################
    
    
  }, selection = 'single', 
  filter = 'top', 
  rownames = FALSE, 
  extensions = 'Buttons',
  options = list(colReorder = TRUE,
                 autoWidth = TRUE,
                 columnDefs = list(list(width = '200px', targets = 1)),
                 scrollX = TRUE,
                 lengthMenu = list(c(10, 15, 20, 25, -1), c('10', '15', '20', '25', 'All')),
                 dom = 'lBrtip',
                 buttons = c('copy', 'csv', 'excel', 'pdf', 'print'))
  )
  
  # selected outlier protein's quant barplot - MSstats_quant_df
  output$selected_outlier_protein_quant_barplot <- renderPlot({
    selected_outlier_protein <- spc_saint_msstats()[input$comparison_outlier_table_row_last_clicked,]$Protein
    quantification_csv <- MSstats_quant_df()
    quantification_csv[,1] <- NULL
    quantification_csv[is.na(quantification_csv)] <- 0
    melt_df <- melt(quantification_csv[quantification_csv$Protein %in% selected_outlier_protein,], id.vars = c("Protein"))
    ggplot(melt_df, aes(variable, value)) +
      geom_bar(aes(fill = variable), position = "dodge", stat="identity") +
      ggtitle(selected_outlier_protein) +
      labs(x = "Run", y = "Quant") +
      theme(plot.title = element_text(hjust=0.5, lineheight=.8, face="bold", size=12), legend.title=element_blank())
    
  })
  
  # Enrichr
  observeEvent(input$enrichr_action, {
    logFC_colname <- colnames(volcano_comparison_df())[grepl("FC", colnames(volcano_comparison_df()))]
    if (input$enrichr_grouping=="Positive"){
      symbol_IDs <- volcano_comparison_df()[volcano_comparison_df()[,"Label"]==input$exp_label_selected & volcano_comparison_df()[,logFC_colname]>0, "Gene ID"]
    }else if (input$enrichr_grouping=="Negative"){
      symbol_IDs <- volcano_comparison_df()[volcano_comparison_df()[,"Label"]==input$exp_label_selected & volcano_comparison_df()[,logFC_colname]<0, "Gene ID"]
    }else if (input$enrichr_grouping=="Both"){
      symbol_IDs <- volcano_comparison_df()[volcano_comparison_df()[,"Label"]==input$exp_label_selected, "Gene ID"]
    }
    enrichr_list <- paste(symbol_IDs, collapse="\n")
    session$sendCustomMessage(type='enrich', message = list(list=enrichr_list, popup=TRUE))
  })
  ##### /Volcano Plot
  
  ##### Heatmap
  # ComplexHeatmap
  output$complexHeatmapPlot <- renderPlot({
    withProgress(message = "Generating the Heatmap", {
      complex_heatmap(input$unidirectional,
                      input$quantCentering,
                      input$heatmap_pvalue_cutoff,
                      quant_experiment_design_df(),
                      MSstats_comparison_df(),
                      MSstats_quant_df(),
                      input$manual_k)
    })
  })
  ##### /Heatmap
  
  ##### PSM-ID
  ## For Common
  # For Pre-selection of columns on PSM table 
  pre_selected_columns = c("file", 
                           "scan", 
                           "charge", 
                           "spectrum precursor m/z", 
                           "percolator q-value", 
                           "protein id", 
                           "sequence", 
                           "ppm")
  

  ## For LFQ
  # PSM Table selection UI: http://stackoverflow.com/questions/36784906/shiny-allowling-users-to-choose-which-columns-to-display // not selected answer is the better solution - simple and effective
  # Actual column selection UI handling for PSM-ID
  output$psm_id_columns_ui <- renderUI({
    checkboxGroupInput("psm_id_columns",
                       "Select Columns",
                       choices=names(psm_df()),
                       selected=pre_selected_columns,
                       inline = T)
  })
  ## For Qual
  output$psm_id_columns_ui_qual <- renderUI({
    checkboxGroupInput("psm_id_columns",
                       "Select Columns",
                       choices=names(psm_df()),
                       selected=pre_selected_columns,
                       inline = T)
  })
  
  # Download PSM table as a CSV file
  ## For LFQ
  output$download_psm_table <- downloadHandler(
    filename = function() {
      "psm-table.csv"
    },
    content = function(file) {
      psm_visible_columns <- pre_selected_columns
      if (!is.null(input$psm_id_columns)){
        psm_visible_columns <- input$psm_id_columns
      }
      write.csv(psm_df()[psm_df()$`percolator q-value`<=input$psm_id_qvalue_cutoff, psm_visible_columns, drop=FALSE], file)
    }
  )
  ## For Qual
  output$download_psm_table_qual <- downloadHandler(
    filename = function() {
      "psm-table.csv"
    },
    content = function(file) {
      psm_visible_columns <- pre_selected_columns
      if (!is.null(input$psm_id_columns)){
        psm_visible_columns <- input$psm_id_columns
      }
      write.csv(psm_df()[psm_df()$`percolator q-value`<=input$psm_id_qvalue_cutoff_qual, psm_visible_columns, drop=FALSE], file)
    }
  )
  
  # PSM Table for PSM-ID - column selection dependent
  ## For LFQ
  output$psm_id_table <- renderDataTable({
    psm_visible_columns <- pre_selected_columns
    if (!is.null(input$psm_id_columns)){
      psm_visible_columns <- input$psm_id_columns
    }
    psm_df()[psm_df()$`percolator q-value`<=input$psm_id_qvalue_cutoff, psm_visible_columns, drop=FALSE]
  }, selection = 'single', 
  filter = 'top', 
  rownames = FALSE, 
  extensions = 'ColReorder',
  options = list(colReorder = TRUE,
                 scrollX = TRUE)
  )
  ## For Qual
  output$psm_id_table_qual <- renderDataTable({
    psm_visible_columns <- pre_selected_columns
    if (!is.null(input$psm_id_columns)){
      psm_visible_columns <- input$psm_id_columns
    }
    psm_df()[psm_df()$`percolator q-value`<=input$psm_id_qvalue_cutoff_qual, psm_visible_columns, drop=FALSE]
  }, selection = 'single', 
  filter = 'top', 
  rownames = FALSE, 
  extensions = 'ColReorder',
  options = list(colReorder = TRUE,
                 scrollX = TRUE)
  )
  
  # Selected Peptide's all PSM across runs
  ## For LFQ
  output$selected_psm_id_table <- renderDataTable({
    filtered_psm_df <- psm_df()[psm_df()$`percolator q-value`<=input$psm_id_qvalue_cutoff, ]
    psm_df()[which(psm_df()$sequence == as.character(filtered_psm_df[row.names(filtered_psm_df[input$psm_id_table_row_last_clicked, ]), ]$sequence)), c("file", "scan")]
  }, selection = 'single', 
  options = list(pageLength=25)
  )
  ## For Qual
  output$selected_psm_id_table_qual <- renderDataTable({
    filtered_psm_df <- psm_df()[psm_df()$`percolator q-value`<=input$psm_id_qvalue_cutoff_qual, ]
    psm_df()[which(psm_df()$sequence == as.character(filtered_psm_df[row.names(filtered_psm_df[input$psm_id_table_qual_row_last_clicked, ]), ]$sequence)), c("file", "scan")]
  }, selection = 'single', 
  options = list(pageLength=25)
  )
  
  # Lorikeet plot - Retrieve data from HDF5
  ## For LFQ
  output$ms2lorikeet <- renderShinylorikeet({
    # initial variable assignments
    filtered_psm_df <- psm_df()[psm_df()$`percolator q-value`<=input$psm_id_qvalue_cutoff, ]
    selected_peptide_index <- row.names(psm_df()[which(psm_df()$sequence == as.character(filtered_psm_df[row.names(filtered_psm_df[input$psm_id_table_row_last_clicked, ]), ]$sequence)), ][input$selected_psm_id_table_row_last_clicked,])
    modSeq <- as.character(psm_df()[selected_peptide_index,]$sequence)
    sequence <- as.character(psm_df()[selected_peptide_index,]$`unmodified sequence`)
    scanNum <- psm_df()[selected_peptide_index,]$scan
    fileName <- as.character(psm_df()[selected_peptide_index,]$file)
    charge <- psm_df()[selected_peptide_index,]$charge
    precursorMz <- psm_df()[selected_peptide_index,]$`spectrum precursor m/z`
    
    # peak processing for lorikeet
    peaks <- h5read(ms2_hdf(), paste(psm_df()[selected_peptide_index,]$file, psm_df()[selected_peptide_index,]$scan-1, sep="/"))
    
    # Obtain modification notation indexes
    left_bracket <- unlist(gregexpr(pattern ='\\[',modSeq))
    right_bracket <- unlist(gregexpr(pattern ='\\]',modSeq))
    
    if(!grepl("\\[", modSeq)){
      variableMods <- ""
    }else{
      # variableMods dataframe setup
      variableMods <- cbind.data.frame(left_bracket, right_bracket)
      variableMods["index"] <- NA
      variableMods["modMass"] <- NA
      variableMods["aminoAcid"] <- NA
      
      # Processing index, aminoacid letter, and modMass
      previous_right_bracket = 0
      previous_index = 0
      for (i in 1:nrow(variableMods)){
        if((variableMods[row.names(variableMods)==i,]$left_bracket-previous_right_bracket)==1){
          #print("Nterminal Mod!")
          variableMods[row.names(variableMods)==i,]$modMass <- as.numeric(substr(modSeq, variableMods[row.names(variableMods)==i,]$left_bracket + 2, variableMods[row.names(variableMods)==i,]$right_bracket - 1))
          variableMods[row.names(variableMods)==i,]$index <- 1
          variableMods[row.names(variableMods)==i,]$aminoAcid <- substr(sequence, variableMods[row.names(variableMods)==i,]$index, variableMods[row.names(variableMods)==i,]$index)
          previous_right_bracket <- variableMods[row.names(variableMods)==i,]$right_bracket
          
        }else if((variableMods[row.names(variableMods)==i,]$left_bracket-previous_right_bracket)==2){
          #print("first aa")
          variableMods[row.names(variableMods)==i,]$modMass <- as.numeric(substr(modSeq, variableMods[row.names(variableMods)==i,]$left_bracket + 2, variableMods[row.names(variableMods)==i,]$right_bracket - 1))
          variableMods[row.names(variableMods)==i,]$index <- 1
          variableMods[row.names(variableMods)==i,]$aminoAcid <- substr(sequence, variableMods[row.names(variableMods)==i,]$index, variableMods[row.names(variableMods)==i,]$index)
          previous_right_bracket <- variableMods[row.names(variableMods)==i,]$right_bracket
          
        }else{
          #print("something in the middle")
          variableMods[row.names(variableMods)==i,]$modMass <- as.numeric(substr(modSeq, variableMods[row.names(variableMods)==i,]$left_bracket + 2, variableMods[row.names(variableMods)==i,]$right_bracket - 1))
          variableMods[row.names(variableMods)==i,]$index <- variableMods[row.names(variableMods)==i,]$left_bracket - previous_right_bracket + previous_index
          if(previous_right_bracket==0 & previous_index==0){
            variableMods[row.names(variableMods)==i,]$index <- variableMods[row.names(variableMods)==i,]$index - 1
          }
          variableMods[row.names(variableMods)==i,]$aminoAcid <- substr(sequence, variableMods[row.names(variableMods)==i,]$index, variableMods[row.names(variableMods)==i,]$index)
          previous_index <- variableMods[row.names(variableMods)==i,]$index - 1
          previous_right_bracket <- variableMods[row.names(variableMods)==i,]$right_bracket
        }
      }
      variableMods$left_bracket <- NULL
      variableMods$right_bracket <- NULL
    }
    shinylorikeet(sequence, scanNum, charge, precursorMz, fileName, variableMods, peaks)
  })
  ## For Qual
  output$ms2lorikeet_qual <- renderShinylorikeet({
    # initial variable assignments
    filtered_psm_df <- psm_df()[psm_df()$`percolator q-value`<=input$psm_id_qvalue_cutoff_qual, ]
    selected_peptide_index <- row.names(psm_df()[which(psm_df()$sequence == as.character(filtered_psm_df[row.names(filtered_psm_df[input$psm_id_table_qual_row_last_clicked, ]), ]$sequence)), ][input$selected_psm_id_table_qual_row_last_clicked,])
    modSeq <- as.character(psm_df()[selected_peptide_index,]$sequence)
    sequence <- as.character(psm_df()[selected_peptide_index,]$`unmodified sequence`)
    scanNum <- psm_df()[selected_peptide_index,]$scan
    fileName <- as.character(psm_df()[selected_peptide_index,]$file)
    charge <- psm_df()[selected_peptide_index,]$charge
    precursorMz <- psm_df()[selected_peptide_index,]$`spectrum precursor m/z`
    
    # peak processing for lorikeet
    peaks <- h5read(ms2_hdf_qual(), paste(psm_df()[selected_peptide_index,]$file, psm_df()[selected_peptide_index,]$scan-1, sep="/"))
    
    # Obtain modification notation indexes
    left_bracket <- unlist(gregexpr(pattern ='\\[',modSeq))
    right_bracket <- unlist(gregexpr(pattern ='\\]',modSeq))
    
    if(!grepl("\\[", modSeq)){
      variableMods <- ""
    }else{
      # variableMods dataframe setup
      variableMods <- cbind.data.frame(left_bracket, right_bracket)
      variableMods["index"] <- NA
      variableMods["modMass"] <- NA
      variableMods["aminoAcid"] <- NA
      
      # Processing index, aminoacid letter, and modMass
      previous_right_bracket = 0
      previous_index = 0
      for (i in 1:nrow(variableMods)){
        if((variableMods[row.names(variableMods)==i,]$left_bracket-previous_right_bracket)==1){
          #print("Nterminal Mod!")
          variableMods[row.names(variableMods)==i,]$modMass <- as.numeric(substr(modSeq, variableMods[row.names(variableMods)==i,]$left_bracket + 2, variableMods[row.names(variableMods)==i,]$right_bracket - 1))
          variableMods[row.names(variableMods)==i,]$index <- 1
          variableMods[row.names(variableMods)==i,]$aminoAcid <- substr(sequence, variableMods[row.names(variableMods)==i,]$index, variableMods[row.names(variableMods)==i,]$index)
          previous_right_bracket <- variableMods[row.names(variableMods)==i,]$right_bracket
          
        }else if((variableMods[row.names(variableMods)==i,]$left_bracket-previous_right_bracket)==2){
          #print("first aa")
          variableMods[row.names(variableMods)==i,]$modMass <- as.numeric(substr(modSeq, variableMods[row.names(variableMods)==i,]$left_bracket + 2, variableMods[row.names(variableMods)==i,]$right_bracket - 1))
          variableMods[row.names(variableMods)==i,]$index <- 1
          variableMods[row.names(variableMods)==i,]$aminoAcid <- substr(sequence, variableMods[row.names(variableMods)==i,]$index, variableMods[row.names(variableMods)==i,]$index)
          previous_right_bracket <- variableMods[row.names(variableMods)==i,]$right_bracket
          
        }else{
          #print("something in the middle")
          variableMods[row.names(variableMods)==i,]$modMass <- as.numeric(substr(modSeq, variableMods[row.names(variableMods)==i,]$left_bracket + 2, variableMods[row.names(variableMods)==i,]$right_bracket - 1))
          variableMods[row.names(variableMods)==i,]$index <- variableMods[row.names(variableMods)==i,]$left_bracket - previous_right_bracket + previous_index
          if(previous_right_bracket==0 & previous_index==0){
            variableMods[row.names(variableMods)==i,]$index <- variableMods[row.names(variableMods)==i,]$index - 1
          }
          variableMods[row.names(variableMods)==i,]$aminoAcid <- substr(sequence, variableMods[row.names(variableMods)==i,]$index, variableMods[row.names(variableMods)==i,]$index)
          previous_index <- variableMods[row.names(variableMods)==i,]$index - 1
          previous_right_bracket <- variableMods[row.names(variableMods)==i,]$right_bracket
        }
      }
      variableMods$left_bracket <- NULL
      variableMods$right_bracket <- NULL
    }
    shinylorikeet(sequence, scanNum, charge, precursorMz, fileName, variableMods, peaks)
  })
  
  ##### PSM-Quant
  # For Pre-selection of columns on PSM table 
  pre_selected_columns_for_quant = c("charge", 
                                     "spectrum precursor m/z", 
                                     "percolator q-value", 
                                     "protein id", 
                                     "sequence", 
                                     "ppm")
  # Actual column selection UI handling for PSM-Quant
  output$psm_quant_columns_ui <- renderUI({
    checkboxGroupInput("psm_quant_columns",
                       "Select Columns",
                       choices=names(psm_df()),
                       selected=pre_selected_columns_for_quant,
                       inline = T)
  })
  
  # Extract psm rows that can be mapped to chromatogram dataframe
  chrom_psm_df <- reactive({
    mapped_psm_df <- psm_df()[psm_df()$`percolator q-value`<=input$psm_quant_qvalue_cutoff & psm_df()$`Skyline Modified Sequence` %in% peptides_for_chromatograms(), , drop=FALSE]
    unique_peptide_df <- mapped_psm_df[!duplicated(mapped_psm_df[,c("charge","sequence")]),]
  })
  # PSM Table for PSM-Quant - column selection dependent AND filtered those have chromatogram!!!!!!!!!!!!! <- need to implement by sequence match!!!!!!!!!!!!
  output$psm_quant_table <- renderDataTable({
    psm_visible_columns <- pre_selected_columns_for_quant
    if (!is.null(input$psm_quant_columns)){
      psm_visible_columns <- input$psm_quant_columns
    }
    chrom_psm_df()[, psm_visible_columns, drop=FALSE]
  }, selection = 'single', 
  filter = 'top', 
  rownames = FALSE, 
  extensions = 'ColReorder',
  options = list(colReorder = TRUE,
                 scrollX = TRUE)
  )
  
  # selected_peptide_quant_barplot
  output$selected_peptide_quant_barplot <- renderPlot({
    selected_peptide_quant_df <- merged_msstats_skyline_chromatograms_df()[merged_msstats_skyline_chromatograms_df()$PeptideModifiedSequence == as.character(chrom_psm_df()[input$psm_quant_table_row_last_clicked,]$`Skyline Modified Sequence`) & merged_msstats_skyline_chromatograms_df()$annotation_QValue<=input$mprophet_qvalue_cutoff,]
    ggplot(selected_peptide_quant_df, aes(Condition, Area)) +   
      geom_bar(aes(fill = FileName), position = "dodge", stat="identity") +
      ggtitle(as.character(chrom_psm_df()[input$psm_quant_table_row_last_clicked,]$`Skyline Modified Sequence`)) +
      theme(plot.title = element_text(hjust=0.5, lineheight=.8, face="bold", size=12))
  })
  
  # Selected Pepetide associated chromaatogram Plot
  output$selected_chromatogram_plot <- renderPlotly({
    # Chromatogram ggplotly
    psm_selected_sequence <- as.character(chrom_psm_df()[input$psm_quant_table_row_last_clicked,]$`Skyline Modified Sequence`)
    psm_selected_charge <- as.character(chrom_psm_df()[input$psm_quant_table_row_last_clicked,]$`charge`)
    filtered_chrom <- merged_msstats_skyline_chromatograms_df()[merged_msstats_skyline_chromatograms_df()$PeptideModifiedSequence==psm_selected_sequence & merged_msstats_skyline_chromatograms_df()$PrecursorCharge==psm_selected_charge,]
    filtered_boundaries <- boundaries_df()[boundaries_df()$PeptideModifiedSequence==psm_selected_sequence & boundaries_df()$PrecursorCharge==psm_selected_charge,]
    # Remove duplicate rows
    setkey(filtered_chrom, NULL)
    setkey(filtered_boundaries, NULL)
    filtered_chrom <- unique(filtered_chrom)
    filtered_boundaries <- unique(filtered_boundaries)
    # Filter mProphet Q-value
    filtered_boundaries <- merge(x=filtered_boundaries, y=msstats_skyline_input_df()[msstats_skyline_input_df()$PeptideModifiedSequence==psm_selected_sequence,c("PeptideModifiedSequence", "FileName", "annotation_QValue")], by=c("PeptideModifiedSequence", "FileName"), all.x=TRUE)
    filtered_boundaries <- filtered_boundaries[filtered_boundaries$annotation_QValue<=input$mprophet_qvalue_cutoff,]
    # Split Time and intensities into rows
    splited_each_chrom <- cSplit(filtered_chrom, c("Times","Intensities"), ",", "long")
    splited_each_chrom$PrecursorMZ_Color <- as.factor(splited_each_chrom$ProductMz)
    # For log2 intensity plot - 1/2
    if (input$chrom_intensity_transform == "Log2"){
      splited_each_chrom[splited_each_chrom == 0] <- 1
    }
    ## Exp mapping with filename
    df <- msstats_skyline_input_df()[which(msstats_skyline_input_df()$annotation_QValue <= input$unnorm_peptide_dist_violin_plot_mprophet_qvalue),c("FileName", "Area", "annotation_QValue")]
    exp_df <- quant_experiment_design_df()
    # Add filename with extension
    exp_df$FilenameWithExt <- paste(exp_df$`Original File Name`, 
                                    tools::file_ext(levels(splited_each_chrom$FileName)[1]), 
                                    sep=".")
    # Add biological condtion with Biorep # & Crux #
    exp_df$cond_rep <- paste(exp_df$`Biological Condition`,exp_df$BioReplicate, exp_df$`Crux File Integer`, sep="_")
    # facet_labels vector
    facet_labels <- exp_df[order(exp_df$`Biological Condition`, exp_df$FilenameWithExt), "cond_rep"]
    names(facet_labels) <- exp_df[order(exp_df$`Biological Condition`, exp_df$FilenameWithExt), "FilenameWithExt"]
    # Give Factor Order by Biological condition 
    splited_each_chrom$FileName <- factor(splited_each_chrom$FileName,
                                          levels = exp_df[order(exp_df$`Biological Condition`, exp_df$FilenameWithExt), "FilenameWithExt"],
                                          ordered = TRUE)
    
    p <- ggplot(splited_each_chrom, aes(x=Times,
                                        y=Intensities,
                                        group=FileName,
                                        color=PrecursorMZ_Color,
                                        label1=MassErrorPPM,
                                        label2=annotation_QValue)) +
      geom_line() +
      facet_wrap(~FileName, labeller=labeller(FileName = facet_labels), nrow=as.numeric(input$grid_row_num)) +
      geom_vline(aes(xintercept=MinStartTime), filtered_boundaries, linetype=4) +
      geom_vline(aes(xintercept=MaxEndTime), filtered_boundaries, linetype=4) +
      ggtitle(paste("Sequence: ", psm_selected_sequence, ", ", "Charge: ", psm_selected_charge))
    # For log2 intensity plot - 2/2
    if (input$chrom_intensity_transform == "Log2"){
      p <- p + scale_y_continuous(trans=log2_trans())
    }
    ggplotly(p, tooltip = c("Times", "FileName", "MassErrorPPM", "annotation_QValue"))
  })
  ##### /PSM-Quant
  
  ##### Excel Report
  exce_report_holder <- reactiveValues(documents = NULL)
  ## For LFQ
  observeEvent(input$generating_excel_report, {
    withProgress(message = "Building the Excel Report...", value = 5, {
      logFC_colname <- colnames(MSstats_comparison_df())[grepl("FC", colnames(MSstats_comparison_df()))]
      outlier_table <- spc_saint_msstats()[ ,c("Protein","Label","DSR_SpC","SAINT.AvgP","adj.pvalue_str","FoldChange_str","AvgSpC","issue","DF","adj.pvalue","pvalue",logFC_colname,"SE","Tvalue","MissingPercentage","ImputationPercentage")]
      exce_report_holder$filepath <- excel_report_generator_v2(filename = "temp/report.xlsx",
                                                               # Volcano Plot, Conditional Format for Quant Table
                                                               pvalue_cutoff = input$excelPvalue,
                                                               log2FC_cutoff = input$excelLogFC,
                                                               # UpSet Plot
                                                               peptide_q_value = input$excelPeptideQvalue,
                                                               protein_q_value = input$excelProteinQvalue,
                                                               # For heatmap
                                                               unidirectional = input$excelUnidirectional,
                                                               quantCentering = input$excelQuantCentering,
                                                               minPvalueFilter = input$excelPvalueHeatMap,
                                                               cluster_k = input$excelManual_k,
                                                               # data from RData
                                                               experiment_design = qual_experiment_design_df(),
                                                               comparison_csv = MSstats_comparison_df(),
                                                               fasta_df = fasta_df(),
                                                               saint_table = protein_saint_table(),
                                                               spc_table = spc_df(),
                                                               psm_table = psm_df(),
                                                               quantification_csv = MSstats_quant_df(),
                                                               outlier_table = outlier_table)
    })
  })
  # Show/Hide Download button based on the availability of the excel file
  output$reportbuilt <-reactive({
    return(!is.null(exce_report_holder$filepath))
  })
  outputOptions(output, 'reportbuilt', suspendWhenHidden=FALSE)
  # Download excel handling
  output$excel_report_download <- downloadHandler(
    filename = function() {
      paste0(input$excelFilename, ".xlsx")
    },
    content = function(file){
      file.copy(exce_report_holder$filepath, file)
      file.remove(exce_report_holder$filepath)
    }
  )

  ## For Qual
  observeEvent(input$generating_excel_report_qual, {
    withProgress(message = "Building the Excel Report...", value = 5, {
      logFC_colname <- colnames(MSstats_comparison_df())[grepl("FC", colnames(MSstats_comparison_df()))]
      exce_report_holder$filepath <- excel_report_generator_qual(filename = "temp/report.xlsx",
                                                               # UpSet Plot
                                                               peptide_q_value = input$excelPeptideQvalue,
                                                               protein_q_value = input$excelProteinQvalue,
                                                               # data from RData
                                                               experiment_design = qual_experiment_design_df(),
                                                               fasta_df = fasta_df(),
                                                               saint_table = protein_saint_table(),
                                                               spc_table = spc_df(),
                                                               psm_table = psm_df())
    })
  })
  # Show/Hide Download button based on the availability of the excel file
  output$reportbuilt_qual <-reactive({
    return(!is.null(exce_report_holder$filepath))
  })
  outputOptions(output, 'reportbuilt_qual', suspendWhenHidden=FALSE)
  # Download excel handling
  output$excel_report_download_qual <- downloadHandler(
    filename = function() {
      paste0(input$excelFilename_qual, ".xlsx")
    },
    content = function(file){
      file.copy(exce_report_holder$filepath, file)
      file.remove(exce_report_holder$filepath)
    }
  )
  ##### /Excel Report
  

  ####################################################################### Galaxy Upload tool server side code v.0.1.5
  ###### Galaxy upload global
  options(shiny.maxRequestSize=50000*1024^2)
  
  galaxy_address<-'192.168.2.102'
  galaxy_API_key<-reactive({USER$api_key})
  
  fileName= c("example-file-1.raw","example-file-2-A.raw","example-file-2-B.raw",NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
  bioReplicate = c(1, 2, 2,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA) 
  #displayName = c("Example 1", "Example 2A", "Example 2B",NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
  #fractionGroupString = c("example-file-1","example-file-2-A","example-file-2-B",NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
  bioConditionName = c("Control", "YFG-KO", "YFG-KO",NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
  controlBool = c(TRUE, FALSE, FALSE,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA) 
  #DF = data.frame(fileName,bioReplicate,displayName,bioConditionName,controlBool,fractionGroupString)
  DF = data.frame(fileName,bioReplicate,bioConditionName,controlBool)
  DF_DIA = data.frame(fileName,bioReplicate,bioConditionName,controlBool)
  #colnames(DF)<-c("File Name","BioReplicate int","Display Name","Condition","Control","FractionGroup String")
  colnames(DF)<-c("File Name","BioReplicate int","Condition","Control")
  colnames(DF_DIA)<-c("File Name","BioReplicate int","Condition","Control")
  ###### /Galaxy upload global
  
  values <- reactiveValues()
  
  output$diawindowfileReceived <-reactive({
    return(!is.null(input$diawindowfile))
  })
  outputOptions(output, 'diawindowfileReceived', suspendWhenHidden=FALSE)

  output$skylinefileReceived <-reactive({
    return(!is.null(input$skylinefile))
  })
  outputOptions(output, 'skylinefileReceived', suspendWhenHidden=FALSE)
  
  output$fastafileReceived <-reactive({
    return(!is.null(input$fastafile))
  })
  outputOptions(output, 'fastafileReceived', suspendWhenHidden=FALSE)
  
  output$datafilesReceived <-reactive({
    return(!is.null(input$files))
  })
  outputOptions(output, 'datafilesReceived', suspendWhenHidden=FALSE)



  output$skylinefileReceivedDIADDA <-reactive({
    return(!is.null(input$skylinefileDIADDA))
  })
  outputOptions(output, 'skylinefileReceivedDIADDA', suspendWhenHidden=FALSE)
  
  output$fastafileReceivedDIADDA <-reactive({
    return(!is.null(input$fastafileDIADDA))
  })
  outputOptions(output, 'fastafileReceivedDIADDA', suspendWhenHidden=FALSE)
  
  output$datafilesReceivedDIADDA <-reactive({
    return(!is.null(input$filesDIADDA))
  })
  outputOptions(output, 'datafilesReceivedDIADDA', suspendWhenHidden=FALSE)

  output$datafilesReceivedDIADDA_DIA <-reactive({
    return(!is.null(input$filesDIADDA_DIA))
  })
  outputOptions(output, 'datafilesReceivedDIADDA_DIA', suspendWhenHidden=FALSE)
  
  
  #We'll handle the file upload here! (LFQ // QA, just DDA)
  output$uploaded<-renderText({
    inFiles <- input$files
    print(inFiles)
    
    
    if(is.null(inFiles)){
      print("inFiles IS NULL!!!")
      return(NULL)
    }
    else{
      print("inFiles WAS NOT null...")
    }
    if(input$historyName==""){
      print("We will upload the data into the most recent history used...")
    } else {
      history<-input$historyName
      print("We will upload data into history ...")
      print(history)
    }
    
    #Let's write out python script...
    
    python_file_tmp<-file(tempfile(pattern = "file", tmpdir = tempdir()))
    python_file_path<-summary(python_file_tmp)$description
    close(python_file_tmp)
    python_file<-file(python_file_path,"w")
    print(python_file_path)
    print("that's the python file location")
    write('from bioblend.galaxy import GalaxyInstance\nimport sys\nimport os\n',python_file,append=TRUE)
    
    
    
    
    write(paste0("gi = GalaxyInstance(\"",galaxy_address,"\", key=\'",galaxy_API_key(),"\')\n"),python_file,append=TRUE)
    write("history_id=0\n",python_file,append=TRUE)
    write("histories=gi.histories.get_histories()\nhistory_id=0\nfor each_history in histories:\n",python_file,append=TRUE)
    write(paste0("\tif each_history[u\'name\']==\'",history,"\' and not each_history[u\'deleted\']:\n"),python_file,append=TRUE)
    write("\t\thistory_id=each_history[u\'id\']\n\t\tbreak\n",python_file,append=TRUE)
    
    write("try:\n\tif history_id==0:\n",python_file,append=TRUE)
    write(paste0("\t\tnew_hist=gi.histories.create_history(name=\'",history,"\')\n"),python_file,append=TRUE)
    write("\t\thistory_id=new_hist[u\'id\']\n",python_file,append=TRUE)
    write(paste0("except:\n\tif history_id==0:\n\t\thistory_id=histories[0][u\'id\']\n"),python_file,append=TRUE)
    
    write("incoming_files={",python_file,append=TRUE)
    for(i in 1:length(inFiles[,1])){
      write(paste0("\'",inFiles[[i,'datapath']],"\':\'",inFiles$name[i],"\'"),python_file,append=TRUE)
      if(i!=length(inFiles[,1])){
        write(",",python_file,append=TRUE)
      }
      
    }
    write("}\n",python_file,append=TRUE)
    
    
    write("for each_path in incoming_files:\n\ttry:\n\t\tos.symlink(each_path,incoming_files[each_path])\n",python_file,append=TRUE)
    write("\texcept:\n\t\tos.remove(incoming_files[each_path])\n\t\tos.symlink(each_path,incoming_files[each_path])\n",python_file,append=TRUE)
    write("collection_dict={\'collection_type\':\'list\',\'element_identifiers\':[],\'name\':\'RAW Mass Spec Files\'}\n",python_file,append=TRUE)
    close(python_file)
    print("That's the python file!")
    
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    python.load(python_file_path)
    progress$set(message = "Transferring Files to Galaxy...", value = 0)
    if(length(inFiles$name)>=1){
      for(i in 1:length(inFiles[,1])){
        progress$inc(1/length(inFiles[,1]), detail = paste("Uploading file ", i))
        python.exec(paste0("uploaded=gi.tools.upload_file(\"",inFiles$name[i],"\",history_id)"))
        python.exec("collection_dict[\'element_identifiers\'].append({\'id':\'{0}\'.format(uploaded[u\'outputs\'][0][u\'id\']),\'name\':uploaded[u\'outputs\'][0][u\'name\'],\'src\':\'hda\'})")
      }
      python.exec("gi.histories.create_dataset_collection(history_id,collection_dict)")
    }
  })
    #End LFQ //QA file upload


  #We'll handle the file upload here! (DIA+DDA) (DDA PORTION)
  output$uploadedDIADDA<-renderText({
    inFiles <- input$filesDIADDA
    print(inFiles)
    
    
    if(is.null(inFiles)){
      print("inFiles IS NULL!!!")
      return(NULL)
    }
    else{
      print("inFiles WAS NOT null...")
    }
    if(input$historyNameDIADDA==""){
      print("We will upload the data into the most recent history used...")
    } else {
      history<-input$historyNameDIADDA
      print("We will upload data into history ...")
      print(history)
    }
    
    #Let's write out python script...
    
    python_file_tmp<-file(tempfile(pattern = "file", tmpdir = tempdir()))
    python_file_path<-summary(python_file_tmp)$description
    close(python_file_tmp)
    python_file<-file(python_file_path,"w")
    print(python_file_path)
    print("that's the python file location")
    write('from bioblend.galaxy import GalaxyInstance\nimport sys\nimport os\n',python_file,append=TRUE)
    
    
    
    
    write(paste0("gi = GalaxyInstance(\"",galaxy_address,"\", key=\'",galaxy_API_key(),"\')\n"),python_file,append=TRUE)
    write("history_id=0\n",python_file,append=TRUE)
    write("histories=gi.histories.get_histories()\nhistory_id=0\nfor each_history in histories:\n",python_file,append=TRUE)
    write(paste0("\tif each_history[u\'name\']==\'",history,"\' and not each_history[u\'deleted\']:\n"),python_file,append=TRUE)
    write("\t\thistory_id=each_history[u\'id\']\n\t\tbreak\n",python_file,append=TRUE)
    
    write("try:\n\tif history_id==0:\n",python_file,append=TRUE)
    write(paste0("\t\tnew_hist=gi.histories.create_history(name=\'",history,"\')\n"),python_file,append=TRUE)
    write("\t\thistory_id=new_hist[u\'id\']\n",python_file,append=TRUE)
    write(paste0("except:\n\tif history_id==0:\n\t\thistory_id=histories[0][u\'id\']\n"),python_file,append=TRUE)
    
    write("incoming_files={",python_file,append=TRUE)
    for(i in 1:length(inFiles[,1])){
      write(paste0("\'",inFiles[[i,'datapath']],"\':\'",inFiles$name[i],"\'"),python_file,append=TRUE)
      if(i!=length(inFiles[,1])){
        write(",",python_file,append=TRUE)
      }
      
    }
    write("}\n",python_file,append=TRUE)
    
    
    write("for each_path in incoming_files:\n\ttry:\n\t\tos.symlink(each_path,incoming_files[each_path])\n",python_file,append=TRUE)
    write("\texcept:\n\t\tos.remove(incoming_files[each_path])\n\t\tos.symlink(each_path,incoming_files[each_path])\n",python_file,append=TRUE)
    write("collection_dict={\'collection_type\':\'list\',\'element_identifiers\':[],\'name\':\'DDA (ID) RAW Mass Spec Files\'}\n",python_file,append=TRUE)
    close(python_file)
    print("That's the python file!")
    
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    python.load(python_file_path)
    progress$set(message = "Transferring Files to Galaxy...", value = 0)
    if(length(inFiles$name)>=1){
      for(i in 1:length(inFiles[,1])){
        progress$inc(1/length(inFiles[,1]), detail = paste("Uploading file ", i))
        python.exec(paste0("uploaded=gi.tools.upload_file(\"",inFiles$name[i],"\",history_id)"))
        python.exec("collection_dict[\'element_identifiers\'].append({\'id':\'{0}\'.format(uploaded[u\'outputs\'][0][u\'id\']),\'name\':uploaded[u\'outputs\'][0][u\'name\'],\'src\':\'hda\'})")
      }
      python.exec("gi.histories.create_dataset_collection(history_id,collection_dict)")
    }
  })
    #End LFQ DIA+DDA (DDA PORTION)

  #We'll handle the file upload here! (DIA+DDA) (DIA PORTION)
  output$uploadedDIADDA_DIA<-renderText({
    inFiles <- input$filesDIADDA_DIA
    print(inFiles)
    
    
    if(is.null(inFiles)){
      print("inFiles IS NULL!!!")
      return(NULL)
    }
    else{
      print("inFiles WAS NOT null...")
    }
    if(input$historyNameDIADDA==""){
      print("We will upload the data into the most recent history used...")
    } else {
      history<-input$historyNameDIADDA
      print("We will upload data into history ...")
      print(history)
    }
    
    #Let's write out python script...
    
    python_file_tmp<-file(tempfile(pattern = "file", tmpdir = tempdir()))
    python_file_path<-summary(python_file_tmp)$description
    close(python_file_tmp)
    python_file<-file(python_file_path,"w")
    print(python_file_path)
    print("that's the python file location")
    write('from bioblend.galaxy import GalaxyInstance\nimport sys\nimport os\n',python_file,append=TRUE)
    
    
    
    
    write(paste0("gi = GalaxyInstance(\"",galaxy_address,"\", key=\'",galaxy_API_key(),"\')\n"),python_file,append=TRUE)
    write("history_id=0\n",python_file,append=TRUE)
    write("histories=gi.histories.get_histories()\nhistory_id=0\nfor each_history in histories:\n",python_file,append=TRUE)
    write(paste0("\tif each_history[u\'name\']==\'",history,"\' and not each_history[u\'deleted\']:\n"),python_file,append=TRUE)
    write("\t\thistory_id=each_history[u\'id\']\n\t\tbreak\n",python_file,append=TRUE)
    
    write("try:\n\tif history_id==0:\n",python_file,append=TRUE)
    write(paste0("\t\tnew_hist=gi.histories.create_history(name=\'",history,"\')\n"),python_file,append=TRUE)
    write("\t\thistory_id=new_hist[u\'id\']\n",python_file,append=TRUE)
    write(paste0("except:\n\tif history_id==0:\n\t\thistory_id=histories[0][u\'id\']\n"),python_file,append=TRUE)
    
    write("incoming_files={",python_file,append=TRUE)
    for(i in 1:length(inFiles[,1])){
      write(paste0("\'",inFiles[[i,'datapath']],"\':\'",inFiles$name[i],"\'"),python_file,append=TRUE)
      if(i!=length(inFiles[,1])){
        write(",",python_file,append=TRUE)
      }
      
    }
    write("}\n",python_file,append=TRUE)
    
    
    write("for each_path in incoming_files:\n\ttry:\n\t\tos.symlink(each_path,incoming_files[each_path])\n",python_file,append=TRUE)
    write("\texcept:\n\t\tos.remove(incoming_files[each_path])\n\t\tos.symlink(each_path,incoming_files[each_path])\n",python_file,append=TRUE)
    write("collection_dict={\'collection_type\':\'list\',\'element_identifiers\':[],\'name\':\'DIA (Quant) RAW Mass Spec Files\'}\n",python_file,append=TRUE)
    close(python_file)
    print("That's the python file!")
    
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    python.load(python_file_path)
    progress$set(message = "Transferring Files to Galaxy...", value = 0)
    if(length(inFiles$name)>=1){
      for(i in 1:length(inFiles[,1])){
        progress$inc(1/length(inFiles[,1]), detail = paste("Uploading file ", i))
        python.exec(paste0("uploaded=gi.tools.upload_file(\"",inFiles$name[i],"\",history_id)"))
        python.exec("collection_dict[\'element_identifiers\'].append({\'id':\'{0}\'.format(uploaded[u\'outputs\'][0][u\'id\']),\'name\':uploaded[u\'outputs\'][0][u\'name\'],\'src\':\'hda\'})")
      }
      python.exec("gi.histories.create_dataset_collection(history_id,collection_dict)")
    }
  })
    #End LFQ DIA+DDA (DIA PORTION)

  
  ## Handsontable LFQ
  observe({
    if (!is.null(input$hot)) {
      DF = hot_to_r(input$hot)
    } else {
      if (is.null(values[["DF"]]))
        DF <- DF
      else
        DF <- values[["DF"]]
    }
    values[["DF"]] <- DF
  })

  ## Handsontable LFQ DDA+DIA
  #DDA
  observe({
    if (!is.null(input$hotDIADDA)) {
      DF = hot_to_r(input$hotDIADDA)
    } else {
      if (is.null(values[["DF"]]))
        DF <- DF
      else
        DF <- values[["DF"]]
    }
    values[["DF"]] <- DF
  })

  #DIA
  observe({
    if (!is.null(input$hotDIADDA_DIA)) {
      DF_DIA = hot_to_r(input$hotDIADDA_DIA)
    } else {
      if (is.null(values[["DF_DIA"]]))
        DF_DIA <- DF_DIA
      else
        DF_DIA <- values[["DF_DIA"]]
    }
    values[["DF_DIA"]] <- DF_DIA
  })

  
  ## Fill out the handsontable based on the uploaded file names... (LFQ)
  observe({
    if(!is.null(input$files)){
      inFiles<-input$files
      print(inFiles)
      
      fileName= inFiles$name
      fileName<-unlist(strsplit(fileName, ".(?!.*\\.)", perl = TRUE))
      fileName<-fileName[fileName!=""]
      bioReplicate = integer(length(fileName))
      bioConditionName <- rep("",times=length(fileName))
      controlBool <- rep(FALSE,times=length(fileName))
	  DF = data.frame(fileName,bioReplicate,bioConditionName,controlBool,stringsAsFactors = FALSE)
	  colnames(DF)<-c("File Name","BioReplicate int","Condition","Control")
      values[["DF"]]<-DF
    }
  })  #END HANDSONTABLE UPADTE CODE LFQ

  ## Fill out the handsontable based on the uploaded file names... (LFQ DIA+DDA)

  observe({ #DDA
    if(!is.null(input$filesDIADDA)){
      inFiles<-input$filesDIADDA
      print(inFiles)
      
      fileName= inFiles$name
      fileName<-unlist(strsplit(fileName, ".(?!.*\\.)", perl = TRUE))
      fileName<-fileName[fileName!=""]
      bioReplicate = integer(length(fileName))
      bioConditionName <- rep("",times=length(fileName))
      controlBool <- rep(FALSE,times=length(fileName))
	  DF = data.frame(fileName,bioReplicate,bioConditionName,controlBool,stringsAsFactors = FALSE)
	  colnames(DF)<-c("File Name","BioReplicate int","Condition","Control")
      values[["DF"]]<-DF
    }
  })  #END HANDSONTABLE UPADTE CODE LFQ DIA+DDA (DDA)

  
  observe({ #DIA
    if(!is.null(input$filesDIADDA_DIA)){
      inFiles<-input$filesDIADDA_DIA
      print(inFiles)
      
      fileName= inFiles$name
      fileName<-unlist(strsplit(fileName, ".(?!.*\\.)", perl = TRUE))
      fileName<-fileName[fileName!=""]
      bioReplicate = integer(length(fileName))
      bioConditionName <- rep("",times=length(fileName))
      controlBool <- rep(FALSE,times=length(fileName))
	  DF_DIA = data.frame(fileName,bioReplicate,bioConditionName,controlBool,stringsAsFactors = FALSE)
	  colnames(DF_DIA)<-c("File Name","BioReplicate int","Condition","Control")
      values[["DF_DIA"]]<-DF_DIA
    }
  })  #END HANDSONTABLE UPADTE CODE LFQ DIA+DDA (DIA)

  
  #Handle FASTA file upload... LFQ
  output$uploadedFASTA<-renderText({
    inFiles <- input$fastafile
    history<-input$historyName
    if(is.null(inFiles)){
      print("inFiles IS NULL!!!")
      return(NULL)
    }
    else{
      print("inFiles WAS NOT null... Uploading FASTA!")
    }
    
    python_file_tmp<-file(tempfile(pattern = "file", tmpdir = tempdir()))
    python_file_path<-summary(python_file_tmp)$description
    close(python_file_tmp)
    python_file<-file(python_file_path,"w")
    print(python_file_path)
    print("that's the python file location for writing the output....")
    write('from bioblend.galaxy import GalaxyInstance\nimport sys\nimport os\n',python_file,append=TRUE)
    
    write(paste0("gi = GalaxyInstance(\"",galaxy_address,"\", key=\'",galaxy_API_key(),"\')\n"),python_file,append=TRUE)
    write("histories=gi.histories.get_histories()\nhistory_id=0\nfor each_history in histories:\n",python_file,append=TRUE)
    write(paste0("\tif each_history[u\'name\']==\'",history,"\' and not each_history[u\'deleted\']:\n"),python_file,append=TRUE)
    write("\t\thistory_id=each_history[u\'id\']\n\t\tbreak\n",python_file,append=TRUE)
    
    write("if True:\n\tif history_id==0:\n",python_file,append=TRUE)##### <------------------------------------------------------ CHANGED FROM TRY TO IF TRUE
    write(paste0("\t\tnew_hist=gi.histories.create_history(name=\'",history,"\')\n"),python_file,append=TRUE)
    write("\t\thistory_id=new_hist[u\'id\']\n",python_file,append=TRUE)
    
    #Here is the code to add in the tags!
    write(paste0("\t\tgi.histories.create_history_tag(history_id,tag=\"pi:",gsub(" ","_",toupper(input$pifirstName)),"_",gsub(" ","_",toupper(input$pilastName)),"\")\n"),python_file,append=TRUE)
    write(paste0("\t\tgi.histories.create_history_tag(history_id,tag=\"cc:",gsub(" ","_",toupper(input$sampleContactName)),"\")\n"),python_file,append=TRUE)
    
    write("incoming_files={",python_file,append=TRUE)
    for(i in 1:length(inFiles[,1])){
      write(paste0("\'",inFiles[[i,'datapath']],"\':\'",inFiles$name[i],"\'"),python_file,append=TRUE)
      if(i!=length(inFiles[,1])){
        write(",",python_file,append=TRUE)
      }
      
    }
    write("}\n",python_file,append=TRUE)
    
    
    write(paste0("for each_path in incoming_files:\n\ttry:\n\t\tprint each_path,incoming_files[each_path],os.chdir(os.path.dirname(each_path))\n\t\tos.symlink(each_path,incoming_files[each_path])\n"),python_file,append=TRUE)
    
    write("\texcept:\n\t\tos.remove(incoming_files[each_path])\n\t\tos.symlink(each_path,incoming_files[each_path])\n",python_file,append=TRUE)
    
    write("for each_file in incoming_files.values():\n\tuploaded=gi.tools.upload_file(each_file,history_id)\n",python_file,append=TRUE)
    
    close(python_file)
    python.load(python_file_path)
    
    print("We're done uploading the FASTA file.")
    
  }) #END LFQ DDA FASTA UPLOAD


  #Handle FASTA file upload... LFQ
  output$uploadedFASTADIADDA<-renderText({
    inFiles <- input$fastafileDIADDA
    history<-input$historyNameDIADDA
    if(is.null(inFiles)){
      print("inFiles IS NULL!!!")
      return(NULL)
    }
    else{
      print("inFiles WAS NOT null... Uploading FASTA!")
    }
    
    python_file_tmp<-file(tempfile(pattern = "file", tmpdir = tempdir()))
    python_file_path<-summary(python_file_tmp)$description
    close(python_file_tmp)
    python_file<-file(python_file_path,"w")
    print(python_file_path)
    print("that's the python file location for writing the output....")
    write('from bioblend.galaxy import GalaxyInstance\nimport sys\nimport os\n',python_file,append=TRUE)
    
    write(paste0("gi = GalaxyInstance(\"",galaxy_address,"\", key=\'",galaxy_API_key(),"\')\n"),python_file,append=TRUE)
    write("histories=gi.histories.get_histories()\nhistory_id=0\nfor each_history in histories:\n",python_file,append=TRUE)
    write(paste0("\tif each_history[u\'name\']==\'",history,"\' and not each_history[u\'deleted\']:\n"),python_file,append=TRUE)
    write("\t\thistory_id=each_history[u\'id\']\n\t\tbreak\n",python_file,append=TRUE)
    
    write("if True:\n\tif history_id==0:\n",python_file,append=TRUE)##### <------------------------------------------------------ CHANGED FROM TRY TO IF TRUE
    write(paste0("\t\tnew_hist=gi.histories.create_history(name=\'",history,"\')\n"),python_file,append=TRUE)
    write("\t\thistory_id=new_hist[u\'id\']\n",python_file,append=TRUE)
    
    #Here is the code to add in the tags!
    write(paste0("\t\tgi.histories.create_history_tag(history_id,tag=\"pi:",gsub(" ","_",toupper(input$pifirstNameDIADDA)),"_",gsub(" ","_",toupper(input$pilastNameDIADDA)),"\")\n"),python_file,append=TRUE)
    write(paste0("\t\tgi.histories.create_history_tag(history_id,tag=\"cc:",gsub(" ","_",toupper(input$sampleContactNameDIADDA)),"\")\n"),python_file,append=TRUE)
    
    write("incoming_files={",python_file,append=TRUE)
    for(i in 1:length(inFiles[,1])){
      write(paste0("\'",inFiles[[i,'datapath']],"\':\'",inFiles$name[i],"\'"),python_file,append=TRUE)
      if(i!=length(inFiles[,1])){
        write(",",python_file,append=TRUE)
      }
      
    }
    write("}\n",python_file,append=TRUE)
    
    
    write(paste0("for each_path in incoming_files:\n\ttry:\n\t\tprint each_path,incoming_files[each_path],os.chdir(os.path.dirname(each_path))\n\t\tos.symlink(each_path,incoming_files[each_path])\n"),python_file,append=TRUE)
    
    write("\texcept:\n\t\tos.remove(incoming_files[each_path])\n\t\tos.symlink(each_path,incoming_files[each_path])\n",python_file,append=TRUE)
    
    write("for each_file in incoming_files.values():\n\tuploaded=gi.tools.upload_file(each_file,history_id)\n",python_file,append=TRUE)
    
    close(python_file)
    python.load(python_file_path)
    
    print("We're done uploading the FASTA file.")
    
  }) #END LFQ DIA+DDA FASTA UPLOAD


  
  #Handle Skyline file upload... LFQ DDA
  output$uploadedSkyline<-renderText({
    inFiles <- input$skylinefile
    history<-input$historyName
    if(is.null(inFiles)){
      print("inFiles IS NULL!!!")
      return(NULL)
    }
    else{
      print("inFiles WAS NOT null... Uploading Skyline File!")
    }
    
    python_file_tmp<-file(tempfile(pattern = "file", tmpdir = tempdir()))
    python_file_path<-summary(python_file_tmp)$description
    close(python_file_tmp)
    python_file<-file(python_file_path,"w")
    print(python_file_path)
    print("that's the python file location for writing the output....")
    write('from bioblend.galaxy import GalaxyInstance\nimport sys\nimport os\n',python_file,append=TRUE)
    
    write(paste0("gi = GalaxyInstance(\"",galaxy_address,"\", key=\'",galaxy_API_key(),"\')\n"),python_file,append=TRUE)
    #write("history_id=0\n",python_file,append=TRUE)
    write("histories=gi.histories.get_histories()\nhistory_id=0\nfor each_history in histories:\n",python_file,append=TRUE)
    write(paste0("\tif each_history[u\'name\']==\'",history,"\' and not each_history[u\'deleted\']:\n"),python_file,append=TRUE)
    write("\t\thistory_id=each_history[u\'id\']\n\t\tbreak\n",python_file,append=TRUE)
    
    write("try:\n\tif history_id==0:\n",python_file,append=TRUE)
    write(paste0("\t\tnew_hist=gi.histories.create_history(name=\'",history,"\')\n"),python_file,append=TRUE)
    write("\t\thistory_id=new_hist[u\'id\']\n",python_file,append=TRUE)
    write(paste0("except:\n\tif history_id==0:\n\t\thistory_id=histories[0][u\'id\']\n"),python_file,append=TRUE)
    
    
    write("incoming_files={",python_file,append=TRUE)
    for(i in 1:length(inFiles[,1])){
      write(paste0("\'",inFiles[[i,'datapath']],"\':\'",inFiles$name[i],"\'"),python_file,append=TRUE)
      if(i!=length(inFiles[,1])){
        write(",",python_file,append=TRUE)
      }
      
    }
    write("}\n",python_file,append=TRUE)
    
    
    write("for each_path in incoming_files:\n\ttry:\n\t\tos.symlink(each_path,incoming_files[each_path])\n",python_file,append=TRUE)
    write("\texcept:\n\t\tos.remove(incoming_files[each_path])\n\t\tos.symlink(each_path,incoming_files[each_path])\n",python_file,append=TRUE)
    
    write("for each_file in incoming_files.values():\n\tuploaded=gi.tools.upload_file(each_file,history_id)\n",python_file,append=TRUE)
    
    close(python_file)
    python.load(python_file_path)
    
    print("We're done uploading the Skyline file.")
    
  })

  #Handle Skyline file upload... LFQ DIA+DDA
  output$uploadedSkylineDIADDA<-renderText({
    inFiles <- input$skylinefileDIADDA
    history<-input$historyNameDIADDA
    if(is.null(inFiles)){
      print("inFiles IS NULL!!!")
      return(NULL)
    }
    else{
      print("inFiles WAS NOT null... Uploading Skyline File!")
    }
    
    python_file_tmp<-file(tempfile(pattern = "file", tmpdir = tempdir()))
    python_file_path<-summary(python_file_tmp)$description
    close(python_file_tmp)
    python_file<-file(python_file_path,"w")
    print(python_file_path)
    print("that's the python file location for writing the output....")
    write('from bioblend.galaxy import GalaxyInstance\nimport sys\nimport os\n',python_file,append=TRUE)
    
    write(paste0("gi = GalaxyInstance(\"",galaxy_address,"\", key=\'",galaxy_API_key(),"\')\n"),python_file,append=TRUE)
    write("histories=gi.histories.get_histories()\nhistory_id=0\nfor each_history in histories:\n",python_file,append=TRUE)
    write(paste0("\tif each_history[u\'name\']==\'",history,"\' and not each_history[u\'deleted\']:\n"),python_file,append=TRUE)
    write("\t\thistory_id=each_history[u\'id\']\n\t\tbreak\n",python_file,append=TRUE)
    
    write("try:\n\tif history_id==0:\n",python_file,append=TRUE)
    write(paste0("\t\tnew_hist=gi.histories.create_history(name=\'",history,"\')\n"),python_file,append=TRUE)
    write("\t\thistory_id=new_hist[u\'id\']\n",python_file,append=TRUE)
    write(paste0("except:\n\tif history_id==0:\n\t\thistory_id=histories[0][u\'id\']\n"),python_file,append=TRUE)
    
    
    write("incoming_files={",python_file,append=TRUE)
    for(i in 1:length(inFiles[,1])){
      write(paste0("\'",inFiles[[i,'datapath']],"\':\'",inFiles$name[i],"\'"),python_file,append=TRUE)
      if(i!=length(inFiles[,1])){
        write(",",python_file,append=TRUE)
      }
      
    }
    write("}\n",python_file,append=TRUE)
    
    
    write("for each_path in incoming_files:\n\ttry:\n\t\tos.symlink(each_path,incoming_files[each_path])\n",python_file,append=TRUE)
    write("\texcept:\n\t\tos.remove(incoming_files[each_path])\n\t\tos.symlink(each_path,incoming_files[each_path])\n",python_file,append=TRUE)
    
    write("for each_file in incoming_files.values():\n\tuploaded=gi.tools.upload_file(each_file,history_id)\n",python_file,append=TRUE)
    
    close(python_file)
    python.load(python_file_path)
    
    print("We're done uploading the Skyline file.")
    
  }) #END SKYLINE FILE DIA+DDA

  #Handle DIAUmpire Window file upload...
  output$uploadedDIAwindowfile<-renderText({
    inFiles <- input$diawindowfile
    history<-input$historyName
    if(is.null(inFiles)){
      print("inFiles IS NULL!!!")
      return(NULL)
    }
    else{
      print("inFiles WAS NOT null... Uploading DIAUmpire Window File!")
    }
    
    python_file_tmp<-file(tempfile(pattern = "file", tmpdir = tempdir()))
    python_file_path<-summary(python_file_tmp)$description
    close(python_file_tmp)
    python_file<-file(python_file_path,"w")
    print(python_file_path)
    print("that's the python file location for writing the output....")
    write('from bioblend.galaxy import GalaxyInstance\nimport sys\nimport os\n',python_file,append=TRUE)
    
    write(paste0("gi = GalaxyInstance(\"",galaxy_address,"\", key=\'",galaxy_API_key(),"\')\n"),python_file,append=TRUE)
    write("histories=gi.histories.get_histories()\nhistory_id=0\nfor each_history in histories:\n",python_file,append=TRUE)
    write(paste0("\tif each_history[u\'name\']==\'",history,"\' and not each_history[u\'deleted\']:\n"),python_file,append=TRUE)
    write("\t\thistory_id=each_history[u\'id\']\n\t\tbreak\n",python_file,append=TRUE)
    
    write("try:\n\tif history_id==0:\n",python_file,append=TRUE)
    write(paste0("\t\tnew_hist=gi.histories.create_history(name=\'",history,"\')\n"),python_file,append=TRUE)
    write("\t\thistory_id=new_hist[u\'id\']\n",python_file,append=TRUE)
    write(paste0("except:\n\tif history_id==0:\n\t\thistory_id=histories[0][u\'id\']\n"),python_file,append=TRUE)
    
    
    write("incoming_files={",python_file,append=TRUE)
    for(i in 1:length(inFiles[,1])){
      write(paste0("\'",inFiles[[i,'datapath']],"\':\'",inFiles$name[i],"\'"),python_file,append=TRUE)
      if(i!=length(inFiles[,1])){
        write(",",python_file,append=TRUE)
      }
      
    }
    write("}\n",python_file,append=TRUE)
    
    
    write("for each_path in incoming_files:\n\ttry:\n\t\tos.symlink(each_path,incoming_files[each_path])\n",python_file,append=TRUE)
    write("\texcept:\n\t\tos.remove(incoming_files[each_path])\n\t\tos.symlink(each_path,incoming_files[each_path])\n",python_file,append=TRUE)
    
    write("for each_file in incoming_files.values():\n\tuploaded=gi.tools.upload_file(each_file,history_id)\n",python_file,append=TRUE)
    
    close(python_file)
    python.load(python_file_path)
    
    print("We're done uploading the DIA Window file.")
    
  })

  
  output$hot <- renderRHandsontable({
    DF <- values[["DF"]]
    if (!is.null(DF))
      rhandsontable(DF, useTypes = TRUE, stretchH = "all",overflow='visible')
  })

    output$hotDIADDA <- renderRHandsontable({
    DF <- values[["DF"]]
    if (!is.null(DF))
      rhandsontable(DF, useTypes = TRUE, stretchH = "all",overflow='visible')
  })

    output$hotDIADDA_DIA <- renderRHandsontable({
    DF_DIA <- values[["DF_DIA"]]
    if (!is.null(DF_DIA))
      rhandsontable(DF_DIA, useTypes = TRUE, stretchH = "all",overflow='visible')
  })

  output$hot_tmt <- renderRHandsontable({
    DF <- values[["DF"]]
    if (!is.null(DF))
      rhandsontable(DF, useTypes = TRUE, stretchH = "all",overflow='visible')
  })

  
  ## Save LFQ DDA
  observeEvent(input$save, {
    finalDF <- isolate(values[["DF"]])
    history<-input$historyName
    #saveRDS(finalDF, file=file.path(outdir, sprintf("%s.rds", outfilename)))
    design_file_tmp<-file(tempfile(pattern = "Experimental_Design", tmpdir = tempdir()))
    design_file_path<-summary(design_file_tmp)$description
    close(design_file_tmp)
    
    design_file<-file(design_file_path,"w")
    
    for (i in 1:length(finalDF[,'File Name'])){
      if(finalDF[i,'Control']){
        #write(paste(finalDF[i,'File Name'],finalDF[i,'FractionGroup String'],finalDF[i,'Condition'],"C",finalDF[i,'BioReplicate int'],sep="___"),design_file,sep="\n")
		write(paste(finalDF[i,'File Name'],finalDF[i,'File Name'],finalDF[i,'Condition'],"C",finalDF[i,'BioReplicate int'],sep="___"),design_file,sep="\n")
      }else{
        #write(paste(finalDF[i,'File Name'],finalDF[i,'FractionGroup String'],finalDF[i,'Condition'],"T",finalDF[i,'BioReplicate int'],sep="___"),design_file,sep="\n")
		write(paste(finalDF[i,'File Name'],finalDF[i,'File Name'],finalDF[i,'Condition'],"T",finalDF[i,'BioReplicate int'],sep="___"),design_file,sep="\n")
      }
    }
    close(design_file)
    print(design_file_path)
    print("above is where the design file was written...")
    
    python_file_tmp<-file(tempfile(pattern = "file", tmpdir = tempdir()))
    python_file_path<-summary(python_file_tmp)$description
    close(python_file_tmp)
    python_file<-file(python_file_path,"w")
    print(python_file_path)
    print("that's the python file location gor writing the output....")
    write('from bioblend.galaxy import GalaxyInstance\nimport sys\nimport os\n',python_file,append=TRUE)
    
    write(paste0("gi = GalaxyInstance(\"",galaxy_address,"\", key=\'",galaxy_API_key(),"\')\n"),python_file,append=TRUE)
    #write("history_id=0\n",python_file,append=TRUE)
    write("histories=gi.histories.get_histories()\nhistory_id=0\nfor each_history in histories:\n",python_file,append=TRUE)
    write(paste0("\tif each_history[u\'name\']==\'",history,"\' and not each_history[u\'deleted\']:\n"),python_file,append=TRUE)
    write("\t\thistory_id=each_history[u\'id\']\n\t\tbreak\n",python_file,append=TRUE)
    
    write("try:\n\tif history_id==0:\n",python_file,append=TRUE)
    write(paste0("\t\tnew_hist=gi.histories.create_history(name=\'",history,"\')\n"),python_file,append=TRUE)
    write("\t\thistory_id=new_hist[u\'id\']\n",python_file,append=TRUE)
    write(paste0("except:\n\tif history_id==0:\n\t\thistory_id=histories[0][u\'id\']\n"),python_file,append=TRUE)
    
    write(paste0("uploaded=gi.tools.upload_file(\"",design_file_path,"\",history_id)\n"),python_file,append=TRUE)
    
    close(python_file)
    python.load(python_file_path)
    
    showNotification("The table has been uploaded to Galaxy!")
    
  })

  ## Save LFQ DIA+DDA (DDA)
  observeEvent(input$saveDIADDA, {
    finalDF <- isolate(values[["DF"]])
    history<-input$historyNameDIADDA
    #saveRDS(finalDF, file=file.path(outdir, sprintf("%s.rds", outfilename)))
    design_file_tmp<-file(tempfile(pattern = "DDA_Experimental_Design", tmpdir = tempdir()))
    design_file_path<-summary(design_file_tmp)$description
    close(design_file_tmp)
    
    design_file<-file(design_file_path,"w")
    
    for (i in 1:length(finalDF[,'File Name'])){
      if(finalDF[i,'Control']){
        #write(paste(finalDF[i,'File Name'],finalDF[i,'FractionGroup String'],finalDF[i,'Condition'],"C",finalDF[i,'BioReplicate int'],sep="___"),design_file,sep="\n")
		write(paste(finalDF[i,'File Name'],finalDF[i,'File Name'],finalDF[i,'Condition'],"C",finalDF[i,'BioReplicate int'],sep="___"),design_file,sep="\n")
      }else{
        #write(paste(finalDF[i,'File Name'],finalDF[i,'FractionGroup String'],finalDF[i,'Condition'],"T",finalDF[i,'BioReplicate int'],sep="___"),design_file,sep="\n")
		write(paste(finalDF[i,'File Name'],finalDF[i,'File Name'],finalDF[i,'Condition'],"T",finalDF[i,'BioReplicate int'],sep="___"),design_file,sep="\n")
      }
    }
    close(design_file)
    print(design_file_path)
    print("above is where the design file was written...")
    
    python_file_tmp<-file(tempfile(pattern = "file", tmpdir = tempdir()))
    python_file_path<-summary(python_file_tmp)$description
    close(python_file_tmp)
    python_file<-file(python_file_path,"w")
    print(python_file_path)
    print("that's the python file location gor writing the output....")
    write('from bioblend.galaxy import GalaxyInstance\nimport sys\nimport os\n',python_file,append=TRUE)
    
    write(paste0("gi = GalaxyInstance(\"",galaxy_address,"\", key=\'",galaxy_API_key(),"\')\n"),python_file,append=TRUE)
    #write("history_id=0\n",python_file,append=TRUE)
    write("histories=gi.histories.get_histories()\nhistory_id=0\nfor each_history in histories:\n",python_file,append=TRUE)
    write(paste0("\tif each_history[u\'name\']==\'",history,"\' and not each_history[u\'deleted\']:\n"),python_file,append=TRUE)
    write("\t\thistory_id=each_history[u\'id\']\n\t\tbreak\n",python_file,append=TRUE)
    
    write("try:\n\tif history_id==0:\n",python_file,append=TRUE)
    write(paste0("\t\tnew_hist=gi.histories.create_history(name=\'",history,"\')\n"),python_file,append=TRUE)
    write("\t\thistory_id=new_hist[u\'id\']\n",python_file,append=TRUE)
    write(paste0("except:\n\tif history_id==0:\n\t\thistory_id=histories[0][u\'id\']\n"),python_file,append=TRUE)
    
    write(paste0("uploaded=gi.tools.upload_file(\"",design_file_path,"\",history_id)\n"),python_file,append=TRUE)
    
    close(python_file)
    python.load(python_file_path)
    
    showNotification("The DDA table has been uploaded to Galaxy!")
    
  })

  ## Save LFQ DIA+DDA (DIA)
  observeEvent(input$saveDIADDA_DIA, {
    finalDF <- isolate(values[["DF_DIA"]])
    history<-input$historyNameDIADDA
    #saveRDS(finalDF, file=file.path(outdir, sprintf("%s.rds", outfilename)))
    design_file_tmp<-file(tempfile(pattern = "DIA_Experimental_Design", tmpdir = tempdir()))
    design_file_path<-summary(design_file_tmp)$description
    close(design_file_tmp)
    
    design_file<-file(design_file_path,"w")
    
    write(paste("Crux File Integer","Original File Name","Fractionation Group ID String","Fractionation Group Name","Biological Condition","Test or Control","BioReplicate",sep="\t"),design_file,sep="\n")
    for (i in 1:length(finalDF[,'File Name'])){
      if(finalDF[i,'Control']){
        #write(paste(finalDF[i,'File Name'],finalDF[i,'FractionGroup String'],finalDF[i,'Condition'],"C",finalDF[i,'BioReplicate int'],sep="___"),design_file,sep="\n")
        write(paste(toString(i-1),finalDF[i,'File Name'],finalDF[i,'File Name'],finalDF[i,'File Name'],finalDF[i,'Condition'],"C",finalDF[i,'BioReplicate int'],sep="\t"),design_file,sep="\n")

	#write(paste(finalDF[i,'File Name'],finalDF[i,'File Name'],finalDF[i,'Condition'],"C",finalDF[i,'BioReplicate int'],sep="___"),design_file,sep="\n")
      }else{
        #write(paste(finalDF[i,'File Name'],finalDF[i,'FractionGroup String'],finalDF[i,'Condition'],"T",finalDF[i,'BioReplicate int'],sep="___"),design_file,sep="\n")
	#write(paste(finalDF[i,'File Name'],finalDF[i,'File Name'],finalDF[i,'Condition'],"T",finalDF[i,'BioReplicate int'],sep="___"),design_file,sep="\n")
        write(paste(toString(i-1),finalDF[i,'File Name'],finalDF[i,'File Name'],finalDF[i,'File Name'],finalDF[i,'Condition'],"T",finalDF[i,'BioReplicate int'],sep="\t"),design_file,sep="\n")

      }
    }
    close(design_file)
    print(design_file_path)
    print("above is where the design file was written...")
    
    python_file_tmp<-file(tempfile(pattern = "file", tmpdir = tempdir()))
    python_file_path<-summary(python_file_tmp)$description
    close(python_file_tmp)
    python_file<-file(python_file_path,"w")
    print(python_file_path)
    print("that's the python file location gor writing the output....")
    write('from bioblend.galaxy import GalaxyInstance\nimport sys\nimport os\n',python_file,append=TRUE)
    
    write(paste0("gi = GalaxyInstance(\"",galaxy_address,"\", key=\'",galaxy_API_key(),"\')\n"),python_file,append=TRUE)
    #write("history_id=0\n",python_file,append=TRUE)
    write("histories=gi.histories.get_histories()\nhistory_id=0\nfor each_history in histories:\n",python_file,append=TRUE)
    write(paste0("\tif each_history[u\'name\']==\'",history,"\' and not each_history[u\'deleted\']:\n"),python_file,append=TRUE)
    write("\t\thistory_id=each_history[u\'id\']\n\t\tbreak\n",python_file,append=TRUE)
    
    write("try:\n\tif history_id==0:\n",python_file,append=TRUE)
    write(paste0("\t\tnew_hist=gi.histories.create_history(name=\'",history,"\')\n"),python_file,append=TRUE)
    write("\t\thistory_id=new_hist[u\'id\']\n",python_file,append=TRUE)
    write(paste0("except:\n\tif history_id==0:\n\t\thistory_id=histories[0][u\'id\']\n"),python_file,append=TRUE)
    
    write(paste0("uploaded=gi.tools.upload_file(\"",design_file_path,"\",history_id)\n"),python_file,append=TRUE)
    
    close(python_file)
    python.load(python_file_path)
    
    showNotification("The DIA table has been uploaded to Galaxy!")
    
  })

  ####################################################################### /Galaxy Upload tool server side code
}
########### /Shiny SERVER

shinyApp(ui = ui, server = server)
