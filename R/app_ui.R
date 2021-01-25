#' The application User-Interface
#' 
#' @param request Internal parameter for `{shiny}`. 
#'     DO NOT REMOVE.
#' @import shiny
#' @import Seurat
#' @import shinythemes
#' @import shinydashboard
#' @import shinycssloaders
#' @import plotly
#' @import ggplot2
#' @import shinyWidgets
#' @import grDevices
#' @import DT
#' @import dplyr
#' @import htmlwidgets
#' @import shinyBS
#' @import heatmaply
#' @import htmltools
#' @noRd
app_ui <- function(request) {
  dashboardPage(
    dashboardHeader(title = "Single Cell UI",uiOutput('saveprojectbutton')),
    dashboardSidebar(
      sidebarMenu(
        id = "tabs",
        menuItem("Import Seurat Object", tabName = "ImportSeuratObject", icon = icon("folder-open", lib = "glyphicon")),
        menuItemOutput("CellBrowserUI"),
        # menuItem("Cell Browser", tabName = "CellBrowser", icon = icon("eye-open", lib = "glyphicon")),
        menuItemOutput("FeaturesUI"),
        # menuItem("Features", tabName = "Features" ,icon = icon("list-alt", lib = "glyphicon")),
        menuItemOutput("HeatmapUI"),
        # menuItem("Heat Map", tabName = "Heatmap", icon = icon("th", lib = "glyphicon")),
        menuItemOutput("ViolinUI"),
        # menuItem("Violin", tabName = "Violin", icon = icon("equalizer", lib = "glyphicon")),
        menuItemOutput("DotPlotUI"),
        # menuItem("Dot Plot", tabName = "DotPlot", icon = icon("braille", lib = "font-awesome")),
        menuItemOutput("Re_clusteringUI"),
        # menuItem("Re-clustering", tabName = "Re-clustering", icon = icon("laptop-code", lib = "font-awesome")),
        menuItem("About", tabName = "About", icon = icon("exclamation-sign", lib = "glyphicon")),
        menuItemOutput("plotmodel")
        )

    ),
    dashboardBody(
      tags$head(tags$style(HTML('
                          /* logo */
                          .skin-blue .main-header .logo {
                          background-color: rgb(21,88,160); color: rgb(255,255,255);
                          font-weight: normal;font-size: 24px;text-align: center;
                          }

                          /* logo when hovered */
                          .skin-blue .main-header .logo:hover {
                          background-color: rgb(21,88,160);
                          }

                          /* navbar (rest of the header) */
                          .skin-blue .main-header .navbar {
                          background-color: rgb(21,88,160);
                          }

                          /* main sidebar */
                          .skin-blue .main-sidebar {
                          background-color: rgb(107,174,214);;
                          }

                          # /* main body */
                          # .skin-blue .main-body {
                          # background-color: rgb(0,144,197);
                          # }

                          /* active selected tab in the sidebarmenu */
                          .skin-blue .main-sidebar .sidebar .sidebar-menu .active a{
                          background-color: rgb(8,48,107);
                          color: rgb(255,255,255);font-weight: bold;font-size: 18px;
                          }

                          /* other links in the sidebarmenu */
                          .skin-blue .main-sidebar .sidebar .sidebar-menu a{
                          background-color: rgb(107,174,214);
                          color: rgb(255,255,255);font-weight: bold;
                          }

                          /* other links in the sidebarmenu when hovered */
                          .skin-blue .main-sidebar .sidebar .sidebar-menu a:hover{
                          background-color: rgb(232,245,251);color: rgb(0,144,197);font-weight: bold;
                          }

                          /* toggle button color  */
                          .skin-blue .main-header .navbar .sidebar-toggle{
                          background-color: rgb(21,88,160);color:rgb(255,255,255);
                          }

                          /* toggle button when hovered  */
                          .skin-blue .main-header .navbar .sidebar-toggle:hover{
                          background-color: rgb(21,88,160);color:rgb(255,255,255);
                          }
                          .wrapper {height: auto !important; position:relative; overflow-x:hidden; overflow-y:hidden}
#                           '))),
      tabItems(
        # ------------------------------------------------------------- #
        tabItem(
          tabName = "ImportSeuratObject",
          # h3("Get started:May include some guide info"),
          # h3(" "),
          # h3(" "),
          includeMarkdown("./mainpage.md"),
          h2("Now, Let't Import Seurat Object Here:"),
          fluidRow(
            mod_import_seurat_obj_ui("SeuratRdataFile", "User data (.Rdata format)")
          ),
          fluidRow(
            shinycssloaders::withSpinner(
              verbatimTextOutput("obj_print")
            )
          )
        ),
        # ------------------------------------------------------------- #
        tabItem(
          tabName = "CellBrowser",
          uiOutput('intermodel'),
          uiOutput('editnamemodel'),
          fluidRow(
            column(width = 10,
                   shinycssloaders::withSpinner(
                     plotlyOutput("cellmap",height = "800px")
                   ),
                   style ='padding-left:1px;padding-top:1px'
            ),
            column(width = 2,
                   uiOutput('genelistpickmodel'),
                   uiOutput('addgenelisenamemodel'),
                   fluidRow(
                     uiOutput('inSelect'),
                     style ='padding-left:1px; padding-right:1px; padding-top:1px; padding-bottom:1px'
                   ),
                   fluidRow(
                     DTOutput('selected_features_0'),
                     style ='height:490px; overflow-y: scroll; padding-left:1px; padding-right:1px; padding-top:1px; padding-bottom:1px;background-color:#E8F5FB'
                   ),
                   uiOutput("ifsplitui"),
                   fluidRow(
                     uiOutput('selectsplitcategory'),
                     style ='background-color:#E8F5FB'
                   ),
                   fluidRow(
                     column(width = 6,
                            actionButton(
                              inputId = "clean_my_gene_list",
                              label = "Delete",
                              style = "unite",
                              width = "100%"
                            ),
                            style ='padding-left:0px; padding-right:0px; padding-top:0px; padding-bottom:0px'
                     ),
                     column(width = 6,
                            downloadButton(
                              outputId = "output_my_gene_list",
                              label = "Output",
                              style = "unite",
                              width = "100%"
                            ),
                            style ='padding-left:0px; padding-right:0px; padding-top:0px; padding-bottom:0px'
                     ),
                     style ='padding-left:1px; padding-right:1px; padding-top:1px; padding-bottom:1px; background-color:#E8F5FB'
                   ),
                   style ='padding-left:1px; padding-right:1px; padding-top:1px; padding-bottom:1px'
                   ),
            style ='padding-left:1px; padding-right:1px; padding-top:1px; padding-bottom:1px'
          ),
          shinycssloaders::withSpinner(
            uiOutput('pick_category',inline = T)
          ),
          shinycssloaders::withSpinner(
            uiOutput('pick_category_all',inline = T)
          ),
          fluidRow(
            column(width = 10,
                   uiOutput('pick_cluster')
                   ),
            column(width = 2,
                   uiOutput('selectcolorforpoints')
                   )
          ),
          uiOutput('output_cell_table_buttom'),
          # shinycssloaders::withSpinner(
          #   uiOutput('my_category_edit')
          # ),
          shinycssloaders::withSpinner(
            verbatimTextOutput("see_click_cluster")
          )
          # shinycssloaders::withSpinner(
          #   verbatimTextOutput("monitor2")
          # )
        ),
        # ------------------------------------------------------------- #
        tabItem(
          tabName = "Features",
          uiOutput('ifusemycategory'),
          uiOutput('categoryforgetmarker'),
          shinycssloaders::withSpinner(
            DTOutput("features_datatable")
          ),
          verbatimTextOutput("monitor3"),
          fluidRow(
            column(width = 4,
              actionButton(
                inputId = "add_list_switch",
                label = 'Add To "My gene list"',
                style = "unite",
                width = "100%",
                style = "color: white; background-color: #1f78b4;"
              )
            ),
            column(width = 4,
              actionButton(
                inputId = "clean_list_switch",
                label = "Clean!",
                style = "unite",
                width = "100%",
                style = "color: white; background-color: #e31a1c;"
              )
            ),
            column(width = 4,
                   downloadButton(
                     outputId = "output_feature_table",
                     label = "Output",
                     style = "unite",
                     width = "100%"
                   ))
          )
        ),
        # ------------------------------------------------------------- #
        tabItem(
          tabName = "Heatmap",
          h1("Heatmap"),
          fluidRow(
            column(width = 6, uiOutput('selecttopormygenelist')),
            column(width = 6, uiOutput('topnforheatmap'))
          ),
          fluidRow(
            column(width = 6, uiOutput('categoryforheatmap')),
            column(width = 6, uiOutput('clusterforheatmap'))
          ),
          fluidRow(
            column(width = 6, uiOutput('checkboxforheatmapcluster')),
            column(width = 6, uiOutput('pickcolorforheatmap'))
          ),
          uiOutput('plotmodel_stat'),
          useShinyjs(),
          actionBttn(
            inputId = "heat_map_run",
            label = "Run!",
            style = "unite",
            color = "primary"
          ),
          uiOutput('all_heat_map_ui')
        ),
        # ------------------------------------------------------------- #
        tabItem(
          tabName = "Violin",
          h1("Violin"),
          fluidRow(
            column(width = 4, uiOutput('selecttopormygenelist_for_violin')),
            column(width = 4, uiOutput('selectviolinfeature')),
            column(width = 4, uiOutput('split_sample_check_box'),
                   style ='padding-left:1px; padding-top:25px')
          ),
          fluidRow(
            column(width = 6, uiOutput('categoryforviolin')),
            column(width = 6, uiOutput('clusterforviolin'))
          ),
          uiOutput('plotmodel_stat_v'),
          useShinyjs(),
          actionBttn(
            inputId = "violin_map_run",
            label = "Run!",
            style = "unite",
            color = "primary"
          ),
          uiOutput('all_violin_plot_ui')
        ),
        #------------------------------------------------------------------- #
        tabItem(
          tabName = "DotPlot",
          h1("DotPlot"),
          fluidRow(
            column(width = 4, uiOutput('selecttopormygenelist_for_dotplot')),
            column(width = 4, uiOutput('selectdotplotfeature')),
            column(width = 4, uiOutput('split_sample_check_box_for_dotplot'),
                   style ='padding-left:1px; padding-top:25px')
          ),
          fluidRow(
            column(width = 4, uiOutput('categoryfordotplot')),
            column(width = 6, uiOutput('clusterfordotplot')),
            column(width = 2, uiOutput('pickcolorfordotplot'))
          ),
          uiOutput('plotmodel_dotplot'),
          actionBttn(
            inputId = "dot_plot_run",
            label = "Run!",
            style = "unite",
            color = "primary"
          ),
          uiOutput('all_dot_plot_ui')
        ),
        #------------------------------------------------------------------- #
        tabItem(
          tabName = "Re_clustering",
          h1("Re-clustering"),
          fluidRow(
            column(width = 3,uiOutput("pick_category_for_recluster")),
            column(width = 9,uiOutput("def_category_msg"))
          ),
          actionBttn(
            inputId = "recluster_run",
            label = "Run!",
            style = "unite",
            color = "primary"
          ),
          uiOutput('running_msg'),
          shinycssloaders::withSpinner(
            verbatimTextOutput('reclustermsg')
          ),
          shinycssloaders::withSpinner(
            uiOutput('recluster_result')
          )
        ),
        # 解决重新聚类的问题
        # 用yong户自定义的细胞，重新聚类
        # 取子集
        # 点run
        # 保存新的对象
        #------------------------------------------------------------------- #
        tabItem(
          tabName = "About",
          includeMarkdown("./about.md")
        )
      )
    )
  )
}

#' Add external Resources to the Application
#' 
#' This function is internally used to add external 
#' resources inside the Shiny application. 
#' 
#' @import shiny
#' @importFrom golem add_resource_path activate_js favicon bundle_resources
#' @noRd
golem_add_external_resources <- function(){
  
  add_resource_path(
    'www', app_sys('app/www')
  )
 
  tags$head(
    favicon(),
    bundle_resources(
      path = app_sys('app/www'),
      app_title = 'scui'
    )
    # Add here other external resources
    # for example, you can add shinyalert::useShinyalert() 
  )
}

