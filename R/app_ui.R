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
#' @import shinyjs
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
          HTML('<h2 id="get-started">Get started</h2>
<p>SCUI is a user-friendly lightweight desktop application for browsing and exploring analysis results of single-cell RNA-seq data generated by Seurat. Biologists can easily get highly customized results from the mouse-interactive UI, instead of any command-line environment.</p>
<h2 id="release-of-scui-v1-0-0">Release of SCUI v1.0.0</h2>
<p>The first version of the SCUI has the following features:</p>
<ul>
<li><strong>CellBrowser</strong> Show the coordinate points of each cell by different categories; Customize the list of genes of interest; Customize categories and clusters; Output figures and data table.</li>
<li><strong>Features</strong> Explore and filter the marker gene table obtained by Seurat; Generate new marker gene table for custom clusters.</li>
<li><strong>Heatmap</strong> Generate a heatmap of the average expression of customized genes in clusters.</li>
<li><strong>Violin</strong> Generate the violin plot of genes of interest in clusters.</li>
<li><strong>DotPlot</strong> Generate the DotPlot result in Seurat.</li>
<li><strong>Re-clustering</strong> The cells in the customized categories can be reclustered; Output a new seurat object with re-clustering results.</li>
</ul>
<p><strong>NOTE:</strong> SCUI requires 4GB ram as a minimum, it may be higher depending on the size of your data. The time elapsed for analysis was influenced by the number of cells.</p>
'),
          h2("Now, Let's Import Seurat Object Here:"),
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
                   # uiOutput('selectcolorforpoints'),
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
            column(width = 12,
                   uiOutput('pick_cluster')
                   )
          ),
          fluidRow(
            column(width = 12,
                   uiOutput('output_cell_table_buttom')
                   )
          )
          # shinycssloaders::withSpinner(
          #   uiOutput('my_category_edit')
          # ),
          # ,shinycssloaders::withSpinner(
          #   verbatimTextOutput("see_click_cluster")
          # )
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
        #------------------------------------------------------------------- #
        tabItem(
          tabName = "About",
          HTML('<h2 id="about">About</h2>
<p>SCUI is an open source software based on MIT license.  We do not guarantee that there will be no problems during the use of the software, but we are happy to accept bug reports and valuable suggestions.</p>
<h2 id="contact-us">Contact Us</h2>
<p>Github issues: https://github.com/listen2099/scui/issues</p>
<p>Email: scui.bug@gmail.com</p>
<h2 id="mit-license">MIT License</h2>
<p>Copyright (c) 2021 Wei Kun</p>
<p>Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the &quot;Software&quot;), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:</p>
<p>The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.</p>
<p>THE SOFTWARE IS PROVIDED &quot;AS IS&quot;, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.</p>')
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

