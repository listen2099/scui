# Building a Prod-Ready, Robust Shiny Application.
# 
# README: each step of the dev files is optional, and you don't have to 
# fill every dev scripts before getting started. 
# 01_start.R should be filled at start. 
# 02_dev.R should be used to keep track of your development during the project.
# 03_deploy.R should be used once you need to deploy your app.
# 
# 
########################################
#### CURRENT FILE: ON START SCRIPT #####
########################################

## Fill the DESCRIPTION ----
## Add meta data about your application
golem::fill_desc(
  pkg_name = "scui", # The Name of the package containing the App 
  pkg_title = "Single-cell analysis user interface", # The Title of the package containing the App 
  pkg_description = "This is a Shiny APP for single cell visual interaction analysis.", # The Description of the package containing the App 
  author_first_name = "Kun", # Your First Name
  author_last_name = "Wei", # Your Last Name
  author_email = "weikun1@genomics.cn", # Your Email
  repo_url = NULL # The URL of the GitHub Repo (optional) 
)     

## Set {golem} options ----
golem::set_golem_options()

## Create Common Files ----
## See ?usethis for more information
usethis::use_mit_license( name = "Wei Kun" )  # You can set another license here
usethis::use_readme_rmd( open = FALSE )
usethis::use_code_of_conduct()
usethis::use_lifecycle_badge( "Experimental" )
usethis::use_news_md( open = FALSE )

## Use git ----
usethis::use_git()

## Init Testing Infrastructure ----
## Create a template for tests
golem::use_recommended_tests()

## Use Recommended Packages ----
golem::use_recommended_deps()
golem::use_recommended_deps(recommended = c("Seurat","shinydashboard","shinycssloaders",
                                            "plotly","shinyWidgets","grDevices",
                                            "DT","dplyr","htmlwidgets","shinyBS",
                                            "heatmaply","htmltools","RColorBrewer","shinyjs")) 
## Favicon ----
# If you want to change the favicon (default is golem's one)
golem::remove_favicon()
golem::use_favicon() # path = "path/to/ico". Can be an online file. 

## Add helper functions ----
golem::use_utils_ui()
golem::use_utils_server()

# You're now set! ----

# go to dev/02_dev.R
rstudioapi::navigateToFile( "dev/02_dev.R" )

# remotes::install_github("satijalab/seurat", ref = "release/4.0.0", lib = .libPaths()[2])
install.packages('RColorBrewer',lib = .libPaths()[2],dependencies = T)
library('RColorBrewer',lib.loc=.libPaths()[2])
remove.packages('reticulate',lib = .libPaths()[2])

BiocManager::install("ComplexHeatmap",lib=.libPaths()[2],dependencies = T)
BiocManager::install("shape",lib=.libPaths()[2],dependencies = T)

devtools::build(path = "C:/myShinyApp/scui")
remove.packages('scui',lib = .libPaths()[length(.libPaths())])
install.packages(
  pkgs = 'C:/myShinyApp/scui/scui_1.0.0.tar.gz',
  lib = .libPaths()[length(.libPaths())],
  repos = NULL, 
  dependencies = T
)
#C:\BGIworkdir\BGIDOC\Single cell cell type study\Datasets\Zheng-PBMCs-2700-hg







