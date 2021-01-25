#' import_seurat_obj UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_import_seurat_obj_ui <- function(id, label = ".Rdata file"){
  ns <- NS(id)
  tagList(
    fileInput(ns("file"),
              label,
              accept = ".Rdata",
              width = "100%")
  )
}
    
#' import_seurat_obj Server Function
#'
#' @noRd 
mod_import_seurat_obj_server <- function(id){
  moduleServer(
    id,
    function(input, output, session){
      userFile <- reactive({
        validate(need(input$file, message = FALSE))
        input$file
      })
      obj_list <- reactive({
        load(userFile()$datapath)
        read_env <- ls(environment())
        if('savevalues' %in% read_env){
          if( 'SCT' %in% Assays(seurat_obj)){
            DefaultAssay(seurat_obj) <- 'SCT'
          }else{
            DefaultAssay(seurat_obj) <- 'RNA'
          }
          res <- list(
            seurat_obj = seurat_obj,
            seurat_obj_markers = seurat_obj_markers,
            seurat_obj_name = seurat_obj_name,
            savevalues = savevalues
          )
          return(res)
        }else{
          if('seurat_obj' %in% read_env){
            for (obj in read_env) {
              if(class(environment()[[obj]]) == "Seurat"){
                seurat_obj <- environment()[[obj]]
                seurat_obj_name <- obj
              }
            }
            if( paste0(seurat_obj_name,'.markers') %in% read_env ){
              seurat_obj_markers <- environment()[[ paste0(seurat_obj_name,'.markers') ]]
              if( 'SCT' %in% Assays(seurat_obj)){
                DefaultAssay(seurat_obj) <- 'SCT'
              }else{
                DefaultAssay(seurat_obj) <- 'RNA'
              }
              res <- list(
                seurat_obj = seurat_obj,
                seurat_obj_markers = seurat_obj_markers,
                seurat_obj_name = seurat_obj_name
              )
              return(res)
            }
          }else{
            seurat_obj <- NULL
            seurat_obj_name <- NULL
            for (obj in read_env) {
              if(class(environment()[[obj]]) == "Seurat"){
                seurat_obj <- environment()[[obj]]
                seurat_obj_name <- obj
              }
            }
            if( paste0(seurat_obj_name,'.markers') %in% read_env ){
              seurat_obj_markers <- environment()[[ paste0(seurat_obj_name,'.markers') ]]
              if( 'SCT' %in% Assays(seurat_obj)){
                DefaultAssay(seurat_obj) <- 'SCT'
              }else{
                DefaultAssay(seurat_obj) <- 'RNA'
              }
              res <- list(
                seurat_obj = seurat_obj,
                seurat_obj_markers = seurat_obj_markers,
                seurat_obj_name = seurat_obj_name
              )
              return(res)
            }
          }
        }
      })
      observe({
        msg <- sprintf("File %s was loaded", userFile()$name)
        cat(msg, "\n")
      })
      return(obj_list)
    }
  )
}
    
## To be copied in the UI
# mod_import_seurat_obj_ui("import_seurat_obj_ui_1")
    
## To be copied in the server
# callModule(mod_import_seurat_obj_server, "import_seurat_obj_ui_1")
 
