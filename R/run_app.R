#' Run the Shiny Application
#'
#' @param ... A series of options to be used inside the app.
#'
#' @export
#' @importFrom shiny shinyApp
#' @importFrom golem with_golem_options
run_app <- function(
  ...
) {
  options(shiny.maxRequestSize=5*1024^3)
  options(width=150)
  with_golem_options(
    app = shinyApp(
      ui = app_ui, 
      server = app_server,
      options = list(port = 18920)
    ), 
    golem_opts = list(...)
  )
}
