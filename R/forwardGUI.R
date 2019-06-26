#' Shiny interface for population model
#'
#' @export
forwardGUI <- function()
{
    shiny::runApp(paste0(system.file("extdata","shiny",package="holoSimCell")))
}

