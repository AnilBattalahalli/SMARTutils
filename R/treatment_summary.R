#' @title Function to get treatment summary from the data generator
#' @param recipe the generated recipe from nocovairiates.from_conditional()
#' @export treatment_summary

treatment_summary <- function(recipe){
  return(recipe$treat_summary)
}
