#' Export diet
#'
#' export the MCMC samples of diet for given species
#' @param myfit a fit as returned by \code{\link{fit_escroc}}
#' @param mydata a list as returned by \code{\link{prepare_data}}
#' @param spec_names a character or vector list containing the species names
#' to be plotted
#'
#' @return a list matrices of the same length as spec_names
#' @importFrom coda varnames
#' @examples
#' #importing data
#' data(signature_data)
#' data(prior_diet_matrix)
#' data(LOQ)
#' data(prior_signature_data)
#' prior_delta <- data.frame(tracer=c("X15N","X13C"),mean=c(3,0),sd=c(1,1))
#'
#' #check that everything is ok
#' mydata <- prepare_data(prior_diet_matrix,signature_data,prior_signature_data,isLeftCensored,
#' LOQ,prior_delta)
#'
#' #build the model
#' mymodel <- building_model(mydata)
#'
#' #fit the model
#' myresults <- fit_escroc(mydata, mymodel)
#'
#'
#' head(export_diet_species(myresults,mydata,"Sole"))
#' @export
export_diet_species <- function(myfit, mydata, spec_names) {
  lapply(spec_names, function(spec_name) {
    if (is.na(match(spec_name, rownames(mydata$prey_id))))
      stop("the species name does not match species from diet_matrix")
    column_names <- paste("diet\\[", spec_name, sep = "")
    subset_fit <-
      as.matrix(myfit[, grep(column_names, varnames(myfit)), drop = FALSE])
    colnames(subset_fit) <-
      sapply(colnames(subset_fit), function(old_name) {
        strsplit(substr(old_name, 1, nchar(old_name) - 1), c("[\\[,]"))[[1]][3]
      })
    return(subset_fit)
  })

}
