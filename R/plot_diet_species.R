#' Plot diet
#'
#' Plot the diet composition of asked species
#' @param mydata a list as returned by function \code{\link{prepare_data}}
#' @param myres a fit as returned by \code{\link{fit_escroc}}
#' @param spec_names a character or vector list containing the species names
#' to be plotted
#'
#' @return a list of the same length as spec_names with a \code{\link[ggplot2]{ggplot}}
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 xlim
#' @importFrom ggplot2 geom_density
#' @importFrom reshape2 melt
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
#' mydata <- prepare_data(prior_diet_matrix,signature_data,
#' LOQ,prior_signature_data,prior_delta)
#'
#' #build the model
#' mymodel <- building_model(mydata)
#'
#' #fit the model
#' myresults <- fit_escroc(mydata, mymodel)
#'
#' #plot diet
#' print(plot_diet_species(mydata, myresults, "Sole"))
#' @export
plot_diet_species <- function(mydata, myres, spec_names) {
  lapply(spec_names, function(spec_name) {
    if (is.na(match(spec_name, rownames(mydata$prey_id))))
      stop("the species name does not match species from diet_matrix")
    column_names <- paste("diet\\[", spec_name, sep = "")
    subset_fit <-
      as.matrix(myres[, grep(column_names, varnames(myres)), drop = FALSE])
    colnames(subset_fit) <-
      sapply(colnames(subset_fit), function(old_name) {
        strsplit(substr(old_name, 1, nchar(old_name) - 1), c("[\\[,]"))[[1]][3]
      })
    subset_fit_long <- melt(subset_fit)
    names(subset_fit_long) = c("Id", "Species", "Proportion")
    ggplot(subset_fit_long,
           aes_string("Proportion", fill = "Species", colour = "Species")) +
      geom_density(alpha = 0.5,aes_string(y="..scaled..")) +
      xlim(0, 1)
  })

}
