#' plot the full diet matrix
#'
#'
#' Plot the full diet matrix. Each line corresponds to a predator and its column to the prey.
#' The color of the circle corresponds to the quantile 50% of the contribution of the prey
#' to the diet of the predator. The radius of the larger circle (respectively smaller circle)
#' corresponds to the quantile 97.5% (respectively 2.5%) of the contribution of the prey
#' to the diet of the predator
#'
#' @param myfit a fit as returned by \code{\link{fit_escroc}}
#' @param mydata a list as returned by \code{\link{prepare_data}}
#'
#' @return nothing
#' @importFrom stats quantile
#' @importFrom corrplot corrplot
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
#' #plot the diet matrix
#' plot_full_diet_matrix(myresults,mydata)
#' @export
plot_full_diet_matrix <- function(myfit, mydata) {
  probs = c(0.025, .5, .975)
  species_name <- rownames(mydata$prey_id)
  q2.5 <-
    q50 <-
    q97.5 <-
    matrix(
      0,
      nrow = length(species_name),
      ncol = length(species_name),
      dimnames = list(species_name, species_name)
    )
  subset_fit <-
    as.matrix(myfit[, grep("diet", varnames(myfit)), drop = FALSE])
  for (pred in species_name) {
    for (prey in species_name) {
      pos <-
        match(paste("diet[", pred, ",", prey, "]", sep = ""),
              colnames(subset_fit))
      if (!is.na(pos)) {
        quant <- quantile(subset_fit[, pos], probs = probs)
        q2.5[pred, prey] <- quant[1]
        q50[pred, prey] <- quant[2]
        q97.5[pred, prey] <- quant[3]
      }

    }
  }
  corrplot(
    q50,
    lowCI.mat	= q2.5,
    uppCI.mat = q97.5,
    method = "circle",
    plotCI = "circle",
    is.corr = FALSE,
    tl.col = "black",
    cl.lim = c(0.000000, 1)
  )
}
