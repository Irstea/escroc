#' Computes the posterior the flow matrix
#'
#' Gives the quantiles of posterior flow matrix of an ecopath model
#'
#' @param mydata a list as returned by function \code{\link{prepare_data}}
#' @param myres a fit as returned by \code{\link{fit_escroc}}
#' @param quant a vector of quantiles
#'
#' @return a  3D array of dim length(quant) x nb species x nb species
#' @export

compute_flow_matrix <- function(mydata, myres, quant=c(.025,.5,.975)){
  if (is.null(mydata$obs_biomass))
    stop("should be an escropath model")
  res_matrix <- as.matrix(myres)
  nb_species <- mydata$nb_species
  comp_names <- rownames(mydata$alpha_diet)
  sol <- array(0, dim = c(length(quant),
                          nb_species + 2,
                          nb_species + 2),
               dimnames = list(paste(quant*100,"%",sep=""),
                               comp_names,
                               comp_names))
  id_flow <- expand.grid(comp_names,comp_names)
  tmp <- mapply(function(i,j){
    diet_name <- paste("diet[", i, ",", j, "]", sep = "")
    if (diet_name %in% colnames(res_matrix)){
      bio_name <- paste("biomass[", i, "]", sep = "")
      conso_name <- paste("consumption_rate[", i, "]", sep = "")
      return(c(i,
               j,
               quantile(
                 apply(res_matrix[, c(diet_name, bio_name, conso_name)],
                       1,
                       prod),
                 quant)))
    } else {
      return(c(i, j, rep(0, 3)))
    }
  }, id_flow$Var1, id_flow$Var2)
  for (q in seq_len(length(quant)))
    sol[q, , ] <- matrix(tmp[q+2, ],
                         nb_species+2,
                         nb_species+2,byrow=FALSE)
  return(sol)
}
