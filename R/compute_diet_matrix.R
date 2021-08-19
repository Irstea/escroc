#' Computes the posterior diet matrix
#'
#' Gives the quantiles of posterior diet
#'
#' @param mydata a list as returned by function \code{\link{prepare_data}}
#' @param myres a fit as returned by \code{\link{fit_escroc}}
#' @param quant a vector of quantiles
#'
#' @return a  3D array of dim length(quant) x nb species x nb species
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
#' #diet
#' compute_diet_matrix(mydata, myresults)
#' @export


compute_diet_matrix <- function(mydata, myres, quant=c(.025,.5,.975)){
  res_matrix <- as.matrix(myres)
  comp_names <- rownames(mydata$alpha_diet)
  sol <- array(0, dim = c(length(quant),
                          length(comp_names),
                          length(comp_names)),
               dimnames = list(paste(quant*100,"%",sep=""),
                               comp_names,
                               comp_names))
  id_flow <- expand.grid(comp_names,comp_names)
  tmp <- mapply(function(i,j){
    diet_name <- paste("diet[", i, ",", j, "]", sep = "")
    if (diet_name %in% colnames(res_matrix)){
      return(c(i,
               j,
               quantile(res_matrix[, c(diet_name)],quant)))
    } else {
      return(c(i, j, rep(0, 3)))
    }
  }, id_flow$Var1, id_flow$Var2)
  for (q in seq_len(length(quant)))
    sol[q, , ] <- matrix(tmp[q+2, ],
                         length(comp_names),
                         length(comp_names),
                         byrow=FALSE)
  return(sol)
}
