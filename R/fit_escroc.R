#' Fit the model
#'
#' This is the core funtion of the packages. It fits a model ESCROC
#' or escropath using \code{\link[runjags]{run.jags}}
#' and returns an \code{\link[coda]{mcmc.list}} object.
#' ESCROC was described in Ballutaud et al. (2019).
#'
#' @param mydata a list as returned by function \code{\link{prepare_data}}
#' @param mymodel a string describing the model as build by \code{\link{building_model}}
#' @param burnin number of burnin iterations
#' @param sample number of iterations that will be saved for inference
#' @param adapt number of iterations for adaptive phase, if not sufficient a warning
#' will be thrown
#' @param method methods to parallize the MCMC chains, see \code{\link[runjags]{run.jags}}
#' for further details
#' @param ... additionnal arguments that will be sent to \code{\link[runjags]{run.jags}}
#'
#' @return a \code{\link[coda]{mcmc.list}} storing the 3 MCMC
#' @references Ballutaud, M., Drouineau, H., Carassou, L., Munoz, G., Chevillot,
#'  X., Labadie, P., Budzinski, H., Lobry, J., 2019. EStimating Contaminants
#'  tRansfers Over Complex food webs (ESCROC): An innovative Bayesian method for
#'   estimating POP’s biomagnification in aquatic food webs. Science of The
#'   Total Environment 658, 638–649.
#'   https://doi.org/10.1016/j.scitotenv.2018.12.058
#' @importFrom runjags run.jags
#' @importFrom coda gelman.diag
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
#' #a summary of the results
#' library(coda)
#' summary(myresults)$qua
#' @export
fit_escroc <- function(mydata,mymodel,burnin=1000,sample=1000,adapt=1000,method="parallel",...){
  escropath <- ! is.null(mydata$obs_biomass)
  if (escropath) {
    monitor <- c("random_effect", "delta","mean_signature",
                            "diet_short", "input_Det", "A",
                            "biomass", "trophic_efficiency", "productivity",
                            "uq", "consumption_rate","export_Det")
  } else {
    monitor <- c("random_effect", "delta","mean_signature","diet")
  }
  myinits <- generate_init(mydata)
  res<-run.jags(model=mymodel,
                monitor = monitor,
                data = mydata,
                n.chains = 3,
                adapt = adapt,
                burnin = burnin,
                sample=sample,
                summarise = FALSE,
                method = method,
                inits=myinits,...)
  myfit <- reformat_results(res,mydata)
  print(gelman.diag(myfit,multivariate=FALSE))
  return(myfit)

}
