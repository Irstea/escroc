#' Fit the model
#' This is the core funtion of the packages. It fits the model ESCROC using \code{\link[runjags]{run.jags}}
#' and returns an \code{\link[coda]{mcmc.list}} object. ESCROC was described in Ballutaud et al. (in press).
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
#' @importFrom runjags run.jags
#' @importFrom coda gelman.diag
#' @importFrom coda varnames
#' @examples
#' #importing data
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
#' @export
fit_escroc <- function(mydata,mymodel,burnin=1000,sample=1000,adapt=1000,method="parallel",...){
  myinits <- generate_init(mydata)
  res<-run.jags(model=mymodel,monitor = c("random_effect", "delta","mean_signature","diet"),data = mydata,n.chains = 3,adapt = adapt,burnin = burnin,sample=sample,
               summarise = FALSE,method = method,inits=myinits,...)
  myfit <- reformat_results(res,mydata)
  gelman.diag(myfit,multivariate=FALSE)
  return(myfit)

}
