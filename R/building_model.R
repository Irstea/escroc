#' Build the jags model code
#'
#' Based on data, this function build the jags model that will be used by runjages
#'
#' @param mydata a list data returned by \code{\link{prepare_data}}
#'
#' @return a string corresponding to the jags model
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
#' @export
building_model <- function(mydata) {
  priors_noise_regression <-
    "for (i in 1:nb_tracer){tau[i]~dgamma(0.01,0.01)}"

  priors_delta_no_prior <-
    "for (i in 1:(nb_tracer-nb_prior_delta)){
  delta_std[id_no_prior_delta[i]]~dnorm(0,.1)
  delta[id_no_prior_delta[i]]<-delta_std[id_no_prior_delta[i]]*sd_tracer[id_no_prior_delta[i]]
}"

  priors_delta_with_prior <-
    "for (i in 1:(nb_prior_delta)){
  delta[id_prior_delta[i]]~dnorm(mu_prior_delta[i],1/pow(sd_prior_delta[i],2))
  delta_std[id_prior_delta[i]]<-delta[id_prior_delta[i]]/sd_tracer[id_prior_delta[i]]
}"

  signature_source_no_prior <-
    "for (i in 1:(nb_combinations-nb_prior_signature)){
  mean_signature_std[combinations[id_no_prior_signature[i],1],combinations[id_no_prior_signature[i],2]]~dnorm(0.01,0.01)
  mean_signature[combinations[id_no_prior_signature[i],1],combinations[id_no_prior_signature[i],2]]<-(mean_signature_std[combinations[id_no_prior_signature[i],1],combinations[id_no_prior_signature[i],2]]*sd_tracer[combinations[id_no_prior_signature[i],2]])+mean_tracer[combinations[id_no_prior_signature[i],2]]
  tau_signature[combinations[id_no_prior_signature[i],1],combinations[id_no_prior_signature[i],2]]~dgamma(0.01,0.01)
  var_signature[combinations[id_no_prior_signature[i],1],combinations[id_no_prior_signature[i],2]]<-1/tau_signature[combinations[id_no_prior_signature[i],1],combinations[id_no_prior_signature[i],2]]
}"


  signature_source_prior <-
    "for (i in 1:nb_prior_signature){
  tau2_signature[i]~dchisq(n_prior_signature[i]-1)
  tau_signature[combinations[id_prior_signature[i],1],combinations[id_prior_signature[i],2]]<-tau2_signature[i]/(n_prior_signature[i]*sd_prior_signature[i])
  var_signature[combinations[id_prior_signature[i],1],combinations[id_prior_signature[i],2]]<-1/tau_signature[combinations[id_prior_signature[i],1],combinations[id_prior_signature[i],2]]
  mean_signature_std[combinations[id_prior_signature[i],1],combinations[id_prior_signature[i],2]]~dnorm(mu_prior_signature[i],n_prior_signature[i]/pow(sd_prior_signature[i],2))
  mean_signature[combinations[id_prior_signature[i],1],combinations[id_prior_signature[i],2]]<-(mean_signature_std[combinations[id_prior_signature[i],1],combinations[id_prior_signature[i],2]]*sd_tracer[combinations[id_prior_signature[i],2]])+mean_tracer[combinations[id_prior_signature[i],2]]
}"

  diet <-
    "for (i in 1:nb_source_species){
  diet[id_source_species[i],1:nb_species]<-rep(0,nb_species)
}
for (i in 1:nb_consumer_single){
  diet[id_consumer_single[i],1]<-1
}
for (i in 1:nb_consumer_multiple){
  diet[id_consumer_multiple[i],1:nb_prey_per_species[id_consumer_multiple[i]]]~ddirich(alpha_diet[id_consumer_multiple[i],1:nb_prey_per_species[id_consumer_multiple[i]]])
}"

  compute_mean_signature <-
    "for (itra in 1:nb_tracer){
  sigma_random[itra]~dunif(0.001,10)
}

for(ispe in c(id_consumer_multiple,id_consumer_single)){
  for (itra in 1:nb_tracer){
    random_effect[ispe,itra]~dnorm(0,1/pow(sigma_random[itra],2)) #random effect
    var_signature[ispe,itra]<-1/tau[itra]+inprod(var_signature[prey_id[ispe,1:nb_prey_per_species[ispe]],itra],diet[ispe,1:nb_prey_per_species[ispe]]*diet[ispe,1:nb_prey_per_species[ispe]]) #variance of the signature for each species and tracer
    mean_signature_std[ispe,itra]<-inprod(diet[ispe,1:nb_prey_per_species[ispe]],mean_signature_std[prey_id[ispe,1:nb_prey_per_species[ispe]],itra])+random_effect[ispe,itra]+delta_std[itra]
    mean_signature[ispe,itra]<-(mean_signature_std[ispe,itra]*sd_tracer[itra])+mean_tracer[itra]
  }
}"

observation_process <-
  "for (indiv in 1:nb_signature){
  for (itra in 1:nb_tracer){
    pred[indiv,itra]<-mean_signature_std[id_species_signature[indiv],itra]
    isLeftCensored[indiv,itra]~dinterval(signature_data[indiv,itra],LOQ[indiv,itra])
    signature_data[indiv,itra]~dnorm(pred[indiv,itra],1/var_signature[id_species_signature[indiv],itra])
  }
}"

return(
  paste(
    "model{",
    priors_noise_regression,
    priors_delta_no_prior,
    ifelse(mydata$nb_prior_delta > 0, priors_delta_with_prior, ""),
    signature_source_no_prior,
    ifelse(mydata$nb_prior_signature > 0, signature_source_prior,""),
    diet,
    compute_mean_signature,
    observation_process,
    "}",
    sep = "\n"
  )
)
}
