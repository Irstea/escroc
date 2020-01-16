#' Build the jags model code
#'
#' Based on data, this function build the jags model that will be used by runjages
#'
#' @param mydata a list data returned by \code{\link{prepare_data}}
#'
#' @return a string corresponding to the jags model
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
#' LOQ,prior_signature_data, prior_delta)
#'
#' #build the model
#' mymodel <- building_model(mydata)
#' @export
building_model <- function(mydata) {
  nb_species <- mydata$nb_species
  priors_noise_regression <-"
#precision of the regression of enrichment (uninformative prior)
for (i in 1:nb_tracer){
  tau[i]~dgamma(0.01,0.01)
}"

  priors_delta_no_prior <-"
#informative prior on enrichment/TMF
for (i in 1:(nb_tracer-nb_prior_delta)){
  delta_std[id_no_prior_delta[i]]~dnorm(0,.1)
  delta[id_no_prior_delta[i]]<-delta_std[id_no_prior_delta[i]]*sd_tracer[id_no_prior_delta[i]]
}"

  priors_delta_with_prior <-"
#mean signature of souces (eat Det or PP) drawn from informative prior
for (i in 1:(nb_prior_delta)){
  delta[id_prior_delta[i]]~dnorm(mu_prior_delta[i],1/pow(sd_prior_delta[i],2))
  delta_std[id_prior_delta[i]]<-delta[id_prior_delta[i]]/sd_tracer[id_prior_delta[i]]
}"

  signature_source_no_prior <-"
#mean signature of souces (eat Det or PP) drawn from uninformative prior
for (i in 1:(nb_combinations-nb_prior_signature)){
  mean_signature_std[combinations[id_no_prior_signature[i],1],combinations[id_no_prior_signature[i],2]]~dnorm(0.01,0.01)
  mean_signature[combinations[id_no_prior_signature[i],1],combinations[id_no_prior_signature[i],2]]<-(mean_signature_std[combinations[id_no_prior_signature[i],1],combinations[id_no_prior_signature[i],2]]*sd_tracer[combinations[id_no_prior_signature[i],2]])+mean_tracer[combinations[id_no_prior_signature[i],2]]
  tau_signature[combinations[id_no_prior_signature[i],1],combinations[id_no_prior_signature[i],2]]~dgamma(0.01,0.01)
  var_signature[combinations[id_no_prior_signature[i],1],combinations[id_no_prior_signature[i],2]]<-1/tau_signature[combinations[id_no_prior_signature[i],1],combinations[id_no_prior_signature[i],2]]
}

for (i in id_source_species){
  for (tra in 1:nb_tracer){
    random_effect[i,tra]<-0
  }
}
"


  signature_source_prior <-"
#prior on the signature of source species (eat Det or PP)
for (i in 1:nb_prior_signature){
  tau2_signature[i]~dchisq(n_prior_signature[i]-1)
  tau_signature[combinations[id_prior_signature[i],1],combinations[id_prior_signature[i],2]]<-tau2_signature[i]/(n_prior_signature[i]*sd_prior_signature[i])
  var_signature[combinations[id_prior_signature[i],1],combinations[id_prior_signature[i],2]]<-1/tau_signature[combinations[id_prior_signature[i],1],combinations[id_prior_signature[i],2]]
  mean_signature_std[combinations[id_prior_signature[i],1],combinations[id_prior_signature[i],2]]~dnorm(mu_prior_signature[i],n_prior_signature[i]/pow(sd_prior_signature[i],2))
  mean_signature[combinations[id_prior_signature[i],1],combinations[id_prior_signature[i],2]]<-(mean_signature_std[combinations[id_prior_signature[i],1],combinations[id_prior_signature[i],2]]*sd_tracer[combinations[id_prior_signature[i],2]])+mean_tracer[combinations[id_prior_signature[i],2]]
}"

  diet <-
    "
# prior for the diet based on dirichlet distribution
#dirichlet drawns should go in contiguous array, therefore we must first
# create a temporary object and then copy it appropriately into the diet matrix
for (i in id_consumer_single){
  diet_short[id_consumer_single,1]<-1
}

for (i in id_consumer_multiple){
  diet_short[i,1:nb_prey_per_species[i]]~ddirich(alpha_diet[i,1:nb_prey_per_species[i]])
}
for (i in 1:nb_species){
  for (p in 1:nb_prey_per_species[i]){
    diet[i,prey_id[i,p]]<-diet_short[i,p]
  }
  for (np in 1:(2+nb_species-nb_prey_per_species[i])){
    diet[i,not_prey_id[i,np]]<-0
  }
}
"


#in mean signature, we should probably add the possibility to have
#    migration ()
#    food from on external species
  compute_mean_signature <- "
#computation of mean signature given signature of the prey
for (itra in 1:nb_tracer){
  sigma_random[itra]~dunif(0.001,10)
}

for(ispe in id_not_source_species){
  for (itra in 1:nb_tracer){
    random_effect[ispe,itra]~dnorm(0,1/pow(sigma_random[itra],2)) #random effect
    var_signature[ispe,itra]<-1/tau[itra]+inprod(var_signature[prey_id[ispe,1:nb_prey_per_species[ispe]],itra],diet[ispe,prey_id[ispe,1:nb_prey_per_species[ispe]]]*diet[ispe,prey_id[ispe,1:nb_prey_per_species[ispe]]]) #variance of the signature for each species and tracer
    mean_signature_std[ispe,itra]<-inprod(diet[ispe,prey_id[ispe,1:nb_prey_per_species[ispe]]],mean_signature_std[prey_id[ispe,1:nb_prey_per_species[ispe]],itra])+random_effect[ispe,itra]+delta_std[itra]
    mean_signature[ispe,itra]<-(mean_signature_std[ispe,itra]*sd_tracer[itra])+mean_tracer[itra]
  }
}"

observation_process_signatures <-"
#observation process of signatures
for (indiv in 1:nb_signature){
  for (itra in 1:nb_tracer){
    pred[indiv,itra]<-mean_signature_std[id_species_signature[indiv],itra]
    isLeftCensored[indiv,itra]~dinterval(signature_data[indiv,itra],LOQ[indiv,itra])
    signature_data[indiv,itra]~dnorm(pred[indiv,itra],1/var_signature[id_species_signature[indiv],itra])
  }
}"


############################We know add ecopath in Escroc
####should we put migration somewhere?
flow_computation <- "
#computation of flows
for (spe in 1:nb_species){
	consumption[spe]<-biomass[spe]*consumption_rate[spe]
	production[spe]<-productivity[spe]*biomass[spe]
	not_assimilated[spe]<-consumption[spe]*uq[spe]
	respiration[spe]<-consumption[spe]-production[spe]-not_assimilated[spe]
}

for (spe in id_not_top_predator){
	other_mortality[spe]<-production[spe]*(1-trophic_efficiency[spe])
}
for (spe in id_top_predator){
	other_mortality[spe]<-production[spe]-landings[spe]-discards[spe]
  trophic_efficiency[spe]<-1-other_mortality[spe]/production[spe]
}


production[id_PP]<-productivity[id_PP]*biomass[id_PP]
other_mortality[id_PP]<-productivity[id_PP]*biomass[id_PP]*(1-trophic_efficiency[id_PP])
not_assimilated[id_PP]<-0
"

###should we put a EE to PP and detN?
priors_ecopath_model <-paste("
###priors on ecopath parameters
for (spe in 1:nb_species){
	#productivity[spe]~dunif(min_prod[spe],min(max_prod[spe],consumption_rate[spe]*
      #(.99999999-uq[spe])))
  A[spe]~dbeta(Aprior[spe,1],Aprior[spe,2])
	consumption_rate[spe]~dunif(min_cons_rate[spe],max_cons_rate[spe])
	uq[spe]~dunif(uq_min[spe],uq_max[spe])
  productivity[spe]<-A[spe]*consumption_rate[spe]*(1-uq[spe])
}

#priors for PP
trophic_efficiency[id_PP]~dbeta(prior_alpha[id_PP],prior_beta[id_PP])
productivity[id_PP] ~ dunif(min_prod[id_PP],max_prod[id_PP])

",
paste("trophic_efficiency[",mydata$id_not_top_predator,"]",
      ifelse(mydata$prior_alpha[mydata$id_not_top_predator] == 0,
             "<-0",
             paste("~dbeta(prior_alpha[",mydata$id_not_top_predator,
                   "],prior_beta[",mydata$id_not_top_predator,"])",sep="")),
             collapse="\n"),sep="")


priors_top_pred <- "
#priors for biomass of top predator
for (i in id_top_predator){
  biomass[i] ~ dunif((landings[i]+discards[i])/productivity[i],bmax)
}
"
#observations of biomass
observations_biomass <- "
#observations of biomass
for (i in 1:(nb_species+1)){
  obs_biomass[i]~dnorm(biomass[i],1/pow(sigma_biomass[i],2))
}
"

cons_top <- "
#top predator are not consumed
for (i in id_top_predator){
  consumed[i] <- 0
}"

cons<-"
###consumption for other species (not top predators)
for (spe in c(id_PP, id_not_top_predator)){
  biomass[spe]<-1/(trophic_efficiency[spe]*productivity[spe])*(discards[spe]+landings[spe]+inprod(diet[pred_id[spe,1:nb_pred_per_species[spe]],spe]*consumption_rate[pred_id[spe,1:nb_pred_per_species[spe]]],biomass[pred_id[spe,1:nb_pred_per_species[spe]]]))
  consumed[spe]<-inprod(diet[pred_id[spe,1:nb_pred_per_species[spe]],spe]*consumption_rate[pred_id[spe,1:nb_pred_per_species[spe]]],biomass[pred_id[spe,1:nb_pred_per_species[spe]]])
}
consumed[id_Det]<-inprod(diet[pred_id[id_Det,1:nb_pred_per_species[id_Det]],id_Det]*consumption_rate[pred_id[id_Det,1:nb_pred_per_species[id_Det]]],biomass[pred_id[id_Det,1:nb_pred_per_species[id_Det]]])
"


#For Detritus, we compute accumulation that should look like observation (input
#or output)
equilibrium <- paste("
#For Detritus, we compute accumulation that should look like observation (input
#or output)
  input_Det<-consumed[id_Det]-sum(not_assimilated+other_mortality+discards)
  input_Det_obs~dnorm(input_Det,tau_idp)
")



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
    observation_process_signatures,
    flow_computation,
    priors_ecopath_model,
    priors_top_pred,
    observations_biomass,
    cons_top,
    cons,
    equilibrium,
    "}",
    sep = "\n"
  )
)
}
