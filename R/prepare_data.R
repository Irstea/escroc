#' preparation of the data
#'
#' prepares the data so that they can be used by jags to fit the model
#'
#' @param prior_diet_matrix a square matrix with similar species as colnames and
#'  rownames. Species should be consistent with Species in signature_data and
#'  prior_signature_data
#' @param signature_data a data frame. First column should be names species and
#' followning columns should correspond to tracers. Data should be transformed
#' first if required. Left censored data and missing data should be denoted NA
#' @param LOQ a data frame with the same number of rows as signature_data and
#' one column per tracer. Gives the limit of detection of limit of
#' quantification of the corresponding measure
#' @param prior_signature_data a data frame that contains data that can be used
#' to build a prior signature for source species. First column should be names
#' Species and second one tracer. Species and tracers should be consistent with
#' signature_data. Then, the three following columns correpond to the mean
#' signature, the standard deviation and the number of samples.
#' @param prior_delta a table that contains prior on magnification/enrichment of
#' tracers per trophic level. First column contains the name of the tracer,
#' second column contains the mean value and third column the standard deviation
#' @param biomass a table with to set priors for biomasses of species and
#'  primary production (PP). First column contains the name of species and PP.
#'  Second column stands for the value of the observed biomass and 3rd
#'  to the standard error arounnd the observation
#'  @param uq a table that stores for each species (first column), the minimum
#' (second column) and maximum values (3rd column)
#' of the prior on the proportion of unassimilated food
#' @param productivity a table that stores for each species (first column),
#' the minimum (second column) and maximum values (3rd column) of the prior on
#' the productivity (production=productivity x biomass)
#' @param consumption rate a table that stores for each species (first column),
#' the minimum (second column) and maximum values (3rd column) of the prior on
#' the consumption rate (consumption=consumption rate x biomass)
#' @param trophic_efficiency rate a table that stores for each species and PP
#' (1st column), a (second column) and b (3rd column) of the beta prior on
#' trophic efficiency (density x^(a-1)*(1-x)^(b-1)/Beta(a,b))
#' @param catch rate a table that stores for each species and PP
#' (1st column), a (second column) with landings and a third with discards#'
#'  @param input_detritus a table with only one row (for Detritus), first column
#'  necessarily "Detritus". Next column indicates mean and sd
#' @return a named list that contains formatted object that will be used by
#' jags to fit the model
#'
#' @importFrom fitdistrplus fitdist
#' @examples
#' #importing data
#' data(prior_diet_matrix)
#' data(LOQ)
#' data(prior_signature_data)
#' prior_delta <- data.frame(tracer=c("X15N","X13C"),mean=c(3,0),sd=c(1,1))
#'
#' #check that everything is ok
#' mydata <- prepare_data(prior_diet_matrix,signature_data,
#' LOQ,prior_signature_data,prior_delta)
#' @export
prepare_data <-
  function(prior_diet_matrix,
           signature_data,
           LOQ,
           prior_signature_data = NULL,
           prior_delta = NULL,
           biomass,
           uq,
           productivity,
           consumption_rate,
           trophic_efficiency,
           catch,
           input_detritus) {
    #check that data are correctly formated
    checking_data(
      prior_diet_matrix,
      signature_data,
      LOQ,
      prior_signature_data,
      prior_delta,
      biomass,
      uq,
      productivity,
      consumption_rate,
      trophic_efficiency,
      catch,
      input_detritus
    )

    #retrieve the tracer list
    nb_tracer <- ncol(signature_data) - 1
    tracer_name <- sort(names(signature_data)[-1])

    #sort the signature_data columns
    signature_data <- signature_data[, c("Species", tracer_name)]
    LOQ <- LOQ[, tracer_name]

    #retrieve the species list
    species_name <- sort(rownames(prior_diet_matrix))
    species_name <- species_name[-which(species_name %in% c("Detritus", "PP"))]
    nb_species <- length(species_name)

    #sort the prior_diet_matrix columns and rows
    prior_diet_matrix <- prior_diet_matrix[c(species_name,
                                             "PP",
                                             "Detritus"), ]
    prior_diet_matrix <- prior_diet_matrix[, c(species_name,
                                               "PP",
                                               "Detritus")]

    id_PP <- which(colnames(prior_diet_matrix) == "PP")
    id_Det <- which(colnames(prior_diet_matrix) == "Detritus")

    #reorganize the diet matrix
    #prey id is a matrix with predator in row and its preys in column
    nb_prey_per_species <- rowSums(prior_diet_matrix > 0)
    prey_id <- alpha_diet <- matrix(NA, ncol(prior_diet_matrix),
                                  ncol(prior_diet_matrix))
    for (i in 1:nrow(prior_diet_matrix)) {
      if (nb_prey_per_species[i] > 0) {
        prey_id[i, 1:nb_prey_per_species[i]] <-
          which(prior_diet_matrix[i, ] > 0)
        alpha_diet[i, 1:nb_prey_per_species[i]] <-
          prior_diet_matrix[i, which(prior_diet_matrix[i, ] > 0)]
      }
    }
    rownames(alpha_diet) <- rownames(prey_id) <- rownames(prior_diet_matrix)

    #we avoid that all alpha_diet per line equals 1 since it can raise to
    #numerical instability in jags
    for (i in seq_len(nb_species)) {
      if (sum(alpha_diet[i,]==1,na.rm=TRUE)==nb_prey_per_species[i])
        alpha_diet[i, 1:nb_prey_per_species[i]] <- .9
    }

    #pred_id is a matrix with prey in row and its predator in column
    nb_pred_per_species <- colSums(prior_diet_matrix > 0)
    pred_id <- matrix(NA, ncol(prior_diet_matrix),
                                  ncol(prior_diet_matrix))
    for (i in 1:ncol(prior_diet_matrix)) {
      if (nb_pred_per_species[i] > 0) {
        pred_id[i, 1:nb_pred_per_species[i]] <- which(prior_diet_matrix[, i] > 0)
      }
    }
    rownames(pred_id) <- rownames(prior_diet_matrix)
    not_prey_id <- as.matrix(prior_diet_matrix*0)
    for (spe in 1:ncol(prior_diet_matrix)){
      not_prey_id[spe,1:(nb_species+2-nb_prey_per_species[spe])]=which(prior_diet_matrix[spe,]==0)
    }


    #top predators: species that are not eaten
    id_top_predator <- which(colSums(prior_diet_matrix) == 0)
    top_predator <- names(id_top_predator)
    nb_top_predator <- length(id_top_predator)

    id_not_top_predator <- (1:nb_species)[-id_top_predator]

    #consumers that eat several species for which diet should be estimated
    id_consumer_multiple <-
      which(apply(prior_diet_matrix, 1, function(x)
        sum(x > 0)) > 1)
    consumer_multiple <- names(id_consumer_multiple)
    nb_consumer_multiple <- length(consumer_multiple)

    #consumers that eat a single species for which diet should not be estimated
    id_consumer_single <-
      which(apply(prior_diet_matrix, 1, function(x)
        sum(x > 0)) == 1)
    consumer_single <- names(id_consumer_single)
    nb_consumer_single <- length(consumer_single)


    #keep the identifier of the species from signature table
    id_species_signature <-
      match(signature_data$Species, species_name)

    #scaling of the signature
    signature_data <- scale(signature_data[, -1])
    nb_signature <- nrow(signature_data)
    mean_tracer <- attr(signature_data, "scaled:center")
    sd_tracer <- attr(signature_data, "scaled:scale")

    #scaling the LOQ
    LOQ <- scale(LOQ, center = mean_tracer, scale = sd_tracer)




    #a table to store all the combination of species and tracer
    #this will be use to check wether we have prior data on signature
    combinations <-
      expand.grid(Species = c(id_PP, id_Det), tracer = 1:nb_tracer)
    nb_combinations <- nrow(combinations)

    if (! is.null(prior_signature_data)) {
      nb_prior_signature <- nrow(prior_signature_data)
      if(!all(prior_signature_data$Species %in% c("PP","Detritus")))
        stop("prior signatures should be either on PP or Detritus")
      id_prior_signature <- match(
        paste(
          prior_signature_data$Species,
          prior_signature_data$tracer,
          sep = "__"
        ),
        paste(colnames(prior_diet_matrix)[combinations$Species],
            tracer_name[combinations$tracer], sep = "__")
      )
      mu_prior_signature <- (prior_signature_data$mean - mean_tracer[
                                    as.character(prior_signature_data$tracer)]) /
        sd_tracer[as.character(prior_signature_data$tracer)]
      sd_prior_signature <- prior_signature_data$sd /
        sd_tracer[as.character(prior_signature_data$tracer)]
      n_prior_signature <- prior_signature_data$n
    } else{
      nb_prior_signature <- 0
    }
    id_no_prior_signature <- 1:nrow(combinations)
    if (nb_prior_signature > 0)
      id_no_prior_signature <- id_no_prior_signature[!id_no_prior_signature %in%
                                                       id_prior_signature]

    #formatting prior delta data. A scaling is required to be constistent with
    # signature_data
    id_no_prior_delta <- 1:nb_tracer
    if (! is.null(prior_delta)) {
      id_prior_delta <- match(prior_delta$tracer, tracer_name)
      id_no_prior_delta <- id_no_prior_delta[!id_no_prior_delta %in%
                                             id_prior_delta]
      mu_prior_delta <- prior_delta$mean
      sd_prior_delta <- prior_delta$sd
      nb_prior_delta <- nrow(prior_delta)
    } else {
      nb_prior_delta <- 0
    }


    #creation of an indicator matrix telling if the data is left censored
    isLeftCensored <- (!(is.na(signature_data) & LOQ)) *1

    ####we now add the required information for ecopath part
    #first we put all tables in the good order
    biomass <- biomass[match(c(species_name, "PP"), biomass[, 1]), ]
    uq <- uq[match(species_name, uq[, 1]), ]
    productivity <- productivity[match(c(species_name, "PP"),
                                       productivity[, 1]), ]
    consumption_rate <- consumption_rate[match(species_name,
                                               consumption_rate[, 1]), ]
    trophic_efficiency <- trophic_efficiency[match(c(species_name, "PP"),
                                                   trophic_efficiency[, 1]), ]
    catch <- catch[match(c(species_name, "PP"),
                         catch[, 1]), ]


    Aprior<-t(sapply(seq_len(nb_species),function(s){
      uqp <- runif(1000, uq[s, 2], uq[s, 3])
      product <- runif(1000, productivity[s, 2], productivity[s, 3])
      cr <- runif(1000, consumption_rate[s, 2], consumption_rate[s, 3])
      A <- product/(cr*(1-uqp))
      A[A > 1] <- NA
      A[A <= 0] <- NA
      if (length(which(is.na(A))) > 0)
        A <- A[-which(is.na(A))]
      fitdist(A,"beta")$estimate
    }))

    bmax <- max(mapply(function(a, b) qnorm(.999,a,b),
                       biomass[, 2],
                       biomass[, 3]))
    mydata <- list(
      signature_data = as.matrix(signature_data),
      mean_tracer = mean_tracer,
      sd_tracer = sd_tracer,
      isLeftCensored = isLeftCensored,
      LOQ = as.matrix(LOQ),
      id_Det =  id_Det,
      id_PP = id_PP,
      id_consumer_multiple = id_consumer_multiple,
      id_not_top_predator = id_not_top_predator,
      id_top_predator=id_top_predator,
      nb_species = nb_species,
      nb_tracer = nb_tracer,
      nb_signature = nb_signature,
      id_species_signature = id_species_signature,
      nb_prey_per_species = nb_prey_per_species,
      prey_id = prey_id,
      not_prey_id=not_prey_id,
      pred_id=pred_id,
      nb_pred_per_species=nb_pred_per_species,
      nb_prior_signature = nb_prior_signature,
      combinations = as.matrix(combinations),
      id_no_prior_signature = id_no_prior_signature,
      nb_combinations = nb_combinations,
      nb_prior_delta = nb_prior_delta,
      id_no_prior_delta = id_no_prior_delta,
      alpha_diet = alpha_diet,
      bmax=bmax,
      obs_biomass=biomass[, 2],
      sigma_biomass=biomass[,3],
      prior_alpha=trophic_efficiency[, 2],
      prior_beta=trophic_efficiency[, 3],
      min_prod=productivity[, 2],
      max_prod=productivity[, 3],
      min_cons_rate=consumption_rate[, 2],
      max_cons_rate=consumption_rate[, 3],
      landings=catch[,2],
      discards=catch[,3],
      uq_min=uq[, 2],
      uq_max=uq[, 3],
      input_Det_obs=input_detritus[1, 2],
      tau_idp=1/input_detritus[1, 3]^2,
      Aprior=Aprior

    )
    if (nb_consumer_single > 0) {
      mydata <- c(mydata,
                  list(id_consumer_single = id_consumer_single))
    }

    if (nb_prior_signature > 0){
      mydata <- c(mydata,
                  list(id_prior_signature = id_prior_signature,
                       mu_prior_signature = mu_prior_signature,
                       sd_prior_signature = sd_prior_signature,
                       n_prior_signature = n_prior_signature))
    }

    if (nb_prior_delta > 0) {
      mydata <- c(mydata,
                  list(mu_prior_delta = mu_prior_delta,
                       sd_prior_delta = sd_prior_delta,
                       id_prior_delta = id_prior_delta))
    }


    return(mydata)

  }
