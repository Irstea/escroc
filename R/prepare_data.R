#' preparation of the data
#'
#' prepares the data so that they can be used by jags to fit the model
#'
#' @param prior_diet_matrix a square matrix with similar species as colnames and rownames. Species should be consistent with Species in signature_data and prior_signature_data
#' @param signature_data a data frame. First column should be names species and followning columns should correspond to tracers. Data should be transformed first if required. Left censored data and missing data should be denoted NA
#' @param LOQ a data frame with the same number of rows as signature_data and one column per tracer. Gives the limit of detection of limit of quantification of the corresponding measure
#' @param prior_signature_data a data frame that contains data that can be used to build a prior signature for source species. First column should be names Species and second one tracer. Species and tracers should be consistent with signature_data. Then, the three following columns correpond to the mean signature, the standard deviation and the number of samples.
#' @param prior_delta a table that contains prior on magnification/enrichment of tracers per trophic level. First column contains the name of the tracer, second column contains the mean value and third column the standard deviation
#'
#' @return a named list that contains formatted object that will be used by jags to fit the model
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
           prior_delta = NULL) {
    #check that data are correctly formated
    checking_data(
      prior_diet_matrix,
      signature_data,
      LOQ,
      prior_signature_data,
      prior_delta
    )

    #retrieve the tracer list
    nb_tracer <- ncol(signature_data) - 1
    tracer_name <- sort(names(signature_data)[-1])

    #sort the signature_data columns
    signature_data <- signature_data[, c("Species", tracer_name)]
    LOQ <- LOQ[, tracer_name]

    #retrieve the species list
    species_name <- sort(rownames(prior_diet_matrix))
    nb_species <- length(species_name)

    #sort the prior_diet_matrix columns and rows
    prior_diet_matrix <- prior_diet_matrix[species_name, ]
    prior_diet_matrix <- prior_diet_matrix[, species_name]

    #reorganize the diet matrix
    nb_prey_per_species <- rowSums(prior_diet_matrix > 0)
    prey_id = alpha_diet = matrix(NA, ncol(prior_diet_matrix), ncol(prior_diet_matrix))
    for (i in 1:nrow(prior_diet_matrix)) {
      if (nb_prey_per_species[i] > 0) {
        prey_id[i, 1:nb_prey_per_species[i]] <- which(prior_diet_matrix[i, ] > 0)
        alpha_diet[i, 1:nb_prey_per_species[i]] <- prior_diet_matrix[i, which(prior_diet_matrix[i, ] > 0)]
      }
    }
    rownames(alpha_diet) <- rownames(prey_id) <- rownames(prior_diet_matrix)

    ###find source species
    id_source_species <- which(rowSums(prior_diet_matrix) == 0)
    source_species <- names(id_source_species)
    nb_source_species <- length(source_species)

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

    id_consumer <- c(id_consumer_multiple,
                     id_consumer_single)

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
      expand.grid(Species = id_source_species, tracer = 1:nb_tracer)
    nb_combinations <- nrow(combinations)
    if (! is.null(prior_signature_data)) {
      nb_prior_signature <- nrow(prior_signature_data)
      id_prior_signature <- match(
        paste(
          prior_signature_data$Species,
          prior_signature_data$tracer,
          sep = "__"
        ),
        paste(species_name[combinations$Species], tracer_name[combinations$tracer], sep =
                "__")
      )
      mu_prior_signature <- (prior_signature_data$mean - mean_tracer[as.character(prior_signature_data$tracer)]) /
        sd_tracer[as.character(prior_signature_data$tracer)]
      sd_prior_signature <- prior_signature_data$sd / sd_tracer[as.character(prior_signature_data$tracer)]
      n_prior_signature <- prior_signature_data$n
    } else {
      nb_prior_signature <- 0
    }

    id_no_prior_signature <- 1:nrow(combinations)
    if (nb_prior_signature > 0)
      id_no_prior_signature <- id_no_prior_signature[!id_no_prior_signature %in%
                                                      id_prior_signature]

    #formatting prior delta data. A scaling is required to be constistent with signature_data
    id_no_prior_delta <- 1:nb_tracer
    if (!is.null(prior_delta)){
      id_prior_delta <- match(prior_delta$tracer, tracer_name)
      id_no_prior_delta <- id_no_prior_delta[!id_no_prior_delta %in% id_prior_delta]
      mu_prior_delta <- prior_delta$mean
      sd_prior_delta <- prior_delta$sd
      nb_prior_delta <- nrow(prior_delta)
    } else {
      nb_prior_delta <- 0
    }


    #creation of an indicator matrix telling if the data is left censored
    isLeftCensored <- (!(is.na(signature_data) & LOQ)) *1


    mydata <- list(
      signature_data = signature_data,
      mean_tracer = mean_tracer,
      sd_tracer = sd_tracer,
      isLeftCensored = isLeftCensored,
      LOQ = LOQ,
      nb_consumer_multiple = nb_consumer_multiple,
      id_consumer_multiple = id_consumer_multiple,
      nb_consumer_single = nb_consumer_single,
      nb_source_species = nb_source_species,
      id_source_species = id_source_species,
      nb_species = nb_species,
      nb_tracer = nb_tracer,
      nb_signature = nb_signature,
      id_species_signature = id_species_signature,
      nb_prey_per_species = nb_prey_per_species,
      prey_id = prey_id,
      nb_prior_signature = nb_prior_signature,
      combinations = as.matrix(combinations),
      id_no_prior_signature = id_no_prior_signature,
      nb_combinations = nb_combinations,
      nb_prior_delta = nb_prior_delta,
      id_no_prior_delta = id_no_prior_delta,
      alpha_diet = alpha_diet,
      id_consumer = id_consumer

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
                       id_prior_delta = id_prior_delta
                  ))
    }

    return(
      mydata
    )

  }
