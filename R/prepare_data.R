#' preparation of the data
#'
#' prepares the data so that they can be used by jags to fit the model
#'
#' @param diet_matrix a square matrix with similar species as colnames and rownames. Species should be consistent with Species in signature_data and prior_signature_data
#' @param signatures_data a data frame. First column should be names species and followning columns should correspond to tracers. Data should be transformed first if required. Left censored data and missing data should be denoted NA
#' @param prior_signatures_data a data frame that contains data that can be used to build a prior signatures for source species. First column should be names Species and second one tracer. Species and tracers should be consistent with signatures_data. Then, the three following columns correpond to the mean signature, the standard deviation and the number of samples.
#' @param isLeftCensored a data frame with the same number of rows as signatures_data and one column per tracer. For each measure from the signatures_data, tells if the observation is left censored (0) or not (1). Put NA for missing data
#' @param LOQ same structure as isLeftCensored. Gives the limit of detection of limit of quantification of the corresponding measure
#'
#' @return a named list that contains formatted object that will be used by jags to fit the model
#' @examples
#' #importing data
#' data(prior_diet_matrix)
#' data(LOQ)
#' data(prior_signatures_data)
#'
#' #check that everything is ok
#' mydata=prepare_data(prior_diet_matrix,signatures_data,prior_signatures_data,isLeftCensored,LOQ)
#' @export
prepare_data <-
  function(diet_matrix,
           signatures_data,
           prior_signatures_data = NULL,
           isLeftCensored,
           LOQ) {
    #check that data are correctly formated
    checking_data(diet_matrix,
                  signatures_data,
                  prior_signatures_data,
                  isLeftCensored,
                  LOQ)

    #retrieve the tracer list
    nb_tracer <- ncol(signatures_data) - 1
    tracer_name <- sort(names(signatures_data)[-1])

    #sort the signatures_data columns
    signatures_data <- signatures_data[, c("Species", tracer_name)]
    LOQ <- LOQ[, tracer_name]
    isLeftCensored <- isLeftCensored[, tracer_name]

    #retrieve the species list
    species_name <- sort(rownames(diet_matrix))
    nb_species <- length(species_name)

    #sort the diet_matrix columns and rows
    diet_matrix <- diet_matrix[species_name,]
    diet_matrix <- diet_matrix[, species_name]

    #reorganize the diet matrix
    nb_prey_per_species <- rowSums(diet_matrix > 0)
    prey_id = matrix(NA, ncol(diet_matrix), ncol(diet_matrix))
    for (i in 1:nrow(diet_matrix)) {
      if (nb_prey_per_species[i] > 0) {
        prey_id[i, 1:nb_prey_per_species[i]] = which(diet_matrix[i,] > 0)
      }
    }
    colnames(prey_id) = rownames(prey_id) = rownames(diet_matrix)

    ###find source species
    id_source_species <- which(rowSums(diet_matrix) == 0)
    source_species <- names(id_source_species)
    nb_source_species <- length(source_species)

    #consumers that eat several species for which diet should be estimated
    id_consumer_multiple <-
      which(apply(diet_matrix, 1, function(x)
        sum(x > 0)) > 1)
    consumer_multiple <- names(id_consumer_multiple)
    nb_consumer_multiple <- length(consumer_multiple)

    #consumers that eat a single species for which diet should not be estimated
    id_consumer_single <-
      which(apply(diet_matrix, 1, function(x)
        sum(x > 0)) == 1)
    consumer_single <- names(id_consumer_single)
    nb_consumer_single <- length(consumer_single)

    #keep the identifier of the species from signature table
    id_species_signatures <-
      match(signatures_data$Species, species_name)

    #scaling of the signatures
    signatures_data <- scale(signatures_data[,-1])
    nb_signatures <- nrow(signatures_data)
    mean_tracer <- attr(signatures_data, "scaled:center")
    sd_tracer <- attr(signatures_data, "scaled:scale")

    #scaling the LOQ
    LOQ <- scale(LOQ, center = mean_tracer, scale = sd_tracer)

    #a table to store all the combination of species and tracer
    #this will be use to check wether we have prior data on signatures
    combinations <-
      expand.grid(Species = id_source_species, tracer = 1:nb_tracer)
    nb_combinations = nrow(combinations)
    nb_prior_signatures = nrow(prior_signatures_data)
    id_prior_signatures = match(
      paste(
        prior_signatures_data$Species,
        prior_signatures_data$tracer,
        sep = "__"
      ),
      paste(species_name[combinations$Species], tracer_name[combinations$tracer], sep =
              "__")
    )
    mu_prior_signatures = prior_signatures_data$mean
    sd_prior_signatures = prior_signatures_data$sd
    n_prior_signatures = prior_signatures_data$sd
    id_no_prior_signatures = 1:nrow(combinations)
    id_no_prior_signatures = id_no_prior_signatures[!id_no_prior_signatures %in%
                                                      id_prior_signatures]

    return(
      list(
        signatures_data = signatures_data,
        mean_tracer = mean_tracer,
        sd_tracer = sd_tracer,
        isLeftCensored = isLeftCensored,
        LOQ = LOQ,
        nb_consumer_multiple = nb_consumer_multiple,
        id_consumer_multiple = id_consumer_multiple,
        nb_consumer_single = nb_consumer_single,
        id_consumer_single = id_consumer_single,
        nb_source_species = nb_source_species,
        id_source_species = id_source_species,
        nb_species = nb_species,
        nb_tracer = nb_tracer,
        nb_signatures = nb_signatures,
        id_species_signatures = id_species_signatures,
        nb_prey_per_species = nb_prey_per_species,
        prey_id = prey_id,
        nb_prior_signatures = nb_prior_signatures,
        combinations = as.matrix(combinations),
        id_no_prior_signatures = id_no_prior_signatures,
        id_prior_signatures = id_prior_signatures,
        mu_prior_signatures = mu_prior_signatures,
        sd_prior_signatures = sd_prior_signatures,
        n_prior_signatures = n_prior_signatures,
        nb_combinations = nb_combinations

      )
    )

  }
