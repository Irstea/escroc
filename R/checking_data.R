#' Checking data format
#'
#' check that the data are correctly formated. This function can be called by the
#' user and will be called by \code{\link{prepare_data}}
#'
#' @param prior_diet_matrix a square matrix with similar species as colnames and rownames. Species should be consistent with Species in signature_data and prior_signature_data
#' @param signature_data a data frame. First column should be names species and followning columns should correspond to tracers. Data should be transformed first if required. Left censored data and missing data should be denoted NA
#' @param isLeftCensored a data frame with the same number of rows as signature_data and one column per tracer. For each measure from the signature_data, tells if the observation is left censored (0) or not (1). Put NA for missing data
#' @param LOQ same structure as isLeftCensored. Gives the limit of detection of limit of quantification of the corresponding measure
#' @param prior_signature_data a data frame that contains data that can be used to build a prior signature for source species. First column should be names Species and second one tracer. Species and tracers should be consistent with signature_data. Then, the three following columns correpond to the mean signature, the standard deviation and the number of samples.
#' @param prior_delta a table that contains prior on magnification/enrichment of tracers per trophic level. First column contains the name of the tracer, second column contains the mean value and third column the standard deviation
#'
#' @return The function does not return any results but raises an error if a missformat is detected
#' @examples
#' #importing data
#' data(signature_data)
#' data(prior_diet_matrix)
#' data(LOQ)
#' data(prior_signature_data)
#'
#' #check that everything is ok
#' checking_data(prior_diet_matrix,signature_data,
#' prior_signature_data,isLeftCensored,LOQ)
#' @export
checking_data <-
  function(prior_diet_matrix,
           signature_data,
           isLeftCensored,
           LOQ,
           prior_signature_data = NULL,
           prior_delta = NULL) {
    #check that prior_diet_matrix has similar species names in row and column
    if (!all(rownames(prior_diet_matrix) == colnames(prior_diet_matrix)))
      stop ("species names are different in row names and col names of the diet matrix")


    #check that signature_data first column is indeed named Species
    if (names(signature_data)[1] != "Species")
      stop("first column of signature_data should be names Species")

    #check that all value are positive
    if (sum(prior_diet_matrix < 0) > 0)
      stop ("elements of diet matrix should all be positive")

    if (nrow(prior_diet_matrix) < 2)
      stop("not enough species in the diet matrix")

    #check the absence of loop in the diet matrix
    if (sum(diag(prior_diet_matrix) != 0))
      stop("diagonal of the prior_diet_matrix should be set to zero")
    power_mat <- prior_diet_matrix
    for (i in 2:nrow(prior_diet_matrix)) {
      power_mat <- power_mat %*% prior_diet_matrix
      if (sum(diag(prior_diet_matrix)) != 0)
        stop(paste("there's a cycle of length", i, "in the diet-matrix"))
    }

    #check that species have similar names in signature_data and prior_diet_matrix
    list_species_diet <- row.names(prior_diet_matrix)
    list_species_signature <- unique(signature_data$Species)
    if (!all(sort(list_species_diet) == sort(list_species_signature)))
      stop ("species from the diet matrix do not match species from the signature data")



    #check that tracers names in the prior_signature_data exist in the signature_data
    list_tracers_prior <- unique(prior_signature_data$tracer)
    if (sum(is.na(match(
      list_tracers_prior, names(signature_data)
    ))) > 0)
      stop ("a tracer name in list_tracers_prior is not consistent with the signature_data table")

    #check that species in the prior data are included in the species list
    if (!is.null(prior_signature_data)) {
      list_species_prior <- unique(prior_signature_data$Species)
      if (sum(is.na(match(
        list_species_prior, list_species_signature
      ))) > 0)
        stop ("some species in the prior table are not present in the signature table")

      #check that signature_data names of the table are correct
      if (!all(names(prior_signature_data) == c("Species", "tracer", "mean", "sd", "n")))
        stop("incorrect column names in prior_signature_data")
    }

    #check that tracers names match in LOQ, is Above, signature and prior
    if (!all(dim(isLeftCensored) == dim(LOQ)))
      stop ("LOQ and isLeftCensored should have similar dimensions")
    if (!all(dim(isLeftCensored) == dim(signature_data[,-1])))
      stop ("inconsistent dimentions between isLeftcensored and signature")
    if (!all(sort(names(isLeftCensored)) == sort(names(names(signature_data)))))
      stop ("names are not consistent between signature_data and isLeftCensored")
    if (!all(sort(names(LOQ)) == sort(names(names(signature_data)))))
      stop ("names are not consistent between signature_data and LOQ")


    #check the format the prior_delta table
    if (!is.null(prior_delta)) {
      if (!(all(names(prior_delta) == c("tracer", "mean", "sd"))))
        stop("incorrect names in prior_delta")
      #check that tracers are correct
      list_tracer_prior <- unique(prior_delta$tracer)
      if (sum(is.na(match(
        list_tracer_prior, names(signature_data)
      ))) > 0)
        stop ("some tracers in prior_delta are not present in the signature table")
    }
  }
