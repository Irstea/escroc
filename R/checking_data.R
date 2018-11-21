#' Checking data format
#' check that the data are correctly formated. This function can be called by the user and will be called by \code{\link{prepare_data}}
#'
#' @param diet_matrix a square matrix with similar species as colnames and rownames. Species should be consistent with Species in signature_data and prior_signature_data
#' @param signatures_data a data frame. First column should be names species and followning columns should correspond to tracers. Data should be transformed first if required. Left censored data and missing data should be denoted NA
#' @param prior_signatures_data a data frame that contains data that can be used to build a prior signatures for source species. First column should be names Species and second one tracer. Species and tracers should be consistent with signatures_data. Then, the three following columns correpond to the mean signature, the standard deviation and the number of samples.
#' @param isLeftCensored a data frame with the same number of rows as signatures_data and one column per tracer. For each measure from the signatures_data, tells if the observation is left censored (0) or not (1). Put NA for missing data
#' @param LOQ same structure as isLeftCensored. Gives the limit of detection of limit of quantification of the corresponding measure
#'
#' @return The function does not return any results but raises an error if a missformat is detected
#' @examples
#' #importing data
#' data(signatures_data)
#' data(prior_diet_matrix)
#' data(LOQ)
#' data(prior_signatures_data)
#'
#' #check that everything is ok
#' checking_data(prior_diet_matrix,signatures_data,prior_signatures_data,isLeftCensored,LOQ)
#' @export
checking_data <-
  function(diet_matrix,
           signatures_data,
           prior_signatures_data = NULL,
           isLeftCensored,
           LOQ) {
    #check that diet_matrix has similar species names in row and column
    if (!all(rownames(diet_matrix) == colnames(diet_matrix)))
      stop ("species names are different in row names and col names of the diet matrix")


    #check that signatures_data first column is indeed named Species
    if (names(signatures_data)[1] != "Species")
      stop("first column of signatures_data should be names Species")

    #check that all value are positive
    if (sum(diet_matrix < 0) > 0)
      stop ("elements of diet matrix should all be positive")

    if (nrow(diet_matrix) < 2)
      stop("not enough species in the diet matrix")

    #check the absence of loop in the diet matrix
    if (sum(diag(diet_matrix) != 0))
      stop("diagonal of the diet_matrix should be set to zero")
    power_mat <- diet_matrix
    for (i in 2:nrow(diet_matrix)) {
      power_mat <- power_mat %*% diet_matrix
      if (sum(diag(diet_matrix)) != 0)
        stop(paste("there's a cycle of length", i, "in the diet-matrix"))
    }

    #check that species have similar names in signatures_data and diet_matrix
    list_species_diet <- row.names(diet_matrix)
    list_species_signature <- unique(signatures_data$Species)
    if (!all(sort(list_species_diet) == sort(list_species_signature)))
      stop ("species from the diet matrix do not match species from the signature data")



    #check that tracers names in the prior_signatures_data exist in the signatures_data
    list_tracers_prior <- unique(prior_signatures_data$tracer)
    if (sum(is.na(match(
      list_tracers_prior, names(signatures_data)
    ))) > 0)
      stop ("a tracer name in list_tracers_prior is not consistent with the signatures_data table")

    #check that species in the prior data are included in the species list
    if (!is.null(prior_signatures_data)) {
      list_species_prior <- unique(prior_signatures_data$Species)
      if (sum(is.na(match(
        list_species_prior, list_species_signature
      ))) > 0)
        stop ("some species in the prior table are not present in the signature table")

      #check that signatures_data names of the table are correct
      if (!all(names(prior_signatures_data) == c("Species", "tracer", "mean", "sd", "n")))
        stop("incorrect column names in prior_signatures_data")
    }

    #check that tracers names match in LOQ, is Above, signatures and prior
    if (!all(dim(isLeftCensored) == dim(LOQ)))
      stop ("LOQ and isLeftCensored should have similar dimensions")
    if (!all(dim(isLeftCensored) == dim(signatures_data[, -1])))
      stop ("inconsistent dimentions between isLeftcensored and signatures")
    if (!all(sort(names(isLeftCensored)) == sort(names(names(signatures_data)))))
      stop ("names are not consistent between signatures_data and isLeftCensored")
    if (!all(sort(names(LOQ)) == sort(names(names(signatures_data)))))
      stop ("names are not consistent between signatures_data and LOQ")
  }
