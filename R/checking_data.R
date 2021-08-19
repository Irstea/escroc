#' Checking data format
#'
#' check that the data are correctly formated. This function can be called by
#' the
#' user and will be called by \code{\link{prepare_data}}. Depending on
#' provided arguments, the function detects whether it is an escropath model
#' or an escroc model
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
#' second column contains the mean value and third column the standard
#' deviation
#' @param biomass a table with to set priors for biomasses of species and
#' primary production (PP). First column contains the name of species and PP.
#' Second column stands for the value of the observed biomass and 3rd
#' to the standard error around the observation (escropath only)
#' @param uq a table that stores for each species (first column), the minimum
#' (second column) and maximum values (3rd column)
#' of the prior on the proportion of unassimilated food (escropath only)
#' @param productivity a table that stores for each species (first column),
#' the minimum (second column) and maximum values (3rd column) of the prior on
#' the productivity (production=productivity x biomass) (escropath only)
#' @param consumption_rate a table that stores for each species (first column),
#' the minimum (second column) and maximum values (3rd column) of the prior on
#' the consumption rate (consumption=consumption rate x biomass) (escropath only)
#' @param trophic_efficiency rate a table that stores for each species and PP
#' (1st column), a (second column) and b (3rd column) of the beta prior on
#' trophic efficiency (density x^(a-1)*(1-x)^(b-1)/Beta(a,b)) (escropath only)
#' @param catch rate a table that stores for each species and PP
#' (1st column), a (second column) with landings and a third with discards
#' @param input_detritus a table with only one row (for Detritus), first column
#' necessarily "Detritus". Next column indicates mean and sd (escropath only)
#'
#' @return The function does not return any results but raises an error if a
#' missformat is detected
#' @examples
#' #importing data
#' data(signature_data)
#' data(prior_diet_matrix)
#' data(LOQ)
#' data(prior_signature_data)
#'
#' #check that everything is ok
#' checking_data(prior_diet_matrix,signature_data,
#' LOQ,prior_signature_data)
#' @export
checking_data <-
  function(prior_diet_matrix,
           signature_data,
           LOQ,
           prior_signature_data = NULL,
           prior_delta = NULL,
           biomass = NULL,
           uq = NULL,
           productivity = NULL,
           consumption_rate = NULL,
           trophic_efficiency = NULL,
           catch = NULL,
           input_detritus = NULL) {
    #check that prior_diet_matrix has similar species names in row and column
    if (!all(rownames(prior_diet_matrix) == colnames(prior_diet_matrix)))
      stop ("species names are different in row names and col names of the
            diet matrix")

    #this is an escropath model
    escropath <- ifelse(is.null(biomass), FALSE, TRUE)
    if (escropath){
      if (is.null(uq))
         stop("missing uq")
      if (is.null(productivity))
        stop("missing productivity")
      if (is.null(consumption_rate))
        stop("missing consumption rate")
      if (is.null(trophic_efficiency))
        stop("missing trophic efficiency")
      if (is.null(catch))
        stop("missing catch")
      if (is.null(input_detritus))
        stop("missing input_detritus")
      if (!("Detritus" %in% rownames(prior_diet_matrix)))
        stop("Detritus should be in the diet network")


      if (!("PP" %in% rownames(prior_diet_matrix)))
        stop("PP should be in the diet network")
    }

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

    #check that species have similar names in signature_data and
    #prior_diet_matrix
    list_species_diet <- row.names(prior_diet_matrix)
    if (escropath)
      list_species_diet <-
        list_species_diet[-which(list_species_diet %in%
                                   c("PP", "Detritus"))]
    list_species_signature <- unique(signature_data$Species)
    if (!all(sort(list_species_diet) == sort(list_species_signature)))
      stop ("species from the diet matrix do not match species from the
            signature data")



    #check that tracers names in the prior_signature_data exist
    #in the signature_data
    list_tracers_prior <- unique(prior_signature_data$tracer)
    if (sum(is.na(match(
      list_tracers_prior, names(signature_data)
    ))) > 0)
      stop ("a tracer name in list_tracers_prior is not consistent with the
            signature_data table")

    #check that species in the prior data are included in the species list
    if (!is.null(prior_signature_data)) {
      list_species_prior <- unique(prior_signature_data$Species)
      if (sum(is.na(match(
        list_species_prior, list_species_signature
      ))) > 0)
        stop ("some species in the prior table are not present in the signature
              table")

      #check that signature_data names of the table are correct
      if (!all(names(prior_signature_data) == c("Species",
                                                "tracer",
                                                "mean",
                                                "sd",
                                                "n")))
        stop("incorrect column names in prior_signature_data")
    }

    #check that LOQ are smaller than available signature_data
    if (sum(LOQ > signature_data[, -1], na.rm = TRUE) > 1)
      stop("some LOQ are greater then corresponding signature")

    #consistency between LOQ and signature_data
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
        stop ("some tracers in prior_delta are not present in the signature
              table")
    }
    if (escropath) {
      #check the format of biomass
      if (nrow(biomass) != (1 + length(list_species_diet)))
        stop("the biomass table should include 1 row per species and one row for
             PP")
      if (!all(sort(biomass[, 1]) == sort(c(list_species_diet, "PP"))))
        stop("in table biomass, first column should include all species + PP")
      if (class(biomass[, 2]) != "numeric" ||
          class(biomass[, 3]) != "numeric")
        stop("in table biomass, 2nd and 3rd column should be numeric")

      #check uq
      if (nrow(uq) != (length(list_species_diet)))
        stop("the uq table should include 1 row per species")
      if (!all(sort(uq[, 1]) == sort(list_species_diet)))
        stop("in table uq, first column should include all species")
      if (class(uq[, 2]) != "numeric" || class(uq[, 3]) != "numeric")
        stop("in table uq, 2nd and 3rd columns should be numeric")

      #check productivity
      if (nrow(productivity) != (length(list_species_diet) + 1))
        stop("the productivity table should include 1 row per species + PP")
      if (!all(sort(productivity[, 1]) == sort(c("PP", list_species_diet))))
        stop("in table productivity, first column should include all species+PP")
      if (class(productivity[, 2]) != "numeric" ||
          class(productivity[, 3]) != "numeric")
        stop("in table productivity, 2nd and 3rd columns should be numeric")


      #check catch
      if (nrow(catch) != (length(list_species_diet) + 1))
        stop("the productivity table should include 1 row per species + PP")
      if (!all(sort(catch[, 1]) == sort(c("PP", list_species_diet))))
        stop("in table productivity, first column should include all species+PP")
      if (class(catch[, 2]) != "numeric"  ||
          class(catch[, 3]) != "numeric")
        stop("in table catch, 2nd and 3rd columns should be numeric")

      #check consumption_rate
      if (nrow(consumption_rate) != (length(list_species_diet)))
        stop("the consumption_rate table should include 1 row per species")
      if (!all(sort(consumption_rate[, 1]) == sort(list_species_diet)))
        stop("in table consumption_rate, first column should include all species")
      if (class(consumption_rate[, 2]) != "numeric" ||
          class(consumption_rate[, 3]) != "numeric")
        stop("in table consumption_rate, 2nd and 3rd columns should be numeric")

      #check the format of trophic_efficiency
      if (nrow(trophic_efficiency) != (1 + length(list_species_diet)))
        stop("the trophic_efficiency table should include 1 row per species
             and one row for PP")
      if (!all(sort(trophic_efficiency[, 1]) == sort(c(list_species_diet, "PP"))))
        stop("in table trophic_efficiency, first column should include all species
             + PP")
      if (class(trophic_efficiency[, 2]) != "numeric" ||
          class(trophic_efficiency[, 3]) != "numeric")
        stop("in table trophic_efficiency, 2nd and 3rd columns should be numeric")

      #check the format of detritus
      if (input_detritus[, 1] != "Detritus")
        stop("in table input_detritus, first column should be Detritus")
      if (class(input_detritus[, 2]) != "numeric" ||
          class(input_detritus[, 3]) != "numeric")
        stop("in table input_detritus, 2nd and 3rd columns should be numeric")
    }
  }
