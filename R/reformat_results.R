#' Result formating
#'
#' This function is called by \code{\link{fit_escroc}} in order to format the
#' results. This function should not be called directly by the user
#'
#' @param myfit a runjags object as returned by \code{\link[runjags]{run.jags}}
#' @param mydata a data object as returned by \code{\link{prepare_data}}
#'
#' @return a formatted \code{\link[coda]{mcmc.list}}
#' @importFrom coda as.mcmc.list
#' @importFrom coda varnames
#' @export
reformat_results <- function(myfit, mydata) {
  escropath <- ! is.null(mydata$obs_biomass)
  diet_name <- ifelse(escropath,
                      "diet_short",
                      "diet")
  res <- as.mcmc.list(myfit)
  kept_diet <- which(!is.na(mydata$prey_id), arr.ind = TRUE)
  names_keep <-
    mapply(function(i, j)
      paste0(diet_name,"[", i, ",", j, "]"),
      kept_diet[, 1],
      kept_diet[, 2])
  names_variables <- varnames(res)
  keep <-
    c(
      grep("random_effect", names_variables),
      grep("delta", names_variables),
      grep("mean_signature", names_variables),
      match(names_keep, names_variables),
      grep("consumption_rate", names_variables),
      grep("trophic_efficiency", names_variables),
      grep("uq", names_variables),
      grep("consumption_rate", names_variables),
      grep("productivity", names_variables),
      grep("input_Det", names_variables),
      grep("export_Det", names_variables),
      grep("A", names_variables),
      grep("biomass", names_variables))

  res <- res[, keep, drop = FALSE]

  names_variables <- varnames(res)

  species_name <- row.names(mydata$prey_id)
  tracer_name <- colnames(mydata$signature_data)

  ####renaming of the varnames
  new_names <- sapply(names_variables, function(old_name) {
    if (length(grep("\\[",old_name))>0){
      pieces <-
        strsplit(substr(old_name, 1, nchar(old_name) - 1), c("[\\[,]"))[[1]]
      if (startsWith(old_name, "random") ||
          startsWith(old_name,"mean_signature")) {
        (paste0(pieces[1],
               "[",
               species_name[as.integer(pieces[2])],
               ",",
               tracer_name[as.integer(pieces[3])], "]"
               ))
      } else if (startsWith(old_name, "delta")) {
        (paste(pieces[1], "[", tracer_name[as.integer(pieces[2])], "]", sep =
                       ""))
      } else if (startsWith(old_name,"diet")){
        idspec = as.integer(pieces[2])
        (paste(pieces[1],
               "[",
               species_name[idspec],
               ",",
               species_name[mydata$prey_id[idspec, as.integer(pieces[3])]], "]",
               sep = ""))
      } else {
        paste(pieces[1], "[", species_name[as.integer(pieces[2])], "]", sep =
                       "")

      }
    } else{
      old_name
    }})
  coda::varnames(res) <- new_names
  if (escropath)
    coda::varnames(res) <- gsub("diet_short", "diet", coda::varnames(res))
  return(res)
}

