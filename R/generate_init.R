#' Generate init
#'
#' This function is called by \code{\link{fit_escroc}} to
#' generate relevant initials values for the 3 Markov chains.
#'
#' @param mydata a list generated by \code{\link{prepare_data}}
#'
#' @return a list of lists of initial values for MCMC chains
#' @importFrom stats runif
#' @importFrom stats rnorm
#' @export
generate_init <- function(mydata) {
  escropath <- ! is.null(mydata$obs_biomass)
  nb_species <- mydata$nb_species
  lapply(1:3, function(chain) {
    delta <- delta_std <- rep(NA,mydata$nb_tracer)
    delta[mydata$id_prior_delta] <- mapply(function(mu,sd) rnorm(1,mu,sd),
                                           mydata$mu_prior_delta,
                                           mydata$sd_prior_delta)
    delta_std [!((1:mydata$nb_tracer) %in% mydata$id_prior_delta)] <-
      runif(mydata$nb_tracer-mydata$nb_prior_delta,.1,1)
    if (escropath) {
      diet_short <- matrix(NA, nb_species, max(mydata$nb_prey_per_species))
      for (i in 1:(nb_species+2)) {
        if (mydata$nb_prey_per_species[i] > 1) {
          alpha <- runif(mydata$nb_prey_per_species[i])
          diet_short[i, 1:mydata$nb_prey_per_species[i]] <-
            alpha / sum(alpha)
        }
      }
    } else {
      diet <- matrix(NA, mydata$nb_species, mydata$nb_species)
      for (i in 1:mydata$nb_species) {
        if (mydata$nb_prey_per_species[i] > 1) {
          alpha <- runif(mydata$nb_prey_per_species[i])
          diet[i, 1:mydata$nb_prey_per_species[i]] <-
            alpha / sum(alpha)
        }
      }
    }
    mean_signature_std <-
      matrix(NA, mydata$nb_species + ifelse(escropath, 2, 0), mydata$nb_tracer)
    if (escropath) {
      mean_signature_std[c(mydata$id_Det, mydata$id_PP),] <-
        runif(2 * mydata$nb_tracer,-3, 1)
    } else {
      mean_signature_std[mydata$id_source_species,] <-
      runif(mydata$nb_source_species * mydata$nb_tracer,-3, 1)
    }
    if (escropath)  {
      random_effect <-
        matrix(
          runif(nb_species * mydata$nb_tracer),
          nb_species,
          mydata$nb_tracer
        )
    } else {
      random_effect <-
        matrix(
          runif(mydata$nb_prey_per_species * mydata$nb_tracer),
          mydata$nb_species,
          mydata$nb_tracer
        )
      random_effect[mydata$id_source_species,] <- NA
    }
    signature_data <-
      matrix(NA,
             nrow(mydata$signature_data),
             ncol(mydata$signature_data))
    random_signature <-
      mydata$LOQ + matrix(runif(prod(dim(signature_data)),-10,-.1),
                   nrow(signature_data),
                   ncol(signature_data))
    signature_data[is.na(mydata$signature_data)] <-
      random_signature[is.na(mydata$signature_data)]
  if (escropath) {
    trophic_efficiency <- mapply(function(a,b){
      if (a==0){
        return(NA)
      } else {
        return(max(0.01,min(0.99,rbeta(1, a, b))))
      }
    }, mydata$prior_alpha, mydata$prior_beta)
    trophic_efficiency[mydata$id_top_predator] <- NA
    consumption_rate <- mapply(function(a,b) runif(1,a,b),
                               mydata$min_cons_rate, mydata$max_cons_rate)

    productivity <- c(rep(NA, nb_species),
                      runif(1,
                            mydata$min_prod[mydata$id_PP],
                            mydata$max_prod[mydata$id_PP]))
    A <- mapply(function(a,b) rbeta(1,a,b),
                mydata$Aprior[,1], mydata$Aprior[,1])

    uq <- mapply(function(a,b) runif(1,a,b),
                 mydata$uq_min, mydata$uq_max)
    tmp <- c(A*consumption_rate*(1-uq), productivity[mydata$id_PP])


    biomass <- mapply(function(is_top,min) ifelse(is_top,
                                                  runif(1,min,mydata$bmax),
                                                  NA),
                      seq_len(nb_species+1) %in% mydata$id_top_predator,
                      (mydata$landings+mydata$discards)/tmp)


  }
    res <- list(
             mean_signature_std = mean_signature_std,
             random_effect = random_effect,
             signature_data = signature_data,
             delta = delta,
             delta_std =delta_std,
             .RNG.seed = chain,
             .RNG.name = c(
               "base::Super-Duper",
               "base::Wichmann-Hill",
               "base::Mersenne-Twister"
             )[chain]
    )
    if (escropath) {
      res <- c(res,
               list(
               diet_short = diet_short,
               trophic_efficiency=trophic_efficiency,
               consumption_rate=consumption_rate,
               A=A,
               uq=uq,
               productivity=productivity,
               biomass=biomass)[chain]
             )
    } else {
       res <- c(res,
                list(diet = diet))
    }
    res
  })
}
