#' compute ENA indices
#'
#' Computes several ENA indices for an escropath model. ENA indices are based
#' on formula of Guesnet et al. This function requires enaR package that is
#' not available on cran anymore but can be found on github
#' https://github.com/SEELab/enaR
#' @param mydata a list as returned by function \code{\link{prepare_data}}
#' @param myres a fit as returned by \code{\link{fit_escroc}}
#' @param quant quantiles of interest
#' @param biomassDet biomass of detritus
#' @param minTL minimum trophic level
#' @return a named array of indices
#' @references
#' @export

compute_ena <- function(mydata,
                        myres,
                        quant = NULL,
                        biomassDet = 0,
                        minTL = 1) {
  if (!requireNamespace("enaR", quietly = TRUE)) {
    stop("Package enaR needed for this function to work. Please install it.",
         call. = FALSE)
  }

  myfit_mat <- as.matrix(myres)
  if (nrow(myfit_mat)>1000)
    myfit_mat <- myfit_mat[seq(1, nrow(myfit_mat), length.out=1000),]
  nb_species <- mydata$nb_species
  species_names <- rownames(mydata$alpha_diet)[seq_len(nb_species)]
  id_flow <- expand.grid(species_names,
                         rownames(mydata$alpha_diet))
  ena_indices <- apply(myfit_mat, 1, function(iter){
    diet_name <- paste("diet[",id_flow$Var1,",",id_flow$Var2,"]",sep="")
    bio_name <- paste("biomass[",id_flow$Var1,"]",sep="")
    conso_name <- paste("consumption_rate[",id_flow$Var1,"]",sep="")


    consump <- iter[paste("consumption_rate[",species_names,"]",sep="")] *
      iter[paste("biomass[",species_names,"]",sep="")]
    input <- iter["input_Det"]
    export <- iter["export_Det"]
    production <- iter[paste("productivity[", species_names,"]", sep="")] *
      iter[paste("biomass[",species_names,"]",sep="")]
    production_PP <- iter["productivity[PP]"] *
      iter["biomass[PP]"]
    uq <- iter[paste("uq[", species_names,"]", sep="")]
    notassimilated <- consump*uq
    #mortality <- production * (1 -
    #iter[paste("trophic_efficiency[", species_names, "]", sep="")])
    bio <- iter[grep("biomass", names(iter))]

    mortality_PP <- production_PP * (1 - iter["trophic_efficiency[PP]"])
    respiration <- consump- #consumption
      production-#productiib
      notassimilated#nourriture non assimilee

    #matrix with the flow from prey (col) to predator (row)
    flow <- matrix(ifelse(diet_name %in% names(iter),
                          iter[diet_name],
                          0) *
                     iter[bio_name]*
                     iter[conso_name],
                   nb_species,
                   nb_species+2,
                   byrow=FALSE)
    flow <- rbind(flow,
                  matrix(0, 2, nb_species+2))
    consumed <- colSums(flow)

    mortality<-production-mydata$landings[seq_len(nb_species)]-
      mydata$discards[seq_len(nb_species)]-consumed[seq_len(nb_species)]

    rownames(flow) <- colnames(flow) <- rownames(mydata$alpha_diet)
    ##we now add flows towards detritus
    #productions which is not consumed and not assimilated food go to detritus
    flow[mydata$id_Det, seq_len(nb_species) ] <- mortality+
      + notassimilated + mydata$discards[seq_len(nb_species)]
    flow[mydata$id_Det, mydata$id_PP ] <- mortality_PP+
      mydata$discards[mydata$id_PP]

    ####flow is transposed to have from in row to to in cols as
    ####in enatool or enaR
    flow <- t(flow)
    exports=c(mydata$landings,export)
    biomassDet=0
    inputs = c(rep(0, nb_species),production_PP,
               input)
    mynetwork <- enaR::pack(flow,
                            respiration = c(respiration, 0, 0),
                            input = inputs,
                            storage=c(bio, biomassDet),
                            export = exports,
                            living=c(rep(TRUE, nb_species + 1), FALSE))
    etl <- enaR::enaTroAgg(mynetwork)$ETL
    k <- which(etl >= minTL)
    mtl <- sum((etl[k[-mydata$id_Det]] *
                  ((mydata$discards+mydata$landings)[k[-mydata$id_Det]]))/
                 sum((mydata$discards+mydata$landings)[k[-mydata$id_Det]]))

    TST <- sum(flow) + sum(respiration) +
      sum(input) + sum(exports)+ sum(mydata$landings) +
      sum(mydata$landings)

    APL <- TST/(sum(exports) + sum(respiration))

    ###not work
    Ti <- rowSums(flow) #outflow
    Tj <- colSums(flow) #inflow
    S <- Tj + c(rep(0, mydata$nb_species + 1), #total inflow
                input)
    #total outflow
    SS <- Ti + c(mydata$landings, export) + c(respiration, 0, 0)

    G <- sweep(flow, 1, S, "/")
    G[is.nan(G) | is.infinite(G)] <- 0
    Leontief <- solve(diag(nrow(G))-G)
    FCI <- sum(S/TST *(diag(Leontief)-1)/diag(Leontief))*100

    CCI <- 1.142 * FCI

    #AI
    scaled_Flow <- flow*TST /(matrix(SS,nrow(flow),ncol(flow))*
                                matrix(S,nrow(flow),ncol(flow),byrow = TRUE))
    Ai <- sum(ifelse(flow>0,
                     flow*log(scaled_Flow),
                     0)) * 1.442695

    #C capacity
    Ci <- - sum (ifelse(flow==0,
                        0,
                        flow*log(flow/TST))) * 1.442695

    Ai_Ci <- Ai/Ci

    res <- c(mtl, TST, APL,
             FCI, CCI, Ai, Ci, Ai_Ci)
    names(res) <- c("mtl",
                    "TST",
                    "APL",
                    "FCI",
                    "CCI",
                    "Ai",
                    "Ci",
                    "Ai_Ci")
    res

  })
  t(ena_indices)
}




