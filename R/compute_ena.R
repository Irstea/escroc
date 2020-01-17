  require(network)
  require(enaR)

compute_ena <- function(myfit, mydata, quant=c(.025,.5,.975)){
  myfit_mat <- as.matrix(myfit)
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
    production <- iter[paste("productivity[", species_names,"]", sep="")] *
      iter[paste("biomass[",species_names,"]",sep="")]
    production_PP <- iter["productivity[PP]"] *
      iter["biomass[PP]"]
    uq <- iter[paste("uq[", species_names,"]", sep="")]
    notassimilated <- consump*uq
    #mortality <- production * (1 -
                #iter[paste("trophic_efficiency[", species_names, "]", sep="")])

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
    ##we know add flows towards detritus
    ####productions which is not consumed and not assimilated food got to detritus
    flow[mydata$id_Det, seq_len(nb_species) ] <- mortality+
      + notassimilated + mydata$discards[seq_len(nb_species)]
    flow[mydata$id_Det, mydata$id_PP ] <- mortality_PP+
      mydata$discards[mydata$id_PP]


    exports=c(mydata$landings,0)

    mynetwork <- pack(t(flow),
                   respiration = c(respiration, 0, 0),
                   input = c(rep(0, nb_species),production_PP,
                             input),
                   export = exports,
                   living=c(rep(TRUE, nb_species + 1), FALSE),
                   storage = rep(0, nb_species + 2))
    tmp <- enaFlow(mynetwork)$ns
    indices=as.vector(tmp)
    names(indices)=colnames(tmp)
    return(indices)})
  apply(ena_indices,1,quantile,probs=quant,na.rm=TRUE)
}

