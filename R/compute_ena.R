  require(network)
  require(enaR)

compute_ena <- function(myfit, mydata, quant=NULL,biomassDet=0){
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
    ##we know add flows towards detritus
    ####productions which is not consumed and not assimilated food got to detritus
    flow[mydata$id_Det, seq_len(nb_species) ] <- mortality+
      + notassimilated + mydata$discards[seq_len(nb_species)]
    flow[mydata$id_Det, mydata$id_PP ] <- mortality_PP+
      mydata$discards[mydata$id_PP]


    exports=c(mydata$landings,ifelse(input>0,0,-input))

    mynetwork <- pack(t(flow),
                   respiration = c(respiration, 0, 0),
                   input = c(rep(0, nb_species),production_PP,
                             ifelse(input>0,input,0)),
                   storage=c(bio, biomassDet),
                   export = exports,
                   living=c(rep(TRUE, nb_species + 1), FALSE))
    mtl <- c(meanTrophicLevel(mynetwork,1))
    names(mtl) <- "MTL"
    enaf <- enaFlow(mynetwork)$ns[1, c("TST","APL","FCI")]
    enaf["CCI"] <- 1.142 * enaf["FCI"]
    enaa <- enaAscendency(mynetwork)[1, ]
    enaa["Ai_Ci"] <- enaa["A.internal"]/enaa["CAP.internal"]
    enaa <- enaa[c("A.internal","CAP.internal","Ai_Ci","ASC","CAP","ASC","OH")]
    tstf <- enaFlow(mynetwork)$ns[1,c("TST","TSTp")]
    others <- c(enaa["CAP"]/tstf["TSTp"],enaAscendency(mynetwork)[1, "AMI"])
    names(others) <- c("Flow_Div","AMI")
    return(c(mtl, enaf, enaa, tstf, others))})
  if (is.null(quant)){
    return(t(ena_indices))
  } else{
    apply(ena_indices,1,quantile,probs=quant,na.rm=TRUE)
  }
}

