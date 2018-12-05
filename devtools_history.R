
library(devtools)
library(roxygen2)
library(usethis)

usethis::use_build_ignore("devtools_history.R")
usethis::use_package("coda")
usethis::use_package("runjags")
usethis::use_package("ggplot2")
usethis::use_package("reshape2")
usethis::use_package("corrplot")

####exporting datasets
marine_wd<-"~/Documents/Bordeaux/equipe estuaire/TMF/Marine/SALMEC/"
signature_data=read.table(file = paste(marine_wd,"indiv_contam_isotop.csv",sep="/"),sep=";",na.strings = "NA",header=TRUE)
signature_data=na.omit(signature_data)
prior_diet_matrix=as.matrix(read.table(paste(marine_wd,"prior_trophic.csv",sep="/"),sep=";",header=TRUE,row.names = 1))
prior_diet_matrix[is.na(prior_diet_matrix)]=0
signature_data$Species=as.factor(sapply(strsplit(as.character(signature_data$species),"_"),function(x) paste(x[-length(x)],collapse ="_")))

signature_data=signature_data[,c(9,2:8)]
signature_data[signature_data==0]=NA
signature_data[,2:6]=log10(signature_data[,2:6])
signature_data=subset(signature_data,!signature_data$Species %in% c("Oyster","Scrobicularia"))



LOQ=c(log10(3*c(0.004,0.004,0.06,0.01,0.008)),c(-100,0))
min_signature=matrix(apply(signature_data[,-1],2,min,na.rm=TRUE),nrow=nrow(signature_data),ncol=ncol(signature_data)-1,byrow=TRUE)-0.1
LOQ=matrix(LOQ,nrow=nrow(signature_data),ncol=ncol(signature_data)-1,byrow=TRUE)
LOQ[!is.na(signature_data[,-1])]=min_signature[!is.na(signature_data[,-1])]
colnames(LOQ)=names(signature_data)[2:8]

data_val<-read.table(paste(marine_wd,"Gam_Cop_sources.csv",sep="/"),sep=";",header = TRUE)
data_val<-subset(data_val,data_val$Zone=="OM")
data_val<-subset(data_val,data_val$Groupe %in% c("Gammares","Copepodes"))
names(data_val)[6:7]=c('X15N',"X13C")
aggregation_prior=aggregate(data_val[,c('X15N',"X13C")],list(data_val$Groupe),function(x) c(mean(x),sd(x),length(x)))
prior_signature_data=data.frame("Species"=factor(levels=c("Copepods","Gammarids")),"tracer"=factor(levels=c("X15N","X13C")),"mean"=vector("numeric"),"sd"=vector("numeric"),"n"=vector("integer"))
prior_signature_data[1,3:5]=aggregation_prior[1,2]
prior_signature_data[1,1:2]=c("Copepods","X15N")
prior_signature_data[2,3:5]=aggregation_prior[1,3]
prior_signature_data[2,1:2]=c("Copepods","X13C")
prior_signature_data[3,3:5]=aggregation_prior[2,2]
prior_signature_data[3,1:2]=c("Gammarids","X15N")
prior_signature_data[4,3:5]=aggregation_prior[2,3]
prior_signature_data[4,1:2]=c("Gammarids","X13C")

usethis::use_data(prior_signature_data,prior_signature_data,overwrite=TRUE)
usethis::use_data(signature_data,signature_data,overwrite=TRUE)
usethis::use_data(LOQ,LOQ,overwrite=TRUE)
usethis::use_data(prior_diet_matrix,prior_diet_matrix,overwrite=TRUE)
usethis::use_vignette("Manual")

library(escrocR)
checking_data(prior_diet_matrix,signature_data,LOQ,prior_signature_data)
prior_delta <- data.frame(tracer=c("X15N","X13C"),mean=c(3,0),sd=c(1,1))

#check that everything is ok
mydata <- prepare_data(prior_diet_matrix,signature_data,
 LOQ,prior_signature_data,prior_delta)
#build the model
mymodel <- building_model(mydata)
#fit the model
myresults <- fit_escroc(mydata, mymodel,burnin=1e6,adapt=10000,sample=50000)
myresults <- window(myresults,thin=150)
usethis::use_data(myresults,myresults,overwrite=TRUE)
usethis::use_readme_md()
usethis::use_github_links()

document()
