#' Contaminants and isotopic measures
#'
#' Measures of concentrations of different perfluoroalkylated
#' substances and isotopic ratios on 138 individuals from 16 species
#' in the Gironde estuary collected (Munoz et al. 2017). NA correspond
#' to left censored observations
#'
#' @format A data frame with 53940 rows and 10 variables:
#' \describe{
#'   \item{Species}{name of the species}
#'   \item{FOSA}{log10 concentration in FOSA}
#'   \item{L.PFOS}{log10 concentration in L.PFOS}
#'   \item{PFNA}{log10 concentration in PFNA}
#'   \item{PFOA}{log10 concentration in PFOA}
#'   \item{PFUnDA}{log10 concentration in PFUnDA}
#'   \item{X13C}{carbon isotopic ratio}
#'   \item{X15N}{nitrogen isotopic ratio}
#' }
#'
#' @references Munoz et al. (2017) Evidence for the Trophic Transfer of Perfluoroalkylated Substances in a Temperate Macrotidal Estuary. Environmental Science & Technology: 51, 8450-8459
#'
#' @source \url{https://pubs.acs.org/doi/abs/10.1021/acs.est.7b02399}
"signature_data"



#' Limit of quantifications
#'
#' Limit of quantifications estimated by Munoz et al. (2017). The number of lines
#' is similar to the number of lines of \code{\link{signature_data}} and gives for
#' each individual and each tracer the limit of quantifications of the measures.
#'
#' @format A data frame with 53940 rows and 10 variables:
#' \describe{
#'   \item{FOSA}{log10 LOQ FOSA}
#'   \item{L.PFOS}{log10 LOQ L.PFOS}
#'   \item{PFNA}{log10 LOQ in PFNA}
#'   \item{PFOA}{log10 LOQ PFOA}
#'   \item{PFUnDA}{log10 LOQ PFUnDA}
#'   \item{X13C}{LOQ is meaningless for Carbon so we put a very low value}
#'   \item{X15N}{LOQ is meaningless for Nitrogen so we put a very low value}
#' }
#'
#' #' @references Munoz et al. (2017) Evidence for the Trophic Transfer of Perfluoroalkylated Substances in a Temperate Macrotidal Estuary. Environmental Science & Technology: 51, 8450-8459
#'
#' @source \url{https://pubs.acs.org/doi/abs/10.1021/acs.est.7b02399}
"LOQ"




#' Measurement from previous studies
#'
#' This table summarizes measurements made by David (2006). Each line gives
#' the mean, standard-deviation and number of samples from measurements made
#' on a give tracer. Values should be log (or log10) transformed if required,
#' similarly to data in \code{\link{signature_data}}
#'
#' @format A data frame with 53940 rows and 10 variables:
#' \describe{
#'   \item{Species}{name of the species}
#'   \item{tracer}{name of the tracer, contaminant of isotopic ratio}
#'   \item{mean}{mean value observed on the sample (no transformation here)}
#'   \item{sd}{standard value observed on the sample (no transformation here)}
#'   \item{n}{number of samples used to estimate mean and sd}
#' }
#'
#' @references David (2006) Dynamique spatio-temporelle du zooplancton dans l'estuaire de la Gironde et implications au sein du réseau trophique planctonique. Ph.D. Thesis, Univ. de Bordeaux 1, 313pp
"prior_signature_data"


#' Prior diet matrix
#'
#' This square matrix is used to provide indicate a predator (in line) can eat a prey (in column).
#' 0 indicates that this is not possible.
#' This matrix came from Ballutaud et al. (in press) and build from analysis
#' made by Pasquaud et al. (2010)
#'
#' @format a square matrix with
#'
#' @references Pasqaud et al. (2010) Determination of fish trophic levels in an estuarine system. Estuarine, Coastal and Shelf Scienc: 86, 237-246
#' @references Ballutaud et al. (in press) EStimating Contaminants tRansfers Over Complex food webs (ESCROC): an innovative Bayesian method for estimating POP's biomagnification in aquatic food webs. Science of the Total Environment
"prior_diet_matrix"



