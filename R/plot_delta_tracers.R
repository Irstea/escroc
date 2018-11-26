#' Plot Magnification / enrichment
#' Plot the posterior distribution of magnification or enrichment for specified tracers
#' @param myfit a fit as returned by \code{\link{fit_escroc}}
#' @param mydata a list as returned by \code{\link{prepare_data}}
#' @param trac_names a character or vector list containing the tracers
#' to be plotted
#' @param transfo either NULL or a list of functions, on per tracer, to carry transformation if required
#'
#' @return a list of the same length as trac_names with a \code{\link[ggplot2]{ggplot}}
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 geom_density
#' @importFrom coda varnames
#' @examples
#' #importing data
#' data(prior_diet_matrix)
#' data(LOQ)
#' data(prior_signature_data)
#' prior_magnification <- data.frame(tracer=c("X15N","X13C"),mean=c(3,0),sd=c(1,1))
#'
#' #check that everything is ok
#' mydata <- prepare_data(prior_diet_matrix,signature_data,prior_signature_data,isLeftCensored,
#' LOQ,prior_magnification)
#'
#' #build the model
#' mymodel <- building_model(mydata)
#'
#' #fit the model
#' myresults <- fit_escroc(mydata, mymodel)
#'
#' #a function to compute TMF ()
#' transfo <- function (x) 10^x
#' #plot diet
#' print(plot_delta_tracers(myresults,mydata,"FOSA",list(transfo)))
#' @export
plot_delta_tracers <- function(myfit, mydata, trac_names,transfo=list(NULL)) {

  mapply(function(trac_name,trans) {
    if (is.na(match(trac_name, colnames(mydata$signature_data))))
      stop("the tracer name does not match tracers from signature_data")
    column_names <- paste("delta\\[", trac_name, sep = "")
    subset_fit <-
      as.matrix(myfit[, grep(column_names, varnames(myfit)), drop = FALSE])
    if (!is.null(trans)) subset_fit[,1] <- trans(subset_fit[,1])
    colnames(subset_fit) <-
      sapply(colnames(subset_fit), function(old_name) {
        strsplit(substr(old_name, 1, nchar(old_name) - 1), c("[\\[,]"))[[1]][2]
      })
    ggplot(as.data.frame(subset_fit),aes_string(trac_name)) +
      geom_density()
  },trac_names,transfo,SIMPLIFY=FALSE)

}
