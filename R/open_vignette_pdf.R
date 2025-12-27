#' Open the prebuilt vignette PDF
#'
#' @export
open_dkosim_vignette_pdf <- function() {
  pdf <- system.file("extdata/vignettes/DKOsimR.pdf", package = "DKOsimR")
  if (pdf == "") stop("Prebuilt PDF vignette not found in the installed package.")
  utils::browseURL(pdf)
}
