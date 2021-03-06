% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spe_meth.R
\name{ts_Tau}
\alias{ts_Tau}
\title{Tau}
\usage{
ts_Tau(x)
}
\arguments{
\item{x}{A numeric vector.}
}
\description{
Function require a vector with expression of one gene in different tissues.
If expression for one tissue is not known, gene specificity for this gene is
NA. Minimum 2 tissues.
}
\examples{
\dontrun{
ts_Tau(c(1,2,3,2,8))
}
}
