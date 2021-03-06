% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spe_meth.R
\name{ts_Counts}
\alias{ts_Counts}
\title{Counts}
\usage{
ts_Counts(x, cutoff)
}
\arguments{
\item{x}{A numeric vector.}

\item{cutoff}{numeric. values under \code{cutoff} will be set to
0(unexpression).}
}
\description{
Function require a vector with expression of one gene in different tissues.
If expression for one tissue is not known, gene specificity for this gene is
NA. Function requires setting of a \code{cutoff} (rpkm/FPKM/TPM).
Transformation: (1-x/N)*(N/(N-1)).
}
\examples{
\dontrun{
ts_Counts(c(0,1,2,3,2,8), cutoff = 1)
}
}
