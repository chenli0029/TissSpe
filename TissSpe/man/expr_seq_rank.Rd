% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/binary_pattern.R
\name{expr_seq_rank}
\alias{expr_seq_rank}
\title{Ranking gene expression according to equal-width-intervals}
\usage{
expr_seq_rank(df, n = 12, min = 0, step = 1)
}
\arguments{
\item{df}{data.frame.}

\item{n}{number of medium intervals to generated (all \code{n}+2). Default 12.}

\item{min}{psi under \code{min} will be set to 0. Default 0.}

\item{step}{interval width. Default 1, LFC: log2(2) = 1.}
}
\value{
data.frame contain ranks from \code{df}.
}
\description{
Function to cut data.frame range into n equal-width-interval points, maximal
point and minimal point, then rank them for their value. This function is
design for processed data, like log2 transformed gene expression
value(RPKM/FPKM/TPM).
}
