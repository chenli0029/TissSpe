% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{ts_psi}
\alias{ts_psi}
\title{Calculate specificity of psi}
\usage{
ts_psi(df, n = 10, min = 0, max = 100, tissues, identifier,
  na.del = TRUE, cutoff = 0.5, mingap = 3)
}
\arguments{
\item{df}{data.frame, which contain psi vaules (0-100). One column is the names of
symbols, like AS events id, etc.}

\item{n}{Integer. \code{n+2} equal-width-intervals be generated of all psi
values. Default 10.}

\item{min}{Numeric. The values under \code{min} will be graded into rank
0(unexpressed). Defaut 0.}

\item{max}{Numeric. The values over \code{max} will be graded into rank
\code{n+1}(highest). Default 100.}

\item{tissues}{Vector of charactors, at leat length of 2. Analysed
Tissues' unique identifier, and must keep away from "Tau", "Gini", "Tsi",
"Counts", "Ee", "Pem", "Hg", "Z", "Spm", "Ib", "Type", "Mean" and "Max",
exactly.}

\item{identifier}{Charactor, length of 1. The colname of unique identifier
for row symbols, like "gene_id", "AS_events", etc.}

\item{na.del}{Logical. Whether NAs will be dropped. Default TRUE.}

\item{cutoff}{Numeric. Values under \code{cutoff} will set to 0(unexpressed)
in Specificity method "Counts".}

\item{mingap}{Integer. Minimal gap of generating binary pattern, Please
refer the paper. Default 3.}
}
\value{
List of 3 data.frame. one of them contains psi values with their
specificty values of 9 methods: "Tau", "Gini", "Tsi", "Counts", "Ee", "Pem",
"Hg", "Z", "Spm"(named "raw"), the second contains rank values and binary
index(named "rank"), the third contains binary pattern values and
index binary "Ib"(named "bin").
}
\description{
Function to find tissue-specific AS-events for a given data.frame.
}
\details{
Function to detect tissue-specific AS events, and return a list with 3
data.frames, which contain psi values, rank values, binary pattern
values and their specificty values of 10 methods: "Tau", "Gini", "Tsi", "Counts",
"Ee", "Pem", "Hg", "Z", "Spm", "Ib". Rows with NA will be dropped. Of the binary
pattern, all values betwwen min and max in the data.frame will be graded into n
equal-width-intervals and then assign the rank 0(unexpressed) to \code{n+1}.
}
\examples{
ts_psi(demo_psi,
       tissues = c("sample_A", "sample_B", "sample_C", "sample_D",
                   "sample_E", "sample_F", "sample_G", "sample_H",
                   "sample_I", "sample_J", "sample_K", "sample_L",
                   "sample_M", "sample_N", "sample_O", "sample_P",
                   "sample_Q"),
                   identifier = "AS_events")
}
