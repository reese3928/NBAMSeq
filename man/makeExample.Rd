% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/makeExample.R
\name{makeExample}
\alias{makeExample}
\title{Make an example NBAMSeqDataSet}
\usage{
makeExample(n = 200, m = 30)
}
\arguments{
\item{n}{number of genes}

\item{m}{number of samples}
}
\value{
a NBAMSeqDataSet object
}
\description{
This function makes an example NBAMSeqDataSet
}
\examples{
gsd = makeExample()
}
\references{
Love, M.I., Huber, W., Anders, S. (2014) Moderated estimation of
fold change and dispersion for RNA-seq data with DESeq2. Genome Biology,
15:550. \url{https://doi.org/10.1186/s13059-014-0550-8}
}
