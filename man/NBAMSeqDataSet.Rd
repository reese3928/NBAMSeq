% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllClasses.R
\name{NBAMSeqDataSet}
\alias{NBAMSeqDataSet}
\title{NBAMSeqDataSet constructor}
\usage{
NBAMSeqDataSet(countData, colData, design, ...)
}
\arguments{
\item{countData}{a matrix or data frame contains gene count}

\item{colData}{a \code{DataFrame} or \code{data.frame}}

\item{design}{a mgcv type design. e.g. \code{~ s(pheno)} or
\code{~ s(pheno) + var1 + var2}}

\item{...}{optional arguments passed to \code{SummarizedExperiment}}
}
\value{
a NBAMSeqDataSet object
}
\description{
NBAMSeqDataSet constructor
}
\examples{
n = 100  ## n stands for number of genes
m = 20   ## m stands for sample size
countData = matrix(rnbinom(n*m, mu=100, size=1/3), ncol = m)
mode(countData) = "integer"
colnames(countData) = paste0("sample", 1:m)
rownames(countData) = paste0("gene", 1:n)
pheno = runif(m, 20, 80)
colData = data.frame(pheno = pheno)
rownames(colData) = paste0("sample", 1:m)
gsd = NBAMSeqDataSet(countData = countData,
colData = colData, design = ~s(pheno))

}
