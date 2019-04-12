
#' @title Make an example NBAMSeqDataSet
#'
#' @description This function makes an example NBAMSeqDataSet
#' @param n number of genes
#' @param m number of samples
#' @importFrom stats rnbinom rnorm runif
#' @export
#' @return a NBAMSeqDataSet object
#' @references Love, M.I., Huber, W., Anders, S. (2014) Moderated estimation of
#' fold change and dispersion for RNA-seq data with DESeq2. Genome Biology,
#' 15:550. \url{https://doi.org/10.1186/s13059-014-0550-8}
#' @examples
#' gsd = makeExample()

makeExample <- function(n = 200, m = 30){
    pheno = runif(m, 20, 80)
    mu = matrix(rep(NA, n*m), nrow = n)
    interceptMean = 3
    interceptSD = 2
    betamat = matrix(rep(NA, n*3), nrow = n)
    betamat[,1] = rnorm(n, 3, 2)
    betamat[,2] = 0.2
    betamat[,3] = -0.0025

    x = cbind(rep(1,m), pheno, pheno^2)
    mu = t(2^(x %*% t(betamat))) + 100
    mumean = apply(mu, 1, mean)
    disp = 10/mumean + 0.001

    countData = matrix(rnbinom(m*n, mu=mu, size=1/disp), ncol=m)
    mode(countData) = "integer"
    colnames(countData) = paste0("sample", seq_len(m))
    rownames(countData) = paste0("gene", seq_len(n))

    colData = data.frame(pheno = pheno)
    NBAMSeqDataSet(countData = countData, colData = colData,design = ~ s(pheno))

}

