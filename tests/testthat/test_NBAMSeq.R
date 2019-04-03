context("Test NBAMSeq")

library(DESeq2)
library(mgcv)

test_that("Test invalid input", {

    nrows = 200
    ncols = 10
    counts = matrix(rnbinom(nrows * ncols, 10, 0.3), nrows)
    pheno = runif(ncols, 20, 80)
    colData = DataFrame(pheno = pheno)
    gsd = NBAMSeqDataSet(countData = counts,
        colData = colData, design = ~s(pheno))
    expect_error(NBAMSeq(gsd, gamma = "unknown"))
    expect_error(NBAMSeq(gsd, gamma = -2),
                 " 'gamma' should be greater or equal to 1.")
    expect_error(NBAMSeq(gsd, gamma = c(1,2)))
    expect_error(NBAMSeq(gsd, parallel = 1))
    expect_error(NBAMSeq(gsd, parallel = c(TRUE, FALSE)))

    dds = makeExampleDESeqDataSet()
    expect_error(NBAMSeq(dds))

})

test_that("Test parallel option", {
    gsd = makeExample(n=30,m=30)
    gsd1 = NBAMSeq(gsd, parallel = TRUE)
    gsd2 = NBAMSeq(gsd, parallel = FALSE)
    expect_identical(gsd1,gsd2)


})

test_that("Test NBAMSeq output", {
    gsd = makeExample(n=30,m=30)
    gsd = NBAMSeq(gsd, parallel = TRUE)
    expect_true("Intercept"%in%names(mcols(gsd)))
    expect_true("edf_pheno"%in%names(mcols(gsd)))
    expect_true("Chisq_pheno"%in%names(mcols(gsd)))
    expect_true("PValue_pheno"%in%names(mcols(gsd)))
    expect_true("df_residual"%in%names(mcols(gsd)))
    expect_true("null_deviance"%in%names(mcols(gsd)))
    expect_true("df_null"%in%names(mcols(gsd)))

    n = 30
    m = 30
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
    colnames(countData) = paste0("sample", 1:m)
    rownames(countData) = paste0("gene", 1:n)

    var1 = rnorm(m)
    var2 = rnorm(m)
    var3 = rnorm(m)
    var4 = as.factor(sample(c(0,1,2), m, replace = TRUE))
    colData = data.frame(pheno = pheno, var1 = var1, var2 = var2,
                         var3 = var3, var4 = var4)
    gsd = NBAMSeqDataSet(countData = countData, colData = colData,
                        design = ~ var4 + s(pheno) + s(var1) + var2 + var3)

    gsd = NBAMSeq(gsd, parallel = TRUE)
    expect_true("Intercept"%in%names(mcols(gsd)))
    expect_true("edf_var1"%in%names(mcols(gsd)))
    expect_true("Chisq_var1"%in%names(mcols(gsd)))
    expect_true("PValue_var1"%in%names(mcols(gsd)))
    expect_true("var2"%in%names(mcols(gsd)))
    expect_true("SE_var2"%in%names(mcols(gsd)))
    expect_true("PValue_var2"%in%names(mcols(gsd)))
    expect_true("smooth_pheno"%in%names(mcols(gsd)))
    expect_true("smooth_var1"%in%names(mcols(gsd)))
    expect_true("df_residual"%in%names(mcols(gsd)))
    expect_true("null_deviance"%in%names(mcols(gsd)))
    expect_true("df_null"%in%names(mcols(gsd)))


})

test_that("Test additional argument for gam", {
    gsd = makeExample(n=30,m=30)
    gsd = NBAMSeq(gsd, parallel = TRUE, control=gam.control(maxit = 1))
    expect_true(all(mcols(gsd)[["innerIter"]]<=1))
})



# test_that("Test ncols(dds) ", {
#    gsd = makeExample()
#    ddsdesign = formula(paste0("~",
#                               paste(all.vars(getDesign(gsd)), collapse= "+")))
#    dds = DESeqDataSetFromMatrix(countData = assay(gsd), colData = colData(gsd),
#                                 design =  ddsdesign)
#    dds = estimateSizeFactors(dds)
#    logsf = log(sizeFactors(dds))
#    gamma = log2(ncol(gsd))

#    dat = data.frame(colData(gsd))
#    dat$logsf = logsf

#    formula_offset = update(getDesign(gsd), y ~ . + offset(logsf))

#    gamFit1 = function(i){
#        dat$y = assay(gsd)[i,]
#        gamfit = gam(formula_offset, family = nb(link = "log"),
#                     gamma = gamma, data = dat)
#        mcols(dds)$dispGeneEst[i] = 1/(gamfit$family$getTheta(TRUE))
#        list(theta = gamfit$family$getTheta(TRUE), sp = gamfit$sp,
#             coef = coef(gamfit), muhat = gamfit$fitted.values)
#    }
#    gamGeneEst = lapply(1:nrow(gsd), gamFit1)
#    expect_true(all(mcols(dds)$dispGeneEst ==
#                        1/vapply(gamGeneEst, function(x) x$theta, 1)))

#})





