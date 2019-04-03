context("Test results")

test_that("Test invalid input", {

    gsd = makeExample()
    expect_error(results(gsd),
        "NBAMSeq function should be run before calling results function")
    gsd = NBAMSeq(gsd, parallel = TRUE)
    expect_error(results(gsd, name = 1))
    expect_error(results(gsd, name = c("pheno","var1")))
    expect_error(results(gsd),
                 "Either 'name' or 'contrast' argument should be provided.")
    expect_error(results(gsd, name = "var10"))

    expect_error(results(gsd, name = "pheno", indepfilter = "unknown"))
    expect_error(results(gsd, name = "pheno", indepfilter = c(TRUE, FALSE)))

    expect_error(results(gsd, name = "pheno", alpha = "unknown"))
    expect_error(results(gsd, name = "pheno", alpha = c(0.5, 0.6)))
    expect_error(results(gsd, name = "pheno", alpha = -0.1))
    expect_error(results(gsd, name = "pheno", alpha = 1.2))
})

test_that("Test result output", {
    n = 200
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
    expect_error(results(gsd, name = "var4"),
                 " 'name' should be a continuous variable. For factors,
                please use 'contrast' argument.")
    expect_error(results(gsd, contrast = c("var4", 1)),
                 " 'contrast' should be a character of length 3.")
    expect_error(results(gsd, contrast = c("var4", 1, 1)),
                 "2nd and 3rd element in constrast should be different.")
    expect_error(results(gsd, contrast = c("var10", 1, 0)),
                 "1st element in contrast should be a variable in colData.")
    expect_error(results(gsd, contrast = c("var2", 1, 0)),
                 "The variable in contrast should be a factor. For continuous
                variables, please use 'name' argument.")
    expect_error(results(gsd, contrast = c("var4", 5, 0)),
                 "2nd element in contrast should be an appropriate level.")
    expect_error(results(gsd, contrast = c("var4", 1, 10)),
                 "3rd element in contrast should be an appropriate level.")

    res1 = results(gsd, name = "pheno")
    expect_true("baseMean"%in%names(res1))
    expect_true("edf"%in%names(res1))
    expect_true("stat"%in%names(res1))
    expect_true("pvalue"%in%names(res1))
    expect_true("padj"%in%names(res1))

    res2 = results(gsd, name = "var3")
    expect_true("baseMean"%in%names(res2))
    expect_true("coef"%in%names(res2))
    expect_true("SE"%in%names(res2))
    expect_true("stat"%in%names(res2))
    expect_true("pvalue"%in%names(res2))
    expect_true("padj"%in%names(res2))

    ## check contrast
    res3 = results(gsd, contrast = c("var4", 1, 2), parallel = TRUE)
    res4 = results(gsd, contrast = c("var4", 2, 1), parallel = TRUE)
    expect_identical(res3[["baseMean"]], res4[["baseMean"]])
    expect_equal(res3[["coef"]],-res4[["coef"]], tolerance=1e-8)
    expect_equal(res3[["SE"]],res4[["SE"]], tolerance=1e-8)
    expect_equal(res3[["stat"]],-res4[["stat"]], tolerance=1e-8)
    expect_equal(res3[["pvalue"]],res4[["pvalue"]], tolerance=1e-8)
    expect_equal(res3[["padj"]],res4[["padj"]], tolerance=1e-8)

    ## check parallel option:
    res7 = results(gsd, contrast = c("var4", 2, 1), parallel = FALSE)
    expect_identical(res4[["baseMean"]], res7[["baseMean"]])
    expect_equal(res4[["coef"]], res7[["coef"]], tolerance=1e-8)
    expect_equal(res4[["SE"]], res7[["SE"]], tolerance=1e-8)
    expect_equal(res4[["stat"]], res7[["stat"]], tolerance=1e-8)
    expect_equal(res4[["pvalue"]], res7[["pvalue"]], tolerance=1e-8)
    expect_equal(res4[["padj"]], res7[["padj"]], tolerance=1e-8)


    res5 = results(gsd, contrast = c("var4", 1, 0))
    res6 = results(gsd, contrast = c("var4", 0, 1))
    expect_identical(res5[["baseMean"]], res6[["baseMean"]])
    expect_identical(res5[["coef"]], -res6[["coef"]])
    expect_identical(res5[["SE"]], res6[["SE"]])
    expect_identical(res5[["stat"]], -res6[["stat"]])
    expect_equal(res5[["pvalue"]],res6[["pvalue"]], tolerance=1e-8)
    expect_equal(res5[["padj"]],res6[["padj"]], tolerance=1e-8)

    expect_true(all(res5[["coef"]]==mcols(gsd)[["var4_1_vs_0"]]))
    expect_true(all(res5[["SE"]]==mcols(gsd)[["SE_var4_1_vs_0"]]))
    expect_true(all(res5[["pvalue"]]==mcols(gsd)[["PValue_var4_1_vs_0"]]))


})





