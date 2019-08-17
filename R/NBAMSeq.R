#' @title Differential expression analysis based on negative binomial additive
#' model
#'
#' @description This function performs differential expression analysis based
#' on negative binomial additive model.
#' @param object a NBAMSeqDataSet object
#' @param gamma a number greater or equal to 1. Increase gamma to create
#' smoother models. Default gamma is 2. See \code{\link[mgcv]{gam}} for
#' details.
#' @param parallel either TRUE or FALSE indicating whether parallel should be
#' used. Default is FALSE
#' @param BPPARAM an argument provided to \code{\link{bplapply}}. See
#' \code{\link[BiocParallel]{register}} for details.
#' @param ... additional arguments provided to \code{\link[mgcv]{gam}}
#' @export
#' @importFrom mgcv gam s nb negbin
#' @import SummarizedExperiment S4Vectors DESeq2
#' @importFrom BiocParallel bplapply bpparam
#' @importFrom methods is new
#' @importFrom stats coef formula update
#' @return a NBAMSeqDataSet object
#' @references Love, M.I., Huber, W., Anders, S. (2014) Moderated estimation of
#' fold change and dispersion for RNA-seq data with DESeq2. Genome Biology,
#' 15:550. \url{https://doi.org/10.1186/s13059-014-0550-8}
#' @examples
#' gsd = makeExample(n = 3, m = 10)
#' gsd = NBAMSeq(gsd)

NBAMSeq <- function(object, gamma = 2, parallel = FALSE,
    BPPARAM = bpparam(), ...){

    ## check input
    stopifnot(is(object, "NBAMSeqDataSet"))
    if(!is.null(gamma)){
        stopifnot(is.numeric(gamma))
        stopifnot(length(gamma)==1)
        if(gamma < 1){
            stop(" 'gamma' should be greater or equal to 1.")
        }
    }
    stopifnot(is.logical(parallel))
    stopifnot(length(parallel)==1)

    ## construct a DESeqDataSet object
    ddsdesign = formula(paste0("~",
        paste(all.vars(getDesign(object)), collapse= "+")) )
    dds = DESeqDataSetFromMatrix(countData = assay(object),
        colData = colData(object), design =  ddsdesign)

    if("sizeFactors"%in%names(colData(dds))){
        logsf = log(colData(object)$sizeFactors)
        object$logsf = logsf    ## save logsf in object
        sizeFactors(dds) = colData(object)$sizeFactors
    }else{
        ## estimate size factors
        dds = estimateSizeFactors(dds)
        logsf = log(sizeFactors(dds))
        colData(object)$sizeFactors = sizeFactors(dds)
        ## save logsf in object
        object$logsf = logsf
    }

    dat = data.frame(colData(object))
    dat$logsf = logsf
    formula_offset = update(getDesign(object), y ~ . + offset(logsf))

    ## gene-wise GAM model
    gamFit1 = function(i){
        dat$y = assay(object)[i,]
        gamfit = tryCatch(
            expr = {
                tmp = gam(formula_offset, family = nb(link = "log"),
                    method = "REML", gamma = gamma, data = dat, ...)
                tmp$gamma = gamma  ## save gamma value in gam object
                tmp
            },
            error = function(e){
                tmp = gam(formula_offset, family = nb(link = "log"),
                    method = "REML", gamma = 1, data = dat, ...)
                tmp$gamma = 1
                tmp
            },
            warning = function(w){
                tmp = gam(formula_offset, family = nb(link = "log"),
                    method = "REML", gamma = 1, data = dat, ...)
                tmp$gamma = 1
                tmp
            }
        )
        list(theta = gamfit$family$getTheta(TRUE), sp = gamfit$sp,
            coef = coef(gamfit), muhat = gamfit$fitted.values,
            outIter = gamfit$outer.info$iter, gamma = gamfit$gamma)
    }

    message("Estimating smoothing parameters and gene-wise dispersions")
    if(parallel){
        gamGeneEst = bplapply(seq_len(nrow(object)),gamFit1,BPPARAM = BPPARAM)
    } else {
        gamGeneEst = lapply(seq_len(nrow(object)), gamFit1)
    }

    ##  get gam gene wise dispersion estimates and save them in mcols(dds)
    mcols(dds)$dispGeneEst = 1/vapply(gamGeneEst, function(x) x$theta, 1)
    ##  bound gene wise dispersion
    maxDisp = pmax(10, ncol(dds))
    mcols(dds)$dispGeneEst = pmin(mcols(dds)$dispGeneEst, maxDisp)

    ##  fit dispersion trend via estimateDispersionsFit function in DESeq2
    message("Estimating dispersion trend")
    dds = estimateDispersionsFit(dds)

    ##  get gam mu estimates and save them in assays(dds)[["mu"]]
    muhat = t(vapply(gamGeneEst, function(x) x$muhat, rep(1, ncol(dds))))
    colnames(muhat) = colnames(dds)
    rownames(muhat) = rownames(dds)
    assays(dds)[["mu"]] = muhat

    ##  MAP dispersion estimates
    message("Estimating MAP dispersion")
    dds = tryCatch(
        expr = {
            estimateDispersionsMAP(dds)
        },
        error = function(e){
            ## avoid possible matrix singular error in DESeq2 C++ code
            assays(dds)[["mu"]] = muhat+1e-6
            estimateDispersionsMAP(dds)
        }
    )
    ind = which(is.na(mcols(dds)$dispMAP))
    if(length(ind)>0){
        mcols(dds)$dispMAP[ind] = mcols(dds)$dispGeneEst[ind]
    }
    
    gamDispMAP = mcols(dds)$dispMAP
    mcols(object) = mcols(dds)

    gamFit2 = function(i){
        dat$y = assay(object)[i,]   ## ith gene count
        start = gamGeneEst[[i]]$coef   ## initial coefficients
        gamFinalFit = gam(formula_offset,
            family = negbin(theta = 1/gamDispMAP[i],link = "log"),
            method = "REML", sp = gamGeneEst[[i]]$sp,
            start = start, data = dat, ...)
        list(paramcoef = summary(gamFinalFit)$p.table[,"Estimate"],
            paramSE = summary(gamFinalFit)$p.table[,"Std. Error"],
            paramZscore = summary(gamFinalFit)$p.table[,"z value"],
            paramPvalue = summary(gamFinalFit)$p.table[,"Pr(>|z|)"],
            smoothedf = summary(gamFinalFit)$s.table[,"edf"],
            smoothChisq = summary(gamFinalFit)$s.table[,"Chi.sq"],
            smoothPvalue = summary(gamFinalFit)$s.table[,"p-value"],
            deviance = gamFinalFit$deviance,
            innerIter = gamFinalFit$iter,  # number of inner iterations
            converged = gamFinalFit$converged,
            residualdf = gamFinalFit$df.residual,
            nulldeviance = gamFinalFit$null.deviance,
            nulldf = gamFinalFit$df.null
            )
    }

    message("Estimating model coefficients")
    if(parallel){
        gamFinal = bplapply(seq_len(nrow(object)), gamFit2, BPPARAM = BPPARAM)
    } else{
        gamFinal = lapply(seq_len(nrow(object)), gamFit2)
    }

    ## process variables
    dat$y = assay(object)[1,]
    gammodel = gam(formula_offset, family = negbin(theta = 3, link = "log"),
        method = "REML", data = dat, fit = FALSE)
    ## process factors
    pterms = vapply(attr(gammodel$pterms,"term.labels"),
                    function(x) is.factor(colData(object)[[x]]), FALSE)
    if(length(pterms)==0){
        pterms = "Intercept"
    } else{
        pterms_name = lapply(seq_along(pterms), function(i){
            if(pterms[i]){
                lv = levels(colData(object)[[names(pterms)[i]]])
                rt = paste0(names(pterms)[i],"_", lv[2:length(lv)],"_vs_",lv[1])
            } else{
                rt = names(pterms)[i]
            }
            rt
        })
        pterms = c("Intercept", unlist(pterms_name))
    }

    sterms = setdiff(attr(gammodel$terms,"term.labels"),
        attr(gammodel$pterms,"term.labels"))
    stopifnot(length(sterms)>=1)

    ## save results in mcols(object)
    n1 = length(gamFinal[[1]]$paramcoef)  ## number of parametric variables
    n2 = length(gamFinal[[1]]$smoothedf)  ## number of nonparametric variables
    nm = names(mcols(object))
    mcols(object) = cbind(mcols(object),
        DataFrame(t(rbind(
            vapply(gamFinal, function(x) x$paramcoef, rep(0, n1)),
            vapply(gamFinal, function(x) x$paramSE, rep(1, n1)),
            vapply(gamFinal, function(x) x$paramPvalue, rep(1, n1)),
            vapply(gamFinal, function(x) x$smoothedf, rep(1, n2)),
            vapply(gamFinal, function(x) x$smoothChisq, rep(1, n2)),
            vapply(gamFinal, function(x) x$smoothPvalue, rep(1, n2)),
            vapply(gamFinal, function(x) x$deviance, 1),
            vapply(gamGeneEst, function(x) x$outIter, 1),
            vapply(gamFinal, function(x) x$innerIter, 1),
            vapply(gamFinal, function(x) x$converged, TRUE),
            vapply(gamGeneEst, function(x) x$sp, rep(1, n2)),
            vapply(gamFinal, function(x) x$residualdf, 0),
            vapply(gamFinal, function(x) x$nulldeviance, 0),
            vapply(gamFinal, function(x) x$nulldf, 0),
            vapply(gamGeneEst, function(x) x$gamma, 1)
            ))
            )
        )

    nm = c(nm, pterms, paste0("SE_",pterms), paste0("PValue_", pterms),
        paste0("edf_", sterms), paste0("Chisq_", sterms),
        paste0("PValue_", sterms), "deviance", "outIter", "innerIter",
        "converged", paste0("smooth_", sterms),
        "df_residual", "null_deviance", "df_null", "gamma")
    colnames(mcols(object)) = nm
    class(mcols(object)[["outIter"]]) = "integer"
    class(mcols(object)[["innerIter"]]) = "integer"
    class(mcols(object)[["converged"]]) = "logical"

    metadata(object)$fitted = TRUE
    message("Done!")
    object
}



