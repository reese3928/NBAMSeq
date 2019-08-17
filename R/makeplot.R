
#' @title Making plots to visualize nonlinear associations
#'
#' @description This function makes plots to visualize nonlinear associations.
#' @param object a NBAMSeqDataSet object
#' @param phenoname the name of nonlinear variable to be visualized
#' @param genename the name of gene to be visualized
#' @param ... additional arguments provided to \code{\link[mgcv]{plot.gam}}
#' @importFrom mgcv gam plot.gam
#' @return the plot made by plot.gam() function
#' @export
#' @examples
#' gsd = makeExample(n = 3, m = 10)
#' gsd = NBAMSeq(gsd)
#' makeplot(gsd, "pheno", "gene3", main = "gene10")


makeplot = function(object, phenoname, genename, ...){
    ## check input
    stopifnot(is(object, "NBAMSeqDataSet"))
    ## check whether NBAMSeq has been run before
    if(!metadata(object)$fitted){
        stop("NBAMSeq function should be run before making plots")
    }
    
    stopifnot(is.character(phenoname))
    stopifnot(length(phenoname)==1)
    stopifnot(phenoname%in%all.vars(getDesign(object)))
    ## check the variable is linear or nonlinear
    flag = grepl(paste0("s\\(",phenoname,"\\)"), 
        as.character(getDesign(object))[2])
    if(!flag){
        stop("Visulization can be made for nonlinear variables only")
    }
    
    stopifnot(is.character(genename))
    stopifnot(length(genename)==1)
    if(!genename%in%rownames(assay(object))){
        stop("genename is not valid")
    }
    ind = which(rownames(assay(object))==genename)
    
    formula_offset = update(getDesign(object), y ~ . + offset(logsf))
    gamDispMAP = mcols(object)$dispMAP
    dat = data.frame(colData(object))
    dat$y = assay(object)[ind,]
    gamfit = gam(formula_offset, family = nb(link = "log"), method = "REML", 
        gamma = mcols(object)[ind,"gamma"], data = dat)
    start = gamfit$coef
    ## get smoothing terms
    sterms = setdiff(attr(gamfit$terms,"term.labels"),
        attr(gamfit$pterms,"term.labels"))
    sps = unlist(mcols(object)[ind, paste0("smooth_",sterms)])
    
    gamFinalFit = gam(formula_offset, 
        family = negbin(theta = 1/gamDispMAP[ind],link = "log"),
        method = "REML", sp = sps, start = start, data = dat)
    
    varind = which(sterms==phenoname)
    plot.gam(gamFinalFit, select = varind, ...)
}



