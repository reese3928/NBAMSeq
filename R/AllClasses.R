
#' NBAMSeqDataSet class
#'
#' \code{NBAMSeqDataSet} is a class inherited from
#' \code{\link{SummarizedExperiment}}.
#' It is used to store the count matrix, colData, and design formula
#' in differential expression analysis.
#' @slot design a mgcv-type design formula
#' @references Martin Morgan, Valerie Obenchain, Jim Hester and 
#' Hervé Pagès (2018). SummarizedExperiment: SummarizedExperiment container. 
#' R package version 1.12.0.
#' @export NBAMSeqDataSet
setClass("NBAMSeqDataSet",
    contains="SummarizedExperiment",
    representation=representation(
        design = "ANY"
    )
)

setValidity("NBAMSeqDataSet", function(object){
    ##########################################################################
    # The following validity checks are taken from DESeq2 package.
    # Author: Michael Love (michaelisaiahlove@gmail.com)
    # Reference: Love, M.I., Huber, W., Anders, S. (2014) Moderated estimation
    # of fold change and dispersion for RNA-seq data with DESeq2. Genome
    # Biology, 15:550.

    if (! ("counts" %in% assayNames(object)) )
        return( "the assays slot must contain a matrix named 'counts'" )
    if ( !is.numeric( assay(object) ) )
        return( "the count data is not numeric" )
    if ( any( is.na( assay(object) ) ) )
        return( "NA values are not allowed in the count matrix" )
    if ( any( assay(object) < 0 ) )
        return( "the count data contains negative values" )

    design <- getDesign(object)
    stopifnot(is(design, "formula"))

    designVars <- all.vars(design)
    if (!all(designVars %in% names(colData(object)))) {
        return("all variables in design formula must be columns in colData")
    }
    designVarsClass <- vapply(designVars, function(v)
        class(colData(object)[[v]]), "numeric")
    if (any(designVarsClass == "character")) {
        return("variables in design formula are character vectors.
convert these columns of colData(object) to factors before including in the
design formula")
    }
    designFactors <- designVars[designVarsClass == "factor"]
    # levels would duplicate after make.names()
    if (any(vapply(designFactors,function(v) {
        factor.lvls <- levels(colData(object)[[v]])
        factor.nms <- make.names(factor.lvls)
        any(duplicated(factor.nms))
    }, FALSE))) {
        return("levels of factors in the design have non-unique level
names after make.names() is applied. best to only use letters and numbers
for levels of factors in the design")
    }
    # levels contain characters other than letters, numbers, and underscore
    if (any(vapply(designFactors,function(v) {
        factor.lvls <- levels(colData(object)[[v]])
        any(!grepl("^[A-Za-z0-9_.]+$",factor.lvls))
    }, FALSE))) {
        # just a warning for now
        message("Note: levels of factors in the design contain characters
other than letters, numbers, '_' and '.'. It is recommended (but not required)
to use only letters, numbers, and delimiters '_' or '.', as these are safe
characters for column names in R.
[This is a message, not an warning or error]")
    }

    ## end: DESeq2 checks
    ##########################################################################

    ## extra checks for design
    if(!grepl("s\\(", as.character(design)[2])){
        return("At least one variable in design should be nonlinear.
If all variables are linear, use DESeq2 instead.")
    }
    sterms = vapply(designVars, function(x)
        grepl(paste0("s\\(",x,"\\)"), as.character(design)[2]), FALSE)
    if(any(sterms&designVarsClass == "factor")){
        return("Non linear term can not be a factor.")
    }

    ## extra checks for counts
    if(all(assay(object)==0)){
        return("countData can not be all 0.")
    }

    TRUE
})


#' NBAMSeqDataSet constructor
#' @param countData a matrix or data frame contains gene count
#' @param colData a \code{DataFrame} or \code{data.frame}
#' @param design a mgcv type design. e.g. \code{~ s(pheno)} or
#' \code{~ s(pheno) + var1 + var2}
#' @param ... optional arguments passed to \code{SummarizedExperiment}
#' @return a NBAMSeqDataSet object
#' @examples
#' n = 100  ## n stands for number of genes
#' m = 20   ## m stands for sample size
#' countData = matrix(rnbinom(n*m, mu=100, size=1/3), ncol = m)
#' mode(countData) = "integer"
#' colnames(countData) = paste0("sample", 1:m)
#' rownames(countData) = paste0("gene", 1:n)
#' pheno = runif(m, 20, 80)
#' colData = data.frame(pheno = pheno)
#' rownames(colData) = paste0("sample", 1:m)
#' gsd = NBAMSeqDataSet(countData = countData,
#' colData = colData, design = ~s(pheno))
#'
NBAMSeqDataSet <- function(countData, colData, design, ...){

    countData = as.matrix(countData)
    colData = DataFrame(colData)

    ##########################################################################
    # The following validity checks are taken from DESeq2 package.
    # Author: Michael Love (michaelisaiahlove@gmail.com)
    # Reference: Love, M.I., Huber, W., Anders, S. (2014) Moderated estimation
    # of fold change and dispersion for RNA-seq data with DESeq2. Genome
    # Biology, 15:550.

    # check that these agree in number
    stopifnot(ncol(countData) == nrow(colData))
    # check if the rownames of colData are simply in different order
    # than the colnames of the countData, if so throw an error
    # as the user probably should investigate what's wrong
    if (!is.null(rownames(colData)) & !is.null(colnames(countData))) {
        if (all(sort(rownames(colData)) == sort(colnames(countData)))) {
            if (!all(rownames(colData) == colnames(countData))) {
                stop(paste("rownames of the colData:
                    ",paste(rownames(colData),collapse=","),"
                    are not in the same order as the colnames of the countData:
                    ",paste(colnames(countData),collapse=",")))
            }
        }
    }
    if (is.null(rownames(colData)) & !is.null(colnames(countData))) {
        rownames(colData) <- colnames(countData)
    }
    ## end: DESeq2 checks
    ##########################################################################

    se = SummarizedExperiment(assays = list(counts = countData),
        colData = colData, ...)
    object = new("NBAMSeqDataSet", se, design = design)

    if(any(rowSums(assay(object))==0)){
        message(sum(rowSums(assay(object))==0),
" genes have all 0 counts for all samples. Consider filtering out these genes
before differential expression analysis." )
    }
    if(any(assay(object)!=round(assay(object)))){
        message("countData contains non-integers, rounded to the nearest
                integer automatically.")
        assay(object) = round(assay(object))
    }
    metadata(object)$fitted = FALSE
    object
}


#' Accessor functions and replace methods for NBAMSeqDataSet object
#' @name NBAMSeq-methods
#' @rdname NBAMSeq-methods
#' @param theObject a NBAMSeqDataSet object
#' @param value the values to be included in the object
NULL

#' @rdname NBAMSeq-methods
#' @export
setGeneric("getDesign", function(theObject) standardGeneric("getDesign"))

#' For \code{getDesign()}: accessor to the design formula
#' @rdname NBAMSeq-methods
#' @examples
#' ## For getDesign() ##
#' gsd = makeExample()
#' design_gsd = getDesign(gsd)
#' @export
#' @return For \code{getDesign()}: design formula
setMethod("getDesign", "NBAMSeqDataSet", function(theObject) theObject@design)


#' @rdname NBAMSeq-methods
#' @export
setGeneric("getsf", function(theObject) standardGeneric("getsf"))

#' For \code{getsf()}: accessor to the size factors
#' @rdname NBAMSeq-methods
#' @examples
#' ## For getsf() ##
#' gsd = makeExample()
#' sf = getsf(gsd)
#' @references Love, M.I., Huber, W., Anders, S. (2014) Moderated estimation of
#' fold change and dispersion for RNA-seq data with DESeq2. Genome Biology,
#' 15:550. \url{https://doi.org/10.1186/s13059-014-0550-8}
#' @return For \code{getsf()}: size factor
#' @export
setMethod("getsf", "NBAMSeqDataSet", function(theObject){
    if (!"sizeFactors" %in% names(colData(theObject)))  return(NULL)
    sf = theObject$sizeFactors
    names(sf) = colnames(theObject)
    sf
})


#' Replace methods for NBAMSeqDataSet object
#' @rdname NBAMSeq-methods
#' @export
setGeneric("setsf<-", function(theObject, value) standardGeneric("setsf<-"))

#' For \code{setsf()}: replace size factors
#' @rdname NBAMSeq-methods
#' @examples
#' ## For setsf() ##
#' n = 100
#' m = 50
#' gsd = makeExample(n = n, m = m)
#' sf = sample(1:5, m, replace = TRUE)
#' setsf(gsd) = sf
#' @return For \code{setsf()}: NBAMSeq object
#' @exportMethod "setsf<-"
setReplaceMethod("setsf", signature = c("NBAMSeqDataSet", "numeric"),
    function(theObject, value){
        stopifnot(is.numeric(value))
        stopifnot(all(!is.na(value)))
        if(length(value)!=ncol(theObject)){
            stop("Size factor length should be the same as number of samples.")
        }
        if(any(value<=0)){
            stop("Size factor cannot be 0 or negative.")
        }
        colData(theObject)$sizeFactors = value
        theObject
})





