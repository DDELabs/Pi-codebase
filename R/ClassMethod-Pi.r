######################################################################
# pNode
######################################################################
#' @title Definition for S3 class \code{pNode}
#' @description \code{pNode} has 7 components: priority, g, SNP, Gene2SNP, nGenes, eGenes and cGenes.
#' @param priority a data frame
#' @param g an 'igraph' object
#' @param SNP a data frame
#' @param Gene2SNP a data frame
#' @param nGenes a data frame
#' @param eGenes a data frame
#' @param cGenes a data frame
#' @return an object of S3 class \code{pNode}
#' @keywords S3 classes
#' @export
#' @examples
#' \dontrun{
#' # Load the library
#' library(Pi)
#' }
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' \dontrun{
#' pNode(evidence, metag)
#' }
pNode <- function(priority, g, SNP, Gene2SNP, nGenes, eGenes, cGenes){
	## integrity checks
	if(!is(priority,'data.frame') | !is(g,'igraph')){
		stop("The S3 class 'pNode' object failed to pass integrity checks!\n")
	}
	value <- list(priority=priority, g=g, SNP=SNP, Gene2SNP=Gene2SNP, nGenes=nGenes, eGenes=eGenes, cGenes=cGenes)
	class(value) <- "pNode"
	return(value)
}
#' @param x an object of class \code{pNode}
#' @param ... other parameters
#' @rdname pNode
#' @export
print.pNode <- function(x, ...) {
	cat(sprintf("An object of S3 class '%s', with %d components:", class(x), length(names(x))), "\n", sep="")
	cat(sprintf("  $g: an igraph object with %d nodes and %d edges", vcount(x$g), ecount(x$g)), "\n", sep="")
	cat(sprintf("  $priority: a data frame of %d rows X %d columns", dim(x$priority)[1],dim(x$priority)[2]), "\n", sep="")
	if(!is.null(x$SNP)){
		cat(sprintf("  $SNP: a data frame of %d rows X %d columns", dim(x$SNP)[1],dim(x$SNP)[2]), "\n", sep="")
	}
	if(!is.null(x$Gene2SNP)){
		cat(sprintf("  $Gene2SNP: a data frame of %d rows X %d columns", dim(x$Gene2SNP)[1],dim(x$Gene2SNP)[2]), "\n", sep="")
	}
	if(!is.null(x$nGenes)){
		cat(sprintf("  $nGenes: a data frame of %d rows X %d columns", dim(x$nGenes)[1],dim(x$nGenes)[2]), "\n", sep="")
	}
	if(!is.null(x$eGenes)){
		cat(sprintf("  $eGenes: a data frame of %d rows X %d columns", dim(x$eGenes)[1],dim(x$eGenes)[2]), "\n", sep="")
	}
	if(!is.null(x$cGenes)){
		cat(sprintf("  $cGenes: a data frame of %d rows X %d columns", dim(x$cGenes)[1],dim(x$cGenes)[2]), "\n", sep="")
	}
	cat("\n--------------------------------------------------\n")
	cat("$priority:\n")
	print(x$priority[1:2,], row.names=FALSE)
	cat("......\n")
	if(!is.null(x$SNP)){
		cat("$SNP:\n")
		print(x$SNP[1:2,], row.names=FALSE)
		cat("......\n")
	}
	if(!is.null(x$Gene2SNP)){
		cat("$Gene2SNP:\n")
		print(x$Gene2SNP[1:2,], row.names=FALSE)
		cat("......\n")
	}
	if(!is.null(x$nGenes)){
		cat("$nGenes:\n")
		print(x$nGenes[1:2,], row.names=FALSE)
		cat("......\n")
	}
	if(!is.null(x$eGenes)){
		cat("$eGenes:\n")
		print(x$eGenes[1:2,], row.names=FALSE)
		cat("......\n")
	}
	if(!is.null(x$cGenes)){
		cat("$cGenes:\n")
		print(x$cGenes[1:2,], row.names=FALSE)
		cat("......\n")
	}
}

######################################################################
# eTarget
######################################################################
#' @title Definition for S3 class \code{eTarget}
#' @description \code{eTarget} has 2 components: evidence and metag.
#' @param evidence a data frame
#' @param metag an 'igraph' object
#' @return an object of S3 class \code{eTarget}
#' @keywords S3 classes
#' @export
#' @examples
#' \dontrun{
#' # Load the library
#' library(Pi)
#' }
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' \dontrun{
#' eTarget(evidence, metag)
#' }
eTarget <- function(evidence, metag){
	## integrity checks
	if(!is(evidence,'data.frame') | !is(metag,'igraph')){
		stop("The S3 class 'eTarget' object failed to pass integrity checks!\n")
	}
	value <- list(evidence=evidence, metag=metag)
	class(value) <- "eTarget"
	return(value)
}
#' @param x an object of class \code{eTarget}
#' @param ... other parameters
#' @rdname eTarget
#' @export
print.eTarget <- function(x, ...) {
	cat(sprintf("An object of S3 class '%s', with %d components:", class(x), length(names(x))), "\n", sep="")
	cat(sprintf("  $metag: an igraph object with %d nodes and %d edges", vcount(x$metag), ecount(x$metag)), "\n", sep="")
	cat(sprintf("  $evidence: a data frame of %d rows X %d columns", dim(x$evidence)[1],dim(x$evidence)[2]), "\n", sep="")
	cat("\n--------------------------------------------------\n")
	cat("$evidence:\n")
	print(x$evidence[1:2,], row.names=FALSE)
	cat("......\n")
}

######################################################################
# dTarget
######################################################################
#' @title Definition for S3 class \code{dTarget}
#' @description \code{dTarget} has 3 components: priority, predictor and metag.
#' @param priority a data frame
#' @param predictor a data frame
#' @param metag an 'igraph' object
#' @return an object of S3 class \code{dTarget}
#' @keywords S3 classes
#' @export
#' @examples
#' \dontrun{
#' # Load the library
#' library(Pi)
#' }
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' \dontrun{
#' dTarget(priority, predictor, metag)
#' }
dTarget <- function(priority, predictor, metag){
	## integrity checks
	if(!is(priority,'data.frame') | !is(predictor,'data.frame') | !is(metag,'igraph')){
		stop("The S3 class 'dTarget' object failed to pass integrity checks!\n")
	}
	value <- list(priority=priority, predictor=predictor, metag=metag)
	class(value) <- "dTarget"
	return(value)
}
#' @param x an object of class \code{dTarget}
#' @param ... other parameters
#' @rdname dTarget
#' @export
print.dTarget <- function(x, ...) {
	cat(sprintf("An object of S3 class '%s', with %d components:", class(x), length(names(x))), "\n", sep="")
	if(!is.null(x$metag)){
		cat(sprintf("  $metag: an igraph object with %d nodes and %d edges", vcount(x$metag), ecount(x$metag)), "\n", sep="")
	}
	cat(sprintf("  $predictor: a data frame of %d rows X %d columns", dim(x$predictor)[1],dim(x$predictor)[2]), "\n", sep="")
	cat(sprintf("  $priority: a data frame of %d rows X %d columns", dim(x$priority)[1],dim(x$priority)[2]), "\n", sep="")
	cat(sprintf("  $list_pNode: a list of %d 'pNode' objects", length(x$list_pNode)), "\n", sep="")
	if(!is.null(x$pPerf)){
		cat(sprintf("  $pPerf: an object of the class 'pPerf'"), "\n", sep="")
	}
	cat("\n--------------------------------------------------\n")
	cat("$priority:\n")
	print(x$priority[1:2,], row.names=FALSE)
	cat("......\n")
	if(!is.null(x$pPerf)){
		cat("$pPerf:\n")
		print(x$pPerf)
		cat("......\n")
	}
}

######################################################################
# sTarget
######################################################################
#' @title Definition for S3 class \code{sTarget}
#' @description \code{sTarget} mush have following components: priority, predictor, performance, importance, evidence.
#' @param priority a data frame
#' @param predictor a data frame
#' @param performance a data frame
#' @param importance a data frame
#' @param evidence an 'eTarget' object
#' @return an object of S3 class \code{sTarget}
#' @keywords S3 classes
#' @export
#' @examples
#' \dontrun{
#' # Load the library
#' library(Pi)
#' }
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' \dontrun{
#' sTarget(priority, predictor, performance, importance, evidence)
#' }
sTarget <- function(priority, predictor, performance, importance, evidence){
	## integrity checks
	if(!is(priority,'data.frame') | !is(evidence,'eTarget')){
		stop("The S3 class 'sTarget' object failed to pass integrity checks!\n")
	}
	value <- list(priority=priority, predictor=predictor, performance=performance, importance=importance, evidence=evidence)
	class(value) <- "sTarget"
	return(value)
}
#' @param x an object of class \code{sTarget}
#' @param ... other parameters
#' @rdname sTarget
#' @export
print.sTarget <- function(x, ...) {
	cat(sprintf("An object of S3 class '%s', with %d components including:", class(x), length(names(x))), "\n", sep="")
	cat(sprintf("  $priority: a data frame of %d rows X %d columns", dim(x$priority)[1],dim(x$priority)[2]), "\n", sep="")
	cat(sprintf("  $predictor: a data frame of %d rows X %d columns", dim(x$predictor)[1],dim(x$predictor)[2]), "\n", sep="")
	cat(sprintf("  $performance: a data frame of %d rows X %d columns", dim(x$performance)[1],dim(x$performance)[2]), "\n", sep="")
	cat(sprintf("  $importance: a data frame of %d rows X %d columns", dim(x$importance)[1],dim(x$importance)[2]), "\n", sep="")
	cat(sprintf("  $evidence: an object of the class 'eTarget'"), "\n", sep="")
	cat(sprintf("  $list_pNode: a list of %d 'pNode' objects", length(x$list_pNode)), "\n", sep="")
	cat("\n--------------------------------------------------\n")
	cat("$priority:\n")
	print(x$priority[1:2,], row.names=FALSE)
	cat("......\n")
	cat("$performance:\n")
	print(x$performance[1:2,], row.names=FALSE)
	cat("......\n")
	cat("$importance:\n")
	print(x$importance[1:2,], row.names=FALSE)
	cat("......\n")
	cat("$evidence:\n")
	print(x$evidence)
	cat("......\n")
}

######################################################################
# cTarget
######################################################################
#' @title Definition for S3 class \code{cTarget}
#' @description \code{cTarget} has 2 components: priority and predictor.
#' @param priority a data frame
#' @param predictor a data frame
#' @return an object of S3 class \code{cTarget}
#' @keywords S3 classes
#' @export
#' @examples
#' \dontrun{
#' # Load the library
#' library(Pi)
#' }
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' \dontrun{
#' cTarget(priority, predictor)
#' }
cTarget <- function(priority, predictor){
	## integrity checks
	if(!is(priority,'data.frame') | !is(predictor,'data.frame')){
		stop("The S3 class 'cTarget' object failed to pass integrity checks!\n")
	}
	value <- list(priority=priority, predictor=predictor)
	class(value) <- "cTarget"
	return(value)
}
#' @param x an object of class \code{cTarget}
#' @param ... other parameters
#' @rdname cTarget
#' @export
print.cTarget <- function(x, ...) {
	cat(sprintf("An object of S3 class '%s', with %d components:", class(x), length(names(x))), "\n", sep="")
	cat(sprintf("  $predictor: a data frame of %d rows X %d columns", dim(x$predictor)[1],dim(x$predictor)[2]), "\n", sep="")
	cat(sprintf("  $priority: a data frame of %d rows X %d columns", dim(x$priority)[1],dim(x$priority)[2]), "\n", sep="")
	cat("\n--------------------------------------------------\n")
	cat("$priority:\n")
	print(x$priority[1:2,], row.names=FALSE)
	cat("......\n")
	cat("$predictor:\n")
	print(x$predictor[1:2,], row.names=FALSE)
	cat("......\n")
}

######################################################################
# sGS
######################################################################
#' @title Definition for S3 class \code{sGS}
#' @description \code{sGS} mush have following components: GSN, GSP, g.
#' @param GSN a vector
#' @param GSP a vector
#' @param g an 'igraph' object
#' @return an object of S3 class \code{sGS}
#' @keywords S3 classes
#' @export
#' @examples
#' \dontrun{
#' # Load the library
#' library(Pi)
#' }
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' \dontrun{
#' sGS(GSN, GSP, g)
#' }
sGS <- function(GSN, GSP, g){
	## integrity checks
	if(!is(GSN,'vector') | !is(g,'igraph')){
		stop("The S3 class 'sGS' object failed to pass integrity checks!\n")
	}
	value <- list(GSN=GSN, GSP=GSP, g=g)
	class(value) <- "sGS"
	return(value)
}
#' @param x an object of class \code{sGS}
#' @param ... other parameters
#' @rdname sGS
#' @export
print.sGS <- function(x, ...) {
	cat(sprintf("An object of S3 class '%s', with %d components including:", class(x), length(names(x))), "\n", sep="")
	cat(sprintf("  $GSN: a vector (%d in total)", length(x$GSN)), "\n", sep="")
	cat(sprintf("  $GSP: a vector (%d in total)", length(x$GSP)), "\n", sep="")
	cat(sprintf("  $g: an igraph object with %d nodes and %d edges", vcount(x$g), ecount(x$g)), "\n", sep="")
	cat("\n--------------------------------------------------\n")
}

######################################################################
# pPerf
######################################################################
#' @title Definition for S3 class \code{pPerf}
#' @description \code{pPerf} mush have following components: PRS, AUROC, Fmax, ROC_perf, PR_perf, Pred_obj.
#' @param PRS a data frame
#' @param AUROC a scalar
#' @param Fmax a scalar
#' @param ROC_perf a ROCR 'performance' object for ROC curve
#' @param PR_perf a ROCR 'performance' object for PR curve
#' @param Pred_obj a ROCR 'prediction' object for other performance measures
#' @return an object of S3 class \code{pPerf}
#' @keywords S3 classes
#' @export
#' @examples
#' \dontrun{
#' # Load the library
#' library(Pi)
#' }
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' \dontrun{
#' pPerf(PRS, AUROC, Fmax, ROC_perf, PR_perf, Pred_obj)
#' }
pPerf <- function(PRS, AUROC, Fmax, ROC_perf, PR_perf, Pred_obj){
	## integrity checks
	if(!is(PRS,'data.frame') | !is(ROC_perf,'performance') | !is(Pred_obj,'prediction')){
		stop("The S3 class 'pPerf' object failed to pass integrity checks!\n")
	}
	value <- list(PRS=PRS, AUROC=AUROC, Fmax=Fmax, ROC_perf=ROC_perf, PR_perf=PR_perf, Pred_obj=Pred_obj)
	class(value) <- "pPerf"
	return(value)
}
#' @param x an object of class \code{pPerf}
#' @param ... other parameters
#' @rdname pPerf
#' @export
print.pPerf <- function(x, ...) {
	cat(sprintf("An object of S3 class '%s', with %d components including:", class(x), length(names(x))), "\n", sep="")
	cat(sprintf("  $PRS: a data frame of %d rows X %d columns", dim(x$PRS)[1],dim(x$PRS)[2]), "\n", sep="")
	cat(sprintf("  $AUROC: %.3f", x$AUROC), "\n", sep="")
	cat(sprintf("  $Fmax: %.3f", x$Fmax), "\n", sep="")
	cat(sprintf("  $ROC_perf: an object of S3 class '%s'", class(x$ROC_perf)), "\n", sep="")
	cat(sprintf("  $PR_perf: an object of S3 class '%s'", class(x$PR_perf)), "\n", sep="")
	cat(sprintf("  $Pred_obj: an object of S3 class '%s'", class(x$Pred_obj)), "\n", sep="")
	cat("\n--------------------------------------------------\n")
}

######################################################################
# eGSEA
######################################################################
#' @title Definition for S3 class \code{eGSEA}
#' @description \code{eGSEA} mush have following components: df_summary, leading, full, cross.
#' @param df_summary a data frame
#' @param leading a list
#' @param full a list
#' @param cross a matrix
#' @return an object of S3 class \code{eGSEA}
#' @keywords S3 classes
#' @export
#' @examples
#' \dontrun{
#' # Load the library
#' library(Pi)
#' }
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' \dontrun{
#' eGSEA(df_summary, leading, full, cross)
#' }
eGSEA <- function(df_summary, leading, full, cross){
	## integrity checks
	if(!is(df_summary,'data.frame')){
		stop("The S3 class 'eGSEA' object failed to pass integrity checks!\n")
	}
	value <- list(df_summary=df_summary, leading=leading, full=full, cross=cross)
	class(value) <- "eGSEA"
	return(value)
}
#' @param x an object of class \code{eGSEA}
#' @param ... other parameters
#' @rdname eGSEA
#' @export
print.eGSEA <- function(x, ...) {
	cat(sprintf("An object of S3 class '%s', with %d components including:", class(x), length(names(x))), "\n", sep="")
	cat(sprintf("  $df_summary: a data frame of %d rows X %d columns", dim(x$df_summary)[1],dim(x$df_summary)[2]), "\n", sep="")
	cat(sprintf("  $leading: a list (%d in total)", length(x$leading)), "\n", sep="")
	cat(sprintf("  $full: a list (%d in total)", length(x$full)), "\n", sep="")
	cat(sprintf("  $cross: a matrix of %d X %d", dim(x$cross)[1], dim(x$cross)[2]), "\n", sep="")
	cat("\n--------------------------------------------------\n")
	cat("$df_summary:\n")
	print(x$df_summary[1:min(2,nrow(x$df_summary)),], row.names=FALSE)
	cat("......\n")
}

######################################################################
# eTerm
######################################################################
#' @title Definition for S3 class \code{eTerm}
#' @description \code{eTerm} mush have following components: term_info, annotation, g, data, background, overlap, fc, zscore, pvalue, adjp, cross.
#' @param term_info a data frame
#' @param annotation a list
#' @param g an 'igraph' object
#' @param data a vector
#' @param background a vector
#' @param overlap a vector
#' @param fc a vector
#' @param zscore a vector
#' @param pvalue a vector
#' @param adjp a vector
#' @param cross a matrix
#' @return an object of S3 class \code{eTerm}
#' @keywords S3 classes
#' @export
#' @examples
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' \dontrun{
#' eTerm(term_info, annotation, g, data, background, overlap, fc, zscore, pvalue, adjp, cross)
#' }
eTerm <- function(term_info, annotation, g, data, background, overlap, fc, zscore, pvalue, adjp, cross){
	## integrity checks
	if(!is(term_info,'data.frame') | !is(g,'igraph')){
		stop("The S3 class 'eTerm' object failed to pass integrity checks!\n")
	}
	value <- list(term_info=term_info, annotation=annotation, g=g, data=data, background=background, overlap=overlap, fc=fc, zscore=zscore, pvalue=pvalue, adjp=adjp, cross=cross)
	class(value) <- "eTerm"
	return(value)
}
#' @param x an object of class \code{eTerm}
#' @param ... other parameters
#' @rdname eTerm
#' @export
#' @method print eTerm
print.eTerm <- function(x, ...) {
	cat(sprintf("An object of S3 class '%s', with %d components including:", class(x), length(names(x))), "\n", sep="")
	cat(sprintf("  $term_info: a data frame of %d rows X %d columns", dim(x$term_info)[1],dim(x$term_info)[2]), "\n", sep="")
	cat(sprintf("  $data: a vector (%d in total)", length(x$data)), "\n", sep="")
	cat(sprintf("  $background: a vector (%d in total)", length(x$background)), "\n", sep="")
	cat(sprintf("  $adjp: a vector (%d in total)", length(x$adjp)), "\n", sep="")
	cat(sprintf("  $cross: a matrix of %d X %d", dim(x$cross)[1], dim(x$cross)[2]), "\n", sep="")
	cat(sprintf("  $g: an 'igraph' object"), "\n", sep="")
	cat(sprintf("  $g$ontology: '%s'", x$g$ontology), "\n", sep="")
	cat("\n--------------------------------------------------\n")
	cat("xEnrichViewer(eTerm):\n")
	print(xEnrichViewer(x), row.names=TRUE)
	cat("......\n")
}

######################################################################
# ls_eTerm
######################################################################
#' @title Definition for S3 class \code{ls_eTerm}
#' @description \code{ls_eTerm} has 3 components: df, mat and gp.
#' @param df a data frame
#' @param mat a matrix
#' @param gp a ggplot object
#' @return an object of S3 class \code{ls_eTerm}
#' @keywords S3 classes
#' @export
#' @examples
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' \dontrun{
#' ls_eTerm(df, mat, gp)
#' }
ls_eTerm <- function(df, mat, gp){
	## integrity checks
	if(!is(df,'data.frame') | !is(mat,'matrix') | is(gp,'ggplot')){
		stop("The S3 class 'ls_eTerm' object failed to pass integrity checks!\n")
	}
	value <- list(df=df, mat=mat, gp=gp)
	class(value) <- "ls_eTerm"
	return(value)
}
#' @param x an object of class \code{ls_eTerm}
#' @param ... other parameters
#' @rdname ls_eTerm
#' @export
#' @method print ls_eTerm
print.ls_eTerm <- function(x, ...) {
	cat(sprintf("An object of S3 class '%s', with %d components:", class(x), length(names(x))), "\n", sep="")
	cat(sprintf("  $df: a data frame of %d rows X %d columns", dim(x$df)[1],dim(x$df)[2]), "\n", sep="")
	cat(sprintf("  $mat: a data matrix of %d rows X %d columns", dim(x$mat)[1],dim(x$mat)[2]), "\n", sep="")
	cat(sprintf("  $gp: a ggplot object"), "\n", sep="")
	cat("\n--------------------------------------------------\n")
	cat("$df:\n")
	print(x$df[1:min(2,nrow(x$df)),1:13], row.names=FALSE)
	cat("......\n")
}

######################################################################
# aOnto
######################################################################
#' @title Definition for S3 class \code{aOnto}
#' @description \code{aOnto} has 2 components: g, anno.
#' @param g an igraph object
#' @param anno a list
#' @return an object of S3 class \code{aOnto}
#' @keywords S3 classes
#' @export
#' @examples
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' \dontrun{
#' aOnto(g, anno)
#' }
aOnto <- function(g, anno){
	## integrity checks
	if(!is(g,'igraph') | !is(anno,'list')){
		stop("The S3 class 'aOnto' object failed to pass integrity checks!\n")
	}
	value <- list(g=g, anno=anno)
	class(value) <- "aOnto"
	return(value)
}
#' @param x an object of class \code{aOnto}
#' @param ... other parameters
#' @rdname aOnto
#' @export
#' @method print aOnto
print.aOnto <- function(x, ...) {
	cat(sprintf("An object of S3 class '%s', with %d components:", class(x), length(names(x))), "\n", sep="")
	cat(sprintf("  $g: an igraph object or NULL"), "\n", sep="")
	cat(sprintf("  $anno: a list with %d length or NULL", length(x$anno)), "\n", sep="")
	
	cat("\n--------------------------------------------------\n")
	cat("$g:\n")
	print(x$g)
	cat("......\n")
}

######################################################################
# GS
######################################################################
#' @title Definition for S3 class \code{GS}
#' @description \code{GS} has 2 components: set_info, gs.
#' @param set_info a data frame
#' @param gs a list
#' @return an object of S3 class \code{GS}
#' @keywords S3 classes
#' @export
#' @examples
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' \dontrun{
#' GS(set_info, gs)
#' }
GS <- function(set_info, gs){
	## integrity checks
	if(!is(set_info,'data.frame') | !is(gs,'list')){
		stop("The S3 class 'GS' object failed to pass integrity checks!\n")
	}
	value <- list(set_info=set_info, gs=gs)
	class(value) <- "GS"
	return(value)
}
#' @param x an object of class \code{GS}
#' @param ... other parameters
#' @rdname GS
#' @export
#' @method print GS
print.GS <- function(x, ...) {
	cat(sprintf("An object of S3 class '%s', with %d components:", class(x), length(names(x))), "\n", sep="")
	cat(sprintf("  $set_info: a data frame of %d rows X %d columns", nrow(x$set_info), ncol(x$set_info)), "\n", sep="")
	cat(sprintf("  $gs: a list with %d length", length(x$gs)), "\n", sep="")
	
	cat("\n--------------------------------------------------\n")
	cat("$set_info:\n")
	print(x$set_info[1:2,], row.names=FALSE)
	cat("......\n")
}

######################################################################
# EG
######################################################################
#' @title Definition for S3 class \code{EG}
#' @description \code{EG} has 1 component: gene_info.
#' @param gene_info a data frame
#' @return an object of S3 class \code{EG}
#' @keywords S3 classes
#' @export
#' @examples
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' \dontrun{
#' EG(gene_info)
#' }
EG <- function(gene_info){
	## integrity checks
	if(!is(gene_info,'data.frame')){
		stop("The S3 class 'EG' object failed to pass integrity checks!\n")
	}
	value <- list(gene_info=gene_info)
	class(value) <- "EG"
	return(value)
}
#' @param x an object of class \code{EG}
#' @param ... other parameters
#' @rdname EG
#' @export
#' @method print EG
print.EG <- function(x, ...) {
	cat(sprintf("An object of S3 class '%s', with %d components:", class(x), length(names(x))), "\n", sep="")
	cat(sprintf("  $gene_info: a data frame of %d rows X %d columns", nrow(x$set_info), ncol(x$set_info)), "\n", sep="")
	
	cat("\n--------------------------------------------------\n")
	cat("$gene_info:\n")
	print(x$gene_info[1:2,], row.names=FALSE)
	cat("......\n")
}

######################################################################
# iSubg
######################################################################
#' @title Definition for S3 class \code{iSubg}
#' @description \code{iSubg} has 2 components: g, ls_subg.
#' @param g an igraph object
#' @param ls_subg a list of igraph objects
#' @return an object of S3 class \code{iSubg}
#' @keywords S3 classes
#' @export
#' @examples
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' \dontrun{
#' iSubg(g, ls_subg)
#' }
iSubg <- function(g, ls_subg){
	## integrity checks
	if(!is(g,'igraph') | !is(ls_subg,'list')){
		stop("The S3 class 'iSubg' object failed to pass integrity checks!\n")
	}
	value <- list(g=g, ls_subg=ls_subg)
	class(value) <- "iSubg"
	return(value)
}
#' @param x an object of class \code{iSubg}
#' @param ... other parameters
#' @rdname iSubg
#' @export
#' @method print iSubg
print.iSubg <- function(x, ...) {
	cat(sprintf("An object of S3 class '%s', with %d components:", class(x), length(names(x))), "\n", sep="")
	cat(sprintf("  $g: an igraph object or NULL"), "\n", sep="")
	cat(sprintf("  $ls_subg: a list with %d length or NULL", length(x$ls_subg)), "\n", sep="")
	
	cat("\n--------------------------------------------------\n")
	cat("$g:\n")
	print(x$g)
	cat("......\n")
}

