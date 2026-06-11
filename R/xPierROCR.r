#' Function to assess the dTarget performance via ROC and Precision-Recall (PR) analysis
#'
#' \code{xPierROCR} is supposed to assess the dTarget performance via Receiver Operating Characteristic (ROC) and Precision-Recall (PR) analysis. It requires three inputs: 1) Gold Standard Positive (GSP) targets; 2) Gold Standard Negative (GSN) targets; 3) dTarget containing predicted targets and predictive scores.
#'
#' @param dTarget a data frame containing dTargets along with predictive scores. It has two columns: 1st column for target, 2nd column for predictive scores (the higher the better). Alternatively, it can be an object of class "pNode" (or "sTarget" or "dTarget") from which a data frame is extracted
#' @param GSP a vector containing Gold Standard Positives (GSP)
#' @param GSN a vector containing Gold Standard Negatives (GSN)
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to TRUE for display
#' @return 
#' an object of the class "dTarget", a list with following components:
#' \itemize{
#'  \item{\code{priority}: a data frame of nGene X 7 containing gene priority (aggregated) information, where nGene is the number of genes, and the 7 columns are "GS" (either 'GSP', or 'GSN', or 'NEW'), "name" (gene names), "rank" (ranks of the priority scores), "pvalue" (the aggregated p-value, converted from empirical cumulative distribution of the probability of being GSP), "fdr" (fdr adjusted from the aggregated p-value), "priority" (-log10(pvalue) but rescaled into the 5-star ratings), "description" (gene description) and seed info including "Overall" for the number of different types of seeds, followed by details on individual type of seeds (that is, "OMIM", "Phenotype", "Function", "nearbyGenes", "eQTL", "HiC")}
#'  \item{\code{predictor}: a data frame containing predictor matrix, with each column/predictor for either priority score, or priorty rank or priority p-value}
#'  \item{\code{metag}: an "igraph" object}
#'  \item{\code{pPerf}: a "pPerf" object, with components "PRS", "AUROC", "Fmax", "ROC_perf", "PR_perf", "Pred_obj"}
#' }
#' @note
#' AUC: the area under ROC
#' F-measure: the maximum of a harmonic mean between precision and recall along PR curve
#' @export
#' @seealso \code{\link{xPierROCR}}
#' @include xPierROCR.r
#' @examples
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' \dontrun{
#' dTarget <- xPierROCR(dTarget, GSP, GSN)
#' gp <- xPredictCompare(dTarget$pPerf)
#' }

xPierROCR <- function(dTarget, GSP, GSN, verbose=TRUE)
{

    if (!is(dTarget,"dTarget")){
    	stop("The function must apply to a 'dTarget' object.\n")
    }

	## pre-process GSP and GSN
	gsp <- unique(GSP)
	gsn <- unique(GSN)
	gsn <- setdiff(gsn, gsp)
	gs_names <- union(gsp, gsn)
	gs_targets <- rep(0, length(gs_names))
	names(gs_targets) <- gs_names
	gs_targets[gsp] <- 1
	### output_gs
	output_gs <- rep('NEW', nrow(dTarget$priority))
	names(output_gs) <- rownames(dTarget$priority)
	ind <- match(rownames(dTarget$priority), names(gs_targets))
	output_gs[!is.na(ind)] <- gs_targets[ind[!is.na(ind)]]
	output_gs[output_gs=='0'] <- 'GSN'
	output_gs[output_gs=='1'] <- 'GSP'
	### apend GS to dTarget$priority
	dTarget$priority <- cbind(GS=output_gs, dTarget$priority)
	
	## add a new component: pPerf
	pPerf <- xPredictROCR(prediction=dTarget, GSP=GSP, GSN=GSN, rescale=FALSE, verbose=verbose)
	
	dTarget$pPerf <- pPerf
	
	invisible(dTarget)
}
