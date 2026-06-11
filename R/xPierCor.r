#' Function to calculate correlation between prioritised genes and user-defined external data
#'
#' \code{xPierCor} is supposed to calculate correlation between prioritised genes and user-defined external data. 
#
#' @param pNode an object of class "pNode" (or "sTarget" or "dTarget"). Alternatively, it can be a data frame with two columns ('name' and 'priority')
#' @param list_vec a named vector containing numeric values for genes (gene symbols). Alternatively it can be a list of named vectors
#' @param method the method used to calcualte correlation. It can be 'pearson' for Pearson's correlation or 'spearman' for Spearman rank correlation
#' @param pvalue.type the type of the p-value calcualted. It can be 'nominal' for nominal p-value or 'empirical' for empirical p-value
#' @param seed an integer specifying the seed
#' @param nperm the number of random permutations
#' @param p.adjust.method the method used to adjust p-values. It can be one of "BH", "BY", "bonferroni", "holm", "hochberg" and "hommel". The first two methods "BH" (widely used) and "BY" control the false discovery rate (FDR: the expected proportion of false discoveries amongst the rejected hypotheses); the last four methods "bonferroni", "holm", "hochberg" and "hommel" are designed to give strong control of the family-wise error rate (FWER). Notes: FDR is a less stringent condition than FWER
#' @param plot logical to indicate whether scatter plot is drawn
#' @return 
#' a list with two componets:
#' \itemize{
#'  \item{\code{df_summary}: a data frame of n x 4, where n is the number of named vectors, and the 4 columns are "name", "cor" (i.e. "correlation"), "pval" (i.e. p-value), "fdr"}
#'  \item{\code{ls_gp}: NULL if the plot is not drawn; otherwise, a list of 'ggplot' objects}
#' }
#' @note none
#' @export
#' @seealso \code{\link{xCorrelation}}
#' @include xPierCor.r
#' @examples
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' \dontrun{
#' # a) provide the seed nodes/genes with the weight info
#' ## load ImmunoBase
#' ImmunoBase <- xRDataLoader(RData.customised='ImmunoBase', RData.location=RData.location)
#' ## get genes within 500kb away from AS GWAS lead SNPs
#' seeds.genes <- ImmunoBase$AS$genes_variants
#' ## seeds weighted according to distance away from lead SNPs
#' data <- 1- seeds.genes/500000
#'
#' # b) perform priority analysis
#' pNode <- xPierGenes(data=data, network="PCommonsDN_medium",restart=0.7, RData.location=RData.location)
#' 
#' # c) do correlation
#' data <- pNode$priority$priority[1:100]
#' name(data) <- pNode$priority$name[1:100]
#' ls_res <- xPierCor(pNode, data, method="pearson", pvalue.type="empirical", nperm=2000, plot=TRUE)
#' }

xPierCor <- function(pNode, list_vec, method=c("pearson","spearman"), pvalue.type=c("nominal","empirical"), seed=825, nperm=2000, p.adjust.method=c("BH","BY","bonferroni","holm","hochberg","hommel"), plot=FALSE)
{
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    method <- match.arg(method)
    p.type <- match.arg(p.type)
    p.adjust.method <- match.arg(p.adjust.method)
    
    if(is(pNode,"pNode")){
        df_priority <- pNode$priority[, c("name","priority")]
    }else if(is(pNode,"sTarget") | is(pNode,"dTarget")){
    	df_priority <- pNode$priority[, c("name","priority")]
    }else if(is(pNode,"data.frame")){
    	df_priority <- pNode[,c(1:2)]
    	colnames(df_priority) <- c("name","priority")
    }else{
    	stop("The function must apply to a 'pNode' or 'sTarget' or 'dTarget' object.\n")
    }
    
    ls_res <- xCorrelation(df=df_priority, list_vec=list_vec, method=method, p.type=p.type, seed=seed, nperm=nperm, p.adjust.method=p.adjust.method, plot=plot)
    
    if(plot){
    	ls_res$ls_gp <- lapply(ls_res$ls_gp, function(gp_curve){
    		gp_curve <- gp_curve + xlab("Pi's 5-star rating")
			
			if(0){
			df <- gp_curve$data
			priority <- name <- NULL
			gp_curve <- gp_curve + scale_x_continuous(limits=c(0,ceiling(max(df$priority)*10)/10)) 
			gp_curve <- gp_curve + scale_y_reverse(limits=c(ceiling(max(df$data)*10)/10, floor(min(df$data)*10)/10))
			gp_curve <- gp_curve + ggrepel::geom_text_repel(aes(label=name), size=2, fontface='bold', box.padding=unit(0.2,"lines"), point.padding=unit(0.2,"lines"), segment.color='grey50', segment.alpha=0.5, arrow=arrow(length=unit(0.01,'npc')), force=0.1)
			}
			
			gp_curve
    	})
    }
    
    invisible(ls_res)
}
