#' Function to prioritise pathways based on enrichment analysis of top prioritised genes
#'
#' \code{xPierPathways} is supposed to prioritise pathways given prioritised genes and the ontology in query. It returns an object of class "eTerm". It is done via enrichment analysis. 
#'
#' @param pNode an object of class "pNode" (or "sTarget" or "dTarget")
#' @param priority.top the number of the top targets used for enrichment analysis. By default, it sets to 100
#' @param background a background vector. It contains a list of Gene Symbols as the test background. If NULL, by default all annotatable are used as background
#' @param ontology the ontology supported currently. By default, it is 'NA' to disable this option. Pre-built ontology and annotation data are detailed in \code{\link{xDefineOntology}}.
#' @param size.range the minimum and maximum size of members of each term in consideration. By default, it sets to a minimum of 10 but no more than 2000
#' @param min.overlap the minimum number of overlaps. Only those terms with members that overlap with input data at least min.overlap (3 by default) will be processed
#' @param which.distance which terms with the distance away from the ontology root (if any) is used to restrict terms in consideration. By default, it sets to 'NULL' to consider all distances
#' @param test the statistic test used. It can be "fisher" for using fisher's exact test, "hypergeo" for using hypergeometric test, or "binomial" for using binomial test. Fisher's exact test is to test the independence between gene group (genes belonging to a group or not) and gene annotation (genes annotated by a term or not), and thus compare sampling to the left part of background (after sampling without replacement). Hypergeometric test is to sample at random (without replacement) from the background containing annotated and non-annotated genes, and thus compare sampling to background. Unlike hypergeometric test, binomial test is to sample at random (with replacement) from the background with the constant probability. In terms of the ease of finding the significance, they are in order: hypergeometric test > binomial test > fisher's exact test. In other words, in terms of the calculated p-value, hypergeometric test < binomial test < fisher's exact test
#' @param background.annotatable.only logical to indicate whether the background is further restricted to the annotatable. By default, it is NULL: if ontology.algorithm is not 'none', it is always TRUE; otherwise, it depends on the background (if not provided, it will be TRUE; otherwise FALSE). Surely, it can be explicitly stated
#' @param p.tail the tail used to calculate p-values. It can be either "two-tails" for the significance based on two-tails (ie both over- and under-overrepresentation)  or "one-tail" (by default) for the significance based on one tail (ie only over-representation)
#' @param p.adjust.method the method used to adjust p-values. It can be one of "BH", "BY", "bonferroni", "holm", "hochberg" and "hommel". The first two methods "BH" (widely used) and "BY" control the false discovery rate (FDR: the expected proportion of false discoveries amongst the rejected hypotheses); the last four methods "bonferroni", "holm", "hochberg" and "hommel" are designed to give strong control of the family-wise error rate (FWER). Notes: FDR is a less stringent condition than FWER
#' @param ontology.algorithm the algorithm used to account for the hierarchy of the ontology. It can be one of "none", "pc", "elim" and "lea". For details, please see 'Note' below
#' @param elim.pvalue the parameter only used when "ontology.algorithm" is "elim". It is used to control how to declare a signficantly enriched term (and subsequently all genes in this term are eliminated from all its ancestors)
#' @param lea.depth the parameter only used when "ontology.algorithm" is "lea". It is used to control how many maximum depth is used to consider the children of a term (and subsequently all genes in these children term are eliminated from the use for the recalculation of the signifance at this term)
#' @param path.mode the mode of paths induced by vertices/nodes with input annotation data. It can be "all_paths" for all possible paths to the root, "shortest_paths" for only one path to the root (for each node in query), "all_shortest_paths" for all shortest paths to the root (i.e. for each node, find all shortest paths with the equal lengths)
#' @param true.path.rule logical to indicate whether the true-path rule should be applied to propagate annotations. By default, it sets to false
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to false for no display
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @param guid a valid (5-character) Global Unique IDentifier for an OSF project. See \code{\link{xRDataLoader}} for details
#' @return 
#' an object of class "eTerm", a list with following components:
#' \itemize{
#'  \item{\code{term_info}: a matrix of nTerm X 4 containing snp/gene set information, where nTerm is the number of terms, and the 4 columns are "id" (i.e. "Term ID"), "name" (i.e. "Term Name"), "namespace" and "distance"}
#'  \item{\code{annotation}: a list of terms containing annotations, each term storing its annotations. Always, terms are identified by "id"}
#'  \item{\code{data}: a vector containing input data in consideration. It is not always the same as the input data as only those mappable are retained}
#'  \item{\code{background}: a vector containing the background data. It is not always the same as the input data as only those mappable are retained}
#'  \item{\code{overlap}: a list of overlapped snp/gene sets, each storing snps overlapped between a snp/gene set and the given input data (i.e. the snps of interest). Always, gene sets are identified by "id"}
#'  \item{\code{zscore}: a vector containing z-scores}
#'  \item{\code{pvalue}: a vector containing p-values}
#'  \item{\code{adjp}: a vector containing adjusted p-values. It is the p value but after being adjusted for multiple comparisons}
#'  \item{\code{call}: the call that produced this result}
#' }
#' @note The interpretation of the algorithms used to account for the hierarchy of the ontology is:
#' \itemize{
#' \item{"none": does not consider the ontology hierarchy at all.}
#' \item{"lea": computers the significance of a term in terms of the significance of its children at the maximum depth (e.g. 2). Precisely, once snps are already annotated to any children terms with a more signficance than itself, then all these snps are eliminated from the use for the recalculation of the signifance at that term. The final p-values takes the maximum of the original p-value and the recalculated p-value.}
#' \item{"elim": computers the significance of a term in terms of the significance of its all children. Precisely, once snps are already annotated to a signficantly enriched term under the cutoff of e.g. pvalue<1e-2, all these snps are eliminated from the ancestors of that term).}
#' \item{"pc": requires the significance of a term not only using the whole snps as background but also using snps annotated to all its direct parents/ancestors as background. The final p-value takes the maximum of both p-values in these two calculations.}
#' \item{"Notes": the order of the number of significant terms is: "none" > "lea" > "elim" > "pc".}
#' }
#' @export
#' @seealso \code{\link{xEnricherGenes}}
#' @include xPierPathways.r
#' @examples
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' \dontrun{
#' # a) provide the SNPs with the significance info
#' ## get lead SNPs reported in AS GWAS and their significance info (p-values)
#' #data.file <- "http://galahad.well.ox.ac.uk/bigdata/AS.txt"
#' #AS <- read.delim(data.file, header=TRUE, stringsAsFactors=FALSE)
#' ImmunoBase <- xRDataLoader(RData.customised='ImmunoBase')
#' gr <- ImmunoBase$AS$variants
#' AS <- as.data.frame(GenomicRanges::mcols(gr)[, c('Variant','Pvalue')])
#'
#' # b) perform priority analysis
#' pNode <- xPierSNPs(data=AS, include.eQTL="JKng_mono", include.HiC='Monocytes', network="PCommonsUN_medium", restart=0.7, RData.location=RData.location)
#' 
#' # c) derive pathway-level priority
#' eTerm <- xPierPathways(pNode=pNode, priority.top=100, ontology="MsigdbC2CP", RData.location=RData.location)
#'
#' # d) view enrichment results for the top significant terms
#' xEnrichViewer(eTerm)
#'
#' # e) save enrichment results to the file called 'Pathways_priority.txt'
#' res <- xEnrichViewer(eTerm, top_num=length(eTerm$adjp), sortBy="adjp", details=TRUE)
#' output <- data.frame(term=rownames(res), res)
#' utils::write.table(output, file="Pathways_priority.txt", sep="\t", row.names=FALSE)
#' }

xPierPathways <- function(pNode, priority.top=100, background=NULL, ontology=NA, size.range=c(10,2000), min.overlap=3, which.distance=NULL, test=c("hypergeo","fisher","binomial"), background.annotatable.only=NULL, p.tail=c("one-tail","two-tails"), p.adjust.method=c("BH", "BY", "bonferroni", "holm", "hochberg", "hommel"), ontology.algorithm=c("none","pc","elim","lea"), elim.pvalue=1e-2, lea.depth=2, path.mode=c("all_paths","shortest_paths","all_shortest_paths"), true.path.rule=FALSE, verbose=TRUE, RData.location="http://galahad.well.ox.ac.uk/bigdata", guid=NULL)
{
    startT <- Sys.time()
    message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=TRUE)
    message("", appendLF=TRUE)
    ####################################################################################
    
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    #ontology <- match.arg(ontology)
    test <- match.arg(test)
    p.adjust.method <- match.arg(p.adjust.method)
    ontology.algorithm <- match.arg(ontology.algorithm)
    path.mode <- match.arg(path.mode)
    p.tail <- match.arg(p.tail)
    
    if(is(pNode,"pNode")){
        df_priority <- pNode$priority[, c("seed","weight","priority")]
    }else if(is(pNode,"sTarget") | is(pNode,"dTarget")){
    	df_priority <- pNode$priority[, c("name","rank","rating")]
    	df_priority$priority <- df_priority$rating
    }else{
    	stop("The function must apply to a 'pNode' or 'sTarget' or 'dTarget' object.\n")
    }
    
	## priority top
	priority.top <- as.integer(priority.top)
    if ( priority.top > nrow(df_priority) ){
        priority.top <- nrow(df_priority)
    }
    
    data <- rownames(df_priority)[1:priority.top]
    if(is.null(background)){
    	background <- rownames(df_priority)
    }
    
    #############################################################################################
    
    if(verbose){
        now <- Sys.time()
        message(sprintf("\n#######################################################"))
        message(sprintf("'xEnricherGenes' is being called (%s):", as.character(now)), appendLF=TRUE)
        message(sprintf("#######################################################"))
    }
    
	eTerm <- xEnricherGenes(data=data, background=background, ontology=ontology, size.range=size.range, min.overlap=min.overlap, which.distance=which.distance, test=test, background.annotatable.only=background.annotatable.only, p.tail=p.tail, p.adjust.method=p.adjust.method, ontology.algorithm=ontology.algorithm, elim.pvalue=elim.pvalue, lea.depth=lea.depth, path.mode=path.mode, true.path.rule=true.path.rule, verbose=verbose, RData.location=RData.location, guid=guid)
	
	if(verbose){
        now <- Sys.time()
        message(sprintf("#######################################################"))
        message(sprintf("'xEnricherGenes' has been finished (%s)!", as.character(now)), appendLF=TRUE)
        message(sprintf("#######################################################\n"))
    }
    
    ####################################################################################
    endT <- Sys.time()
    message(paste(c("\nEnd at ",as.character(endT)), collapse=""), appendLF=TRUE)
    
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=TRUE)
    
    invisible(eTerm)
}
