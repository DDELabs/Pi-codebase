#' Function to prioritise genes given a list of genomic regions
#'
#' \code{xPierGRs} is supposed to prioritise genes given a list of genomic regions with or without the significance level. To prioritise genes, it first defines and scores genes crosslinking to an input list of genomic regions (GR). With seed genes and their scores, it then uses Random Walk with Restart (RWR) to calculate the affinity score of all nodes in the input graph to the seed genes. The priority score is the affinity score. Parallel computing is also supported for Linux-like or Windows operating systems. It returns an object of class "pNode".
#'
#' @param data a named input vector containing the significance level for genomic regions (GR). For this named vector, the element names are GR, in the format of 'chrN:start-end', where N is either 1-22 or X, start (or end) is genomic positional number; for example, 'chr1:13-20', the element values for the significance level (measured as p-value or fdr). Alternatively, it can be a matrix or data frame with two columns: 1st column for GR, 2nd column for the significance level. Also supported is the input with GR only (without the significance level)
#' @param significance.threshold the given significance threshold. By default, it is set to NULL, meaning there is no constraint on the significance level when transforming the significance level of GR into scores. If given, those GR below this are considered significant and thus scored positively. Instead, those above this are considered insignificant and thus receive no score
#' @param score.cap the maximum score being capped. By default, it is set to NULL, meaning that no capping is applied
#' @param build.conversion the conversion from one genome build to another. The conversions supported are "hg38.to.hg19" and "hg18.to.hg19". By default it is NA (no need to do so)
#' @param crosslink the built-in crosslink info with a score quantifying the link of a GR to a gene. See \code{\link{xGR2xGenes}} for details
#' @param crosslink.customised the crosslink info with a score quantifying the link of a GR to a gene. A user-input matrix or data frame with 4 columns: 1st column for genomic regions (formatted as "chr:start-end", genome build 19), 2nd column for Genes, 3rd for crosslink score (crosslinking a genomic region to a gene, such as -log10 significance level), and 4th for contexts (optional; if nor provided, it will be added as 'C'). Alternatively, it can be a file containing these 4 columns. Required, otherwise it will return NULL
#' @param cdf.function a character specifying how to transform the input crosslink score. It can be one of 'original' (no such transformation), and 'empirical'  for looking at empirical Cumulative Distribution Function (cdf; as such it is converted into pvalue-like values [0,1])
#' @param scoring.scheme the method used to calculate seed gene scores under a set of GR (also over Contexts if many). It can be one of "sum" for adding up, "max" for the maximum, and "sequential" for the sequential weighting. The sequential weighting is done via: \eqn{\sum_{i=1}{\frac{R_{i}}{i}}}, where \eqn{R_{i}} is the \eqn{i^{th}} rank (in a descreasing order)
#' @param nearby.distance.max the maximum distance between genes and GR. Only those genes no far way from this distance will be considered as seed genes. This parameter will influence the distance-component weights calculated for nearby GR per gene
#' @param nearby.decay.kernel a character specifying a decay kernel function. It can be one of 'slow' for slow decay, 'linear' for linear decay, and 'rapid' for rapid decay. If no distance weight is used, please select 'constant'
#' @param nearby.decay.exponent a numeric specifying a decay exponent. By default, it sets to 2
#' @param network the built-in network. Currently two sources of network information are supported: the STRING database (version 10) and the Pathway Commons database (version 7). STRING is a meta-integration of undirect interactions from the functional aspect, while Pathways Commons mainly contains both undirect and direct interactions from the physical/pathway aspect. Both have scores to control the confidence of interactions. Therefore, the user can choose the different quality of the interactions. In STRING, "STRING_highest" indicates interactions with highest confidence (confidence scores>=900), "STRING_high" for interactions with high confidence (confidence scores>=700), "STRING_medium" for interactions with medium confidence (confidence scores>=400), and "STRING_low" for interactions with low confidence (confidence scores>=150). For undirect/physical interactions from Pathways Commons, "PCommonsUN_high" indicates undirect interactions with high confidence (supported with the PubMed references plus at least 2 different sources), "PCommonsUN_medium" for undirect interactions with medium confidence (supported with the PubMed references). For direct (pathway-merged) interactions from Pathways Commons, "PCommonsDN_high" indicates direct interactions with high confidence (supported with the PubMed references plus at least 2 different sources), and "PCommonsUN_medium" for direct interactions with medium confidence (supported with the PubMed references). In addition to pooled version of pathways from all data sources, the user can also choose the pathway-merged network from individual sources, that is, "PCommonsDN_Reactome" for those from Reactome, "PCommonsDN_KEGG" for those from KEGG, "PCommonsDN_HumanCyc" for those from HumanCyc, "PCommonsDN_PID" for those froom PID, "PCommonsDN_PANTHER" for those from PANTHER, "PCommonsDN_ReconX" for those from ReconX, "PCommonsDN_TRANSFAC" for those from TRANSFAC, "PCommonsDN_PhosphoSite" for those from PhosphoSite, and "PCommonsDN_CTD" for those from CTD. For direct (pathway-merged) interactions sourced from KEGG, it can be 'KEGG' for all, 'KEGG_metabolism' for pathways grouped into 'Metabolism', 'KEGG_genetic' for 'Genetic Information Processing' pathways, 'KEGG_environmental' for 'Environmental Information Processing' pathways, 'KEGG_cellular' for 'Cellular Processes' pathways, 'KEGG_organismal' for 'Organismal Systems' pathways, and 'KEGG_disease' for 'Human Diseases' pathways. 'REACTOME' for protein-protein interactions derived from Reactome pathways
#' @param STRING.only the further restriction of STRING by interaction type. If NA, no such restriction. Otherwide, it can be one or more of "neighborhood_score","fusion_score","cooccurence_score","coexpression_score","experimental_score","database_score","textmining_score". Useful options are c("experimental_score","database_score"): only experimental data (extracted from BIND, DIP, GRID, HPRD, IntAct, MINT, and PID) and curated data (extracted from Biocarta, BioCyc, GO, KEGG, and Reactome) are used
#' @param weighted logical to indicate whether edge weights should be considered. By default, it sets to false. If true, it only works for the network from the STRING database
#' @param network.customised an object of class "igraph". By default, it is NULL. It is designed to allow the user analysing their customised network data that are not listed in the above argument 'network'. This customisation (if provided) has the high priority over built-in network. If the user provides the "igraph" object with the "weight" edge attribute, RWR will assume to walk on the weighted network
#' @param seeds.inclusive logical to indicate whether non-network seed genes are included for prioritisation. If TRUE (by default), these genes will be added to the netowrk
#' @param normalise the way to normalise the adjacency matrix of the input graph. It can be 'laplacian' for laplacian normalisation, 'row' for row-wise normalisation, 'column' for column-wise normalisation, or 'none'
#' @param restart the restart probability used for Random Walk with Restart (RWR). The restart probability takes the value from 0 to 1, controlling the range from the starting nodes/seeds that the walker will explore. The higher the value, the more likely the walker is to visit the nodes centered on the starting nodes. At the extreme when the restart probability is zero, the walker moves freely to the neighbors at each step without restarting from seeds, i.e., following a random walk (RW)
#' @param normalise.affinity.matrix the way to normalise the output affinity matrix. It can be 'none' for no normalisation, 'quantile' for quantile normalisation to ensure that columns (if multiple) of the output affinity matrix have the same quantiles
#' @param parallel logical to indicate whether parallel computation with multicores is used. By default, it sets to true, but not necessarily does so. Partly because parallel backends available will be system-specific (now only Linux or Mac OS). Also, it will depend on whether these two packages "foreach" and "doMC" have been installed
#' @param multicores an integer to specify how many cores will be registered as the multicore parallel backend to the 'foreach' package. If NULL, it will use a half of cores available in a user's computer. This option only works when parallel computation is enabled
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @param guid a valid (5-character) Global Unique IDentifier for an OSF project. See \code{\link{xRDataLoader}} for details
#' @return
#' an object of class "pNode", a list with following components:
#' \itemize{
#'  \item{\code{priority}: a matrix of nNode X 6 containing node priority information, where nNode is the number of nodes in the input graph, and the 6 columns are "name" (node names), "node" (1 for network genes, 0 for non-network seed genes), "seed" (1 for seeds, 0 for non-seeds), "weight" (weight values),  "priority" (the priority scores that are rescaled to the range [0,1]), "rank" (ranks of the priority scores), "description" (node description)}
#'  \item{\code{g}: an input "igraph" object}
#'  \item{\code{mSeed}: a list with following components 'GR', 'Gene' and 'Link'}
#' }
#' @export
#' @seealso \code{\link{xGR2xGeneScores}}, \code{\link{xPierGenes}}
#' @include xPierGRs.r
#' @examples
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#'
#' \dontrun{
#' # a) provide the seed SNPs with the significance info
#' data(ImmunoBase)
#' ## only AS GWAS SNPs and their significance info (p-values)
#' df <- as.data.frame(ImmunoBase$AS$variant, row.names=NULL)
#' GR <- paste0(df$seqnames,':',df$start,'-',df$end)
#' data <- cbind(GR=GR, Sig=df$Pvalue)
#'
#' # b) perform priority analysis
#' pNode <- xPierGRs(data=data, crosslink="PCHiC_combined", network="STRING_highest", restart=0.7, RData.location=RData.location)
#'
#' # c) save to the file called 'GRs_priority.txt'
#' write.table(pNode$priority, file="GRs_priority.txt", sep="\t", row.names=FALSE)
#' 
#' # d) manhattan plot
#' mp <- xPierManhattan(pNode, top=20, top.label.size=1.5, y.scale="sqrt", RData.location=RData.location)
#' #pdf(file="Gene_manhattan.pdf", height=6, width=12, compress=TRUE)
#' print(mp)
#' #dev.off()
#' }

xPierGRs <- function(data, significance.threshold=NULL, score.cap=NULL, build.conversion=c(NA,"hg38.to.hg19","hg18.to.hg19"), crosslink=c("genehancer","PCHiC_combined","GTEx_V6p_combined","nearby"), crosslink.customised=NULL, cdf.function=c("original","empirical"), scoring.scheme=c("max","sum","sequential"), nearby.distance.max=50000, nearby.decay.kernel=c("rapid","slow","linear","constant"), nearby.decay.exponent=2, network=c("STRING_highest","STRING_high","STRING_medium","STRING_low","PCommonsUN_high","PCommonsUN_medium","PCommonsDN_high","PCommonsDN_medium","PCommonsDN_Reactome","PCommonsDN_KEGG","PCommonsDN_HumanCyc","PCommonsDN_PID","PCommonsDN_PANTHER","PCommonsDN_ReconX","PCommonsDN_TRANSFAC","PCommonsDN_PhosphoSite","PCommonsDN_CTD", "KEGG","KEGG_metabolism","KEGG_genetic","KEGG_environmental","KEGG_cellular","KEGG_organismal","KEGG_disease","REACTOME"), STRING.only=c(NA,"neighborhood_score","fusion_score","cooccurence_score","coexpression_score","experimental_score","database_score","textmining_score")[1], weighted=FALSE, network.customised=NULL, seeds.inclusive=TRUE, normalise=c("laplacian","row","column","none"), restart=0.7, normalise.affinity.matrix=c("none","quantile"), parallel=TRUE, multicores=NULL, verbose=TRUE, RData.location="http://galahad.well.ox.ac.uk/bigdata", guid=NULL)
{

    startT <- Sys.time()
    if(verbose){
        message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=TRUE)
        message("", appendLF=TRUE)
    }
    ####################################################################################
    
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    build.conversion <- match.arg(build.conversion)
    #crosslink <- match.arg(crosslink)
    cdf.function <- match.arg(cdf.function)
    scoring.scheme <- match.arg(scoring.scheme)
    nearby.decay.kernel <- match.arg(nearby.decay.kernel)
    network <- match.arg(network)
    normalise <- match.arg(normalise)
    normalise.affinity.matrix <- match.arg(normalise.affinity.matrix)
    
    ####################################################################################

    if(is.vector(data) && length(data)>1 && is.null(names(data)) && is.character(data)){
    	data <- data.frame(GR=data, Sig=rep(0.1,length(data)), stringsAsFactors=FALSE)
    	significance.threshold <- NULL
    	score.cap <- NULL
    	
		if(verbose){
			message(sprintf("The input has GRs only (without the significance level)"), appendLF=TRUE)
		}
    }

    if(verbose){
        now <- Sys.time()
        message(sprintf("\n#######################################################", appendLF=TRUE))
        message(sprintf("'xGR2xGeneScores' is being called to score seed genes (%s):", as.character(now)), appendLF=TRUE)
        message(sprintf("#######################################################", appendLF=TRUE))
    }
    
    mSeed <- xGR2xGeneScores(data=data, significance.threshold=significance.threshold, score.cap=score.cap, build.conversion=build.conversion, crosslink=crosslink, crosslink.customised=crosslink.customised, cdf.function=cdf.function, scoring.scheme=scoring.scheme, nearby.distance.max=nearby.distance.max, nearby.decay.kernel=nearby.decay.kernel, nearby.decay.exponent=nearby.decay.exponent, verbose=verbose, RData.location=RData.location, guid=guid)
	
	if(verbose){
        now <- Sys.time()
        message(sprintf("#######################################################", appendLF=TRUE))
        message(sprintf("'xGR2xGeneScores' has been finished (%s)!", as.character(now)), appendLF=TRUE)
        message(sprintf("#######################################################\n", appendLF=TRUE))
    }
    
    ######################################################################################
    
    if(verbose){
        now <- Sys.time()
        message(sprintf("\n#######################################################", appendLF=TRUE))
        message(sprintf("'xPierGenes' is being called to prioritise target genes (%s):", as.character(now)), appendLF=TRUE)
        message(sprintf("#######################################################", appendLF=TRUE))
    }
    
    seeds.genes <- mSeed$Gene[,c('Gene','Score')]
    
	pNode <- suppressMessages(xPierGenes(data=seeds.genes, network=network, weighted=weighted, network.customised=network.customised, seeds.inclusive=seeds.inclusive, normalise=normalise, restart=restart, normalise.affinity.matrix=normalise.affinity.matrix, parallel=parallel, multicores=multicores, verbose=verbose, RData.location=RData.location, guid=guid))
	
	if(verbose){
        now <- Sys.time()
        message(sprintf("#######################################################", appendLF=TRUE))
        message(sprintf("'xPierGenes' has been finished (%s)!", as.character(now)), appendLF=TRUE)
        message(sprintf("#######################################################\n", appendLF=TRUE))
    }
    
    #######################
    ## if pNode==NULL, return NULL
    if(is.null(pNode)){
    	return(NULL)
    }
    #######################   
    
	if(verbose){
		now <- Sys.time()
		message(sprintf("A total of %d genes are prioritised, based on:", nrow(pNode$priority)), appendLF=TRUE)
		message(sprintf("\t%d GRs scored positively;", nrow(mSeed$GR)), appendLF=TRUE)
		
		if(is.null(crosslink.customised)){
			if(crosslink=='nearby'){
				message(sprintf("\t%d nearby genes within %d(bp) genomic distance window of %d GRs, scored in a '%s' manner", nrow(mSeed$Gene), nearby.distance.max, nrow(mSeed$GR), nearby.decay.kernel), appendLF=TRUE)
			}
			if(sum(grep("GTEx_V6p_",crosslink,perl=TRUE)) > 0){
				message(sprintf("\t%d eQTL genes with expression modulated by %d GRs", nrow(mSeed$Gene), nrow(mSeed$GR)), appendLF=TRUE)
			}
			if(sum(grep("PCHiC_",crosslink,perl=TRUE)) > 0){
				message(sprintf("\t%d HiC genes physically interacted with %d GRs", nrow(mSeed$Gene), nrow(mSeed$GR)), appendLF=TRUE)
			}

		}
		
		message(sprintf("\t%d genes defined as seeds from %d GRs", nrow(seeds.genes), nrow(mSeed$GR)), appendLF=TRUE)
		message(sprintf("\trandomly walk the network (%d nodes and %d edges) starting from %d seed genes/nodes (with %.2f restarting prob.)", vcount(pNode$g), ecount(pNode$g), nrow(seeds.genes), restart), appendLF=TRUE)
	}
    
    #####
    ## append
    pNode[['mSeed']] <- mSeed
	
    ####################################################################################
    endT <- Sys.time()
    if(verbose){
        message(paste(c("\nFinish at ",as.character(endT)), collapse=""), appendLF=TRUE)
    }
    
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=TRUE)
    
    invisible(pNode)
}
