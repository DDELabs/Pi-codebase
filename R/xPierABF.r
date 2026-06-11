#' Function to prioritise genes based on seed eGenes identified through ABF integrating GWAS and eQTL summary data
#'
#' \code{xPierABF} is supposed to prioritise genes based on seed eGenes identified through ABF integrating GWAS and eQTL summary data. To prioritise genes, it first conducts colocalisation analysis through Wakefield's Approximate Bayes Factor (ABF) integrating GWAS and eQTL summary data to identify and score seed genes (that is, eGenes weighted by posterior probability of the SNP being causal for both GWAS and eQTL traits). It implements Random Walk with Restart (RWR) and calculates the affinity score of all nodes in the graph to the seeds. The priority score is the affinity score. Parallel computing is also supported for Linux-like or Windows operating systems. It returns an object of class "pNode". 
#'
#' @param data a data frame storing GWAS summary data with following required columns 'snp', 'effect' (the effect allele assessed), 'other' (other allele), 'b' (effect size for the allele assessed; log(odds ratio) for a case-control study), 'se' (standard error), 'p' (p-value)
#' @param eqtl context-specific eQTL summary data. It can be one of "Bcell","Blood","CD14","CD4","CD8","IFN","LPS24","LPS2","Monocyte","Neutrophil","NK","shared_CD14","shared_IFN","shared_LPS24","shared_LPS2"
#' @param prior.eqtl the prior probability an eQTL associated with the eQTL trait. The default value is 1e-4
#' @param prior.gwas the prior probability an SNP associated with the GWAS trait. The default value is 1e-4
#' @param prior.both the prior probability an eQTL/SNP associated with both eQTL/GWAS traits. The default value is 1e-5
#' @param cutoff.H4 the H4 cutoff used to define eGenes. This cutoff is based on the posterior probabilities of H4 - one shared causal variant. The default value is 0.8
#' @param cutoff.pgwas the GWAS p-value cutoff that must be met to consider SNPs. The default value is 1e-5
#' @param network the built-in network. Currently two sources of network information are supported: the STRING database (version 10) and the Pathways Commons database (version 7). STRING is a meta-integration of undirect interactions from the functional aspect, while Pathways Commons mainly contains both undirect and direct interactions from the physical/pathway aspect. Both have scores to control the confidence of interactions. Therefore, the user can choose the different quality of the interactions. In STRING, "STRING_highest" indicates interactions with highest confidence (confidence scores>=900), "STRING_high" for interactions with high confidence (confidence scores>=700), "STRING_medium" for interactions with medium confidence (confidence scores>=400), and "STRING_low" for interactions with low confidence (confidence scores>=150). For undirect/physical interactions from Pathways Commons, "PCommonsUN_high" indicates undirect interactions with high confidence (supported with the PubMed references plus at least 2 different sources), "PCommonsUN_medium" for undirect interactions with medium confidence (supported with the PubMed references). For direct (pathway-merged) interactions from Pathways Commons, "PCommonsDN_high" indicates direct interactions with high confidence (supported with the PubMed references plus at least 2 different sources), and "PCommonsUN_medium" for direct interactions with medium confidence (supported with the PubMed references). In addtion to pooled version of pathways from all data sources, the user can also choose the pathway-merged network from individual sources, that is, "PCommonsDN_Reactome" for those from Reactome, "PCommonsDN_KEGG" for those from KEGG, "PCommonsDN_HumanCyc" for those from HumanCyc, "PCommonsDN_PID" for those froom PID, "PCommonsDN_PANTHER" for those from PANTHER, "PCommonsDN_ReconX" for those from ReconX, "PCommonsDN_TRANSFAC" for those from TRANSFAC, "PCommonsDN_PhosphoSite" for those from PhosphoSite, and "PCommonsDN_CTD" for those from CTD. For direct (pathway-merged) interactions sourced from KEGG, it can be 'KEGG' for all, 'KEGG_metabolism' for pathways grouped into 'Metabolism', 'KEGG_genetic' for 'Genetic Information Processing' pathways, 'KEGG_environmental' for 'Environmental Information Processing' pathways, 'KEGG_cellular' for 'Cellular Processes' pathways, 'KEGG_organismal' for 'Organismal Systems' pathways, and 'KEGG_disease' for 'Human Diseases' pathways
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
#'  \item{\code{priority}: a matrix of nNode X 6 containing node priority information, where nNode is the number of nodes in the input graph, and the 5 columns are "name" (node names), "node" (1 for network genes, 0 for non-network seed genes), "seed" (1 for seeds, 0 for non-seeds), "weight" (weight values),  "priority" (the priority scores that are rescaled to the range [0,1]), "rank" (ranks of the priority scores), "description" (node description)}
#'  \item{\code{g}: an input "igraph" object}
#'  \item{\code{evidence}: a data frame storing evidence}
#'  \item{\code{call}: the call that produced this result}
#' }
#' @note The input graph will treat as an unweighted graph if there is no 'weight' edge attribute associated with
#' @export
#' @seealso \code{\link{xRDataLoader}}, \code{\link{xMEabf}}, \code{\link{xPierGenes}}
#' @include xPierABF.r
#' @examples
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' \dontrun{
#' data <- utils::read.delim(file="summary_gwas.RA.txt", header=TRUE, row.names=NULL, stringsAsFactors=FALSE)
#' pNode_abf <- xPierABF(data, eqtl="Blood", network="STRING_high", restart=0.7, RData.location=RData.location)
#' write.table(pNode_abf$priority, file="Genes_priority.ABF.txt", sep="\t", row.names=FALSE)
#' }

xPierABF <- function(data, eqtl=c("CD14","LPS2","LPS24","IFN","Bcell","NK","Neutrophil","CD4","CD8","Blood","Monocyte","shared_CD14","shared_LPS2","shared_LPS24","shared_IFN"), prior.eqtl=1e-4, prior.gwas=1e-4, prior.both=1e-5, cutoff.H4=0.8, cutoff.pgwas=1e-5, network=c("STRING_highest","STRING_high","STRING_medium","STRING_low","PCommonsUN_high","PCommonsUN_medium","PCommonsDN_high","PCommonsDN_medium","PCommonsDN_Reactome","PCommonsDN_KEGG","PCommonsDN_HumanCyc","PCommonsDN_PID","PCommonsDN_PANTHER","PCommonsDN_ReconX","PCommonsDN_TRANSFAC","PCommonsDN_PhosphoSite","PCommonsDN_CTD","KEGG","KEGG_metabolism","KEGG_genetic","KEGG_environmental","KEGG_cellular","KEGG_organismal","KEGG_disease"), STRING.only=c(NA,"neighborhood_score","fusion_score","cooccurence_score","coexpression_score","experimental_score","database_score","textmining_score")[1], weighted=FALSE, network.customised=NULL, seeds.inclusive=TRUE, normalise=c("laplacian","row","column","none"), restart=0.7, normalise.affinity.matrix=c("none","quantile"), parallel=TRUE, multicores=NULL, verbose=TRUE, RData.location="http://galahad.well.ox.ac.uk/bigdata", guid=NULL)
{

    startT <- Sys.time()
    if(verbose){
        message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=TRUE)
        message("", appendLF=TRUE)
    }
    ####################################################################################
    
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    eqtl <- match.arg(eqtl)
    network <- match.arg(network)
    normalise <- match.arg(normalise)
    normalise.affinity.matrix <- match.arg(normalise.affinity.matrix)
    
    pNode <- NULL
    se <- b <- effect <- other <- p <- NULL
    pp_ABF <- p_GWAS <- p_eQTL <- H4 <- NULL
	
	###########################	
	# summary GWAS
	###########################
	if(is(data,'data.frame')){
		if(all(c('snp','effect','other','b','p','se') %in% colnames(data))){
			#summary_gwas <- data[,c("snp","effect","other","b","p","se")] %>% dplyr::filter(!is.na(se), se!=0, b!=0, effect!='', other!='')
			summary_gwas <- data %>% dplyr::filter(!is.na(se), se!=0, b!=0, effect!='', other!='') %>% dplyr::arrange(p)
			################
			# IMPORTANT!: keep the most significant if multiple entries for the same SNP
			summary_gwas <- summary_gwas[!duplicated(summary_gwas$snp),]
			################
			
		}else{
			warnings("The input data.frame does not contain required columns: c('snp','effect','other','b','p','se').\n")
			return(NULL)
		}
	}
	
	###########################
	# built-in summary eQTL
	###########################
	## number of samples analysed
	vec_N_eqtl <- c(414, 261, 322, 367, 286, 245, 101, 293, 283, 5311, 287, 228, 228, 228, 228)
names(vec_N_eqtl) <- c("CD14","LPS2","LPS24","IFN","Bcell","NK","Neutrophil","CD4","CD8","Blood","Monocyte","shared_CD14","shared_LPS2","shared_LPS24","shared_IFN")
	N_eqtl <- vec_N_eqtl[eqtl]
	if(0){
		JK_cohort_xMEdb <- xRDataLoader('JK_cohort_xMEdb', verbose=FALSE, RData.location=RData.location, guid=guid)
		## summary_eqtl: extracted 
		ind <- match(JK_cohort_xMEdb$context, eqtl)
		summary_eqtl <- JK_cohort_xMEdb[!is.na(ind), ]
	}else{
		summary_eqtl <- xRDataLoader(paste0('JK_cohort_xMEdb_',eqtl), verbose=FALSE, RData.location=RData.location, guid=guid)
	}
	
	##########################
	## much quick!
	ind <- match(summary_gwas$snp, summary_eqtl$snps)
	summary_gwas <- summary_gwas[!is.na(ind), ]
	##########################	
	
	## analysis per probe
	ls_gene <- split(x=summary_eqtl, f=summary_eqtl$gene)
	ls_df_output <- lapply(1:length(ls_gene), function(i){
		
		#eqtl <- 'Blood'
		#i <- match(4570768, names(ls_gene))
		
		if(verbose){
			if(i %% 1000 == 0){
				message(sprintf("Analysing %d (%d) (%s)", i, length(ls_gene), as.character(Sys.time())), appendLF=TRUE)
			}
			message(sprintf("Analysing %d (%d) (%s)", i, length(ls_gene), as.character(Sys.time())), appendLF=TRUE)
		}
		
		df_output <- NULL
		
		## SNPs analysed per probe
		df_eqtl <- ls_gene[[i]]

		## common SNPs analysed
		ind <- match(summary_gwas$snp, df_eqtl$snps)
		if(sum(!is.na(ind))>=3){
			df_gwas <- summary_gwas[!is.na(ind), ]
			df_eqtl <- df_eqtl[ind[!is.na(ind)], ]

			## double-check the effect allele in GWAS, fixing the effect allele in eQTL
			df_gwas$b_corrected <- df_gwas$b
			## the effect allele reversed
			ind <- which(df_gwas$effect==df_eqtl$other_allele & df_gwas$other==df_eqtl$effect_allele)
			df_gwas$b_corrected[ind] <- -1 * df_gwas$b[ind]
			## remove SNPs if alleles disagreed
			ind <- df_gwas$effect!=df_eqtl$effect_allele & df_gwas$effect!=df_eqtl$other_allele
			if(sum(ind)>0){
				df_gwas <- df_gwas[!ind, ]
				df_eqtl <- df_eqtl[!ind, ]
			}
			
			if(nrow(df_gwas)>=3){
				
				## eqtl.summary ('beta', 'varbeta', 'N', 'MAF', 'snp') required and extracted from 'df_eqtl'
				eqtl.summary <- list(beta=df_eqtl$beta, varbeta=df_eqtl$se^2, N=N_eqtl, MAF=df_eqtl$effect_maf, type="quant", snp=df_eqtl$snps)
				
				## gwas.summary ('beta', 'varbeta', 'snp') required and extracted from 'df_gwas'
				gwas.summary <- list(beta=df_gwas$b_corrected, varbeta=df_gwas$se^2, type="cc", snp=df_gwas$snp)
		
				## Bayesian colocalisation analysis through ABF
				res <- xMEabf(eqtl.summary, gwas.summary, prior.eqtl=prior.eqtl, prior.gwas=prior.gwas, prior.both=prior.both)
				
				## post-processing to obtain df_output
				df_res <- res$results
				## append df_eqtl and df_gwas
				ind <- match(df_res$snp, df_eqtl$snps)
				df_res <- cbind(df_res, df_eqtl[ind,], df_gwas[ind,])
				## df_output
				df_output <- df_res[,c('context','mode','ProbeID','Symbol','gene_cse','snps','snp_cse','effect_allele','other_allele','b_corrected','beta','p','pvalue','SNP.PP.H4')]
				colnames(df_output) <- c('context','mode','ProbeID','Symbol','gene_cse','snps','snp_cse','A1','A2','b_GWAS','b_eQTL','p_GWAS','p_eQTL','pp_ABF')
				## append b_ABF
				df_output$b_ABF <- df_output$b_GWAS / df_output$b_eQTL
				## sort by pp_ABF (that is, SNP.PP.H4)
				df_output <- df_output %>% dplyr::arrange(-pp_ABF, p_GWAS, p_eQTL)

				## append: nsnps, H0, H1, H2, H3, H4
				df_output$nsnps <- res$summary[1]
				df_output$H0 <- res$summary[2]
				df_output$H1 <- res$summary[3]
				df_output$H2 <- res$summary[4]
				df_output$H3 <- res$summary[5]
				df_output$H4 <- res$summary[6]
			}
		}
		df_output
	})
	df_output <- do.call(rbind, ls_df_output)

	## Keep the top SNP (with the highest pp_ABF, and then the most significant p_GWAS and p_eQTL) per ProbeID
	ind <- which(!duplicated(df_output$ProbeID))
	df <- df_output[ind,]

	## ProbeID-level cutoff: H4 > 0.8
	## ProbeID-SNP level cutoff: p_GWAS > 1e-5
	df <- subset(df, H4>cutoff.H4 & p_GWAS<cutoff.pgwas)
	
	##############
	## df_evidence
	df_evidence <- df %>% dplyr::arrange(-H4, p_GWAS, p_eQTL)
	#utils::write.table(df_evidence, file="df_output_ABF.txt", row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
	##############
	
	## Keep the top SNP (with the highest H4, the most significant p_GWAS and p_eQTL) per Symbol
	df <- df %>% dplyr::arrange(-H4, p_GWAS, p_eQTL)
	ind <- which(!duplicated(df$Symbol))
	df <- df[ind,]
	#sort(df$Symbol)
	#df[,c('Symbol','mode','nsnps')] %>% dplyr::arrange(mode,nsnps)
	
	## the seed gene (eGene) weighted by pp_ABF
	data_subset <- df[,c("Symbol","pp_ABF")]
	if(nrow(data_subset)!=0){
		pNode <- suppressMessages(xPierGenes(data=data_subset, network=network, STRING.only=STRING.only, weighted=weighted, network.customised=network.customised, seeds.inclusive=seeds.inclusive, normalise=normalise, restart=restart, normalise.affinity.matrix=normalise.affinity.matrix, parallel=parallel, multicores=multicores, verbose=verbose, RData.location=RData.location, guid=guid))
		
		pNode$evidence <- df_evidence
    }
    
    ####################################################################################
    endT <- Sys.time()
    if(verbose){
        message(paste(c("\nFinish at ",as.character(endT)), collapse=""), appendLF=TRUE)
    }
    
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=TRUE)
    
    invisible(pNode)
}
