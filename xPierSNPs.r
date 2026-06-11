PierSNPs <- function(data, include.LD=NA, LD.customised=NULL, LD.r2=0.8, significance.threshold=5e-5, score.cap=10, distance.max=2000, decay.kernel=c("slow","constant","linear","rapid"), decay.exponent=2, GR.SNP=c("dbSNP_GWAS","dbSNP_Common","dbSNP_Single"), GR.Gene=c("UCSC_knownGene","UCSC_knownCanonical"), include.TAD=c("none","GM12878","IMR90","MSC","TRO","H1","MES","NPC"), include.eQTL=NA, eQTL.customised=NULL, include.HiC=NA, cdf.function=c("empirical","exponential"), relative.importance=c(1/3,1/3,1/3), scoring.scheme=c("max","sum","sequential"), network=c("STRING_highest","STRING_high","STRING_medium","STRING_low","PCommonsUN_high","PCommonsUN_medium","PCommonsDN_high","PCommonsDN_medium","PCommonsDN_Reactome","PCommonsDN_KEGG","PCommonsDN_HumanCyc","PCommonsDN_PID","PCommonsDN_PANTHER","PCommonsDN_ReconX","PCommonsDN_TRANSFAC","PCommonsDN_PhosphoSite","PCommonsDN_CTD", "KEGG","KEGG_metabolism","KEGG_genetic","KEGG_environmental","KEGG_cellular","KEGG_organismal","KEGG_disease","REACTOME"), STRING.only=c(NA,"neighborhood_score","fusion_score","cooccurence_score","coexpression_score","experimental_score","database_score","textmining_score")[1], weighted=FALSE, network.customised=NULL, seeds.inclusive=TRUE, normalise=c("laplacian","row","column","none"), restart=0.7, normalise.affinity.matrix=c("none","quantile"), parallel=TRUE, multicores=NULL, verbose=TRUE, RData.location="http://galahad.well.ox.ac.uk/bigdata", guid=NULL)
{

    startT <- Sys.time()
    if(verbose){
        message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=TRUE)
        message("", appendLF=TRUE)
    }
    ####################################################################################
    
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    decay.kernel <- match.arg(decay.kernel)
    # include.TAD <- match.arg(include.TAD)
    cdf.function <- match.arg(cdf.function)
    scoring.scheme <- match.arg(scoring.scheme)
    network <- match.arg(network)
    normalise <- match.arg(normalise)
    normalise.affinity.matrix <- match.arg(normalise.affinity.matrix)
    
    ####################################################################################
    
    if(verbose){
        now <- Sys.time()
        message(sprintf("\n#######################################################", appendLF=TRUE))
        message(sprintf("'xSNPscores' is being called to score SNPs (%s):", as.character(now)), appendLF=TRUE)
        message(sprintf("#######################################################", appendLF=TRUE))
    }
    
	df_SNP <- xSNPscores(data=data, include.LD=include.LD, LD.customised=LD.customised, LD.r2=LD.r2, significance.threshold=significance.threshold, score.cap=score.cap, verbose=verbose, RData.location=RData.location, guid=guid)
	
	if(verbose){
        now <- Sys.time()
        message(sprintf("#######################################################", appendLF=TRUE))
        message(sprintf("'xSNPscores' has been finished (%s)!", as.character(now)), appendLF=TRUE)
        message(sprintf("#######################################################\n", appendLF=TRUE))
    }
    
    ####################################################################################
    
    if(verbose){
        now <- Sys.time()
        message(sprintf("\n#######################################################", appendLF=TRUE))
        message(sprintf("'xSNP2nGenes' is being called to define nearby genes (%s):", as.character(now)), appendLF=TRUE)
        message(sprintf("#######################################################", appendLF=TRUE))
    }
    
    
    if(relative.importance[1] != 0){
		df_nGenes <- xSNP2nGenes(data=df_SNP$SNP, distance.max=distance.max, decay.kernel=decay.kernel, decay.exponent=decay.exponent, GR.SNP=GR.SNP, GR.Gene=GR.Gene, include.TAD=include.TAD, verbose=verbose, RData.location=RData.location, guid=guid)
		if(length(include.TAD)>0){
			TAD <- NULL
			df_nGenes <- base::subset(df_nGenes, TAD!='Excluded')
		}
	}else{
		df_nGenes <- NULL
		if(verbose){
			now <- Sys.time()
			message(sprintf("No nearby genes are defined"), appendLF=TRUE)
		}
	}
	
	if(verbose){
        now <- Sys.time()
        message(sprintf("#######################################################", appendLF=TRUE))
        message(sprintf("'xSNP2nGenes' has been finished (%s)!", as.character(now)), appendLF=TRUE)
        message(sprintf("#######################################################\n", appendLF=TRUE))
    }
    
    ####################################################################################
    
    if(verbose){
        now <- Sys.time()
        message(sprintf("\n#######################################################", appendLF=TRUE))
        message(sprintf("'xSNP2eGenes' is being called to define eQTL-containing genes (%s):", as.character(now)), appendLF=TRUE)
        message(sprintf("#######################################################", appendLF=TRUE))
    }
    
    if(relative.importance[2] != 0){
		df_eGenes <- xSNP2eGenes(data=df_SNP$SNP, include.eQTL=include.eQTL, eQTL.customised=eQTL.customised, cdf.function=cdf.function, plot=FALSE, verbose=verbose, RData.location=RData.location, guid=guid)
	}else{
		df_eGenes <- NULL
		
		if(verbose){
			now <- Sys.time()
			message(sprintf("No eQTL genes are defined"), appendLF=TRUE)
		}
	}
	
	if(verbose){
        now <- Sys.time()
        message(sprintf("#######################################################", appendLF=TRUE))
        message(sprintf("'xSNP2eGenes' has been finished (%s)!", as.character(now)), appendLF=TRUE)
        message(sprintf("#######################################################\n", appendLF=TRUE))
    }
    ####################################################################################
    
    ####################################################################################
    
    if(verbose){
        now <- Sys.time()
        message(sprintf("\n#######################################################", appendLF=TRUE))
        message(sprintf("'xSNP2cGenes' is being called to define HiC-captured genes (%s):", as.character(now)), appendLF=TRUE)
        message(sprintf("#######################################################", appendLF=TRUE))
    }
    
    if(relative.importance[3] != 0){
		df_cGenes <- xSNP2cGenes(data=df_SNP$SNP, entity="SNP", include.HiC=include.HiC, GR.SNP=GR.SNP, cdf.function=cdf.function, plot=FALSE, verbose=verbose, RData.location=RData.location, guid=guid)
	}else{
		df_cGenes <- NULL
		
		if(verbose){
			now <- Sys.time()
			message(sprintf("No HiC genes are defined"), appendLF=TRUE)
		}
	}
	
	if(verbose){
        now <- Sys.time()
        message(sprintf("#######################################################", appendLF=TRUE))
        message(sprintf("'xSNP2cGenes' has been finished (%s)!", as.character(now)), appendLF=TRUE)
        message(sprintf("#######################################################\n", appendLF=TRUE))
    }
    ####################################################################################
    
    if(is.null(df_nGenes) & is.null(df_eGenes) & is.null(df_cGenes)){
    	G2S <- NULL
    }else{
    
		## df_SNP df_nGenes df_eGenes df_cGenes
		allGenes <- sort(base::Reduce(base::union, list(df_nGenes$Gene,df_eGenes$Gene,df_cGenes$Gene)))
		allSNPs <- sort(df_SNP$SNP)
	
		## sparse matrix of nGenes X SNPs
		G2S_n <- xSparseMatrix(df_nGenes[,c("Gene","SNP","Weight")], rows=allGenes, columns=allSNPs, verbose=FALSE)
		## sparse matrix of eGenes X SNPs
		G2S_e <- xSparseMatrix(df_eGenes[,c("Gene","SNP","Weight")], rows=allGenes, columns=allSNPs, verbose=FALSE)
		## sparse matrix of cGenes X SNPs
		G2S_c <- xSparseMatrix(df_cGenes[,c("Gene","SNP","Weight")], rows=allGenes, columns=allSNPs, verbose=FALSE)
	
		## combine both sparse matrix
		### wG2S_n
		if(is.null(G2S_n)){
			wG2S_n <- 0
		}else{
			wG2S_n <- G2S_n * relative.importance[1]
		}
		### wG2S_e
		if(is.null(G2S_e)){
			wG2S_e <- 0
		}else{
			wG2S_e <- G2S_e * relative.importance[2]
		}
		### wG2S_c
		if(is.null(G2S_c)){
			wG2S_c <- 0
		}else{
			wG2S_c <- G2S_c * relative.importance[3]
		}
		
		if(is.null(G2S_n) & is.null(G2S_e) & is.null(G2S_c)){
			G2S <- NULL
		}else{
			G2S <- wG2S_n + wG2S_e + wG2S_c
		}
    
    }
    
    #######################
    ## if NULL, return NULL
    if(is.null(G2S)){
    	return(NULL)
    }
    #######################
    
    ## consider SNP scores
    ind <- match(colnames(G2S), df_SNP$SNP)
    ########
    df_SNP <- df_SNP[ind,]
    ########
    SNP_score <- df_SNP$Score
    names(SNP_score) <- colnames(G2S)
    ## convert into matrix
    mat_SNP_score <- matrix(rep(SNP_score,each=nrow(G2S)), nrow=nrow(G2S))
    
    ## calculate genetic influence score for a gene-SNP pair
    G2S_score <- G2S * mat_SNP_score
    
    ## Gene2SNP
    Gene2SNP <- xSM2DF(data=G2S_score, verbose=FALSE)
    colnames(Gene2SNP) <- c('Gene','SNP','Score')
    Gene2SNP <- Gene2SNP[order(Gene2SNP$Gene,-Gene2SNP$Score,decreasing=FALSE),]
	
	ls_gene <- split(x=Gene2SNP$Score, f=Gene2SNP$Gene)
    ## calculate genetic influence score under a set of SNPs for each seed gene
    if(scoring.scheme=='max'){
		seeds.genes <- sapply(ls_gene, max)
		
    }else if(scoring.scheme=='sum'){
		seeds.genes <- sapply(ls_gene, sum)
		
    }else if(scoring.scheme=='sequential'){
		seeds.genes <- sapply(ls_gene, function(x){
			base::sum(x / base::rank(-x,ties.method="min"))
		})
		
    }
	
	if(verbose){
		now <- Sys.time()
		message(sprintf("%d Genes are defined as seeds and scored using '%s' scoring scheme from %d SNPs", length(seeds.genes), scoring.scheme, ncol(G2S_score)), appendLF=TRUE)
	}
    
    ######################################################################################
    
    if(verbose){
        now <- Sys.time()
        message(sprintf("\n#######################################################", appendLF=TRUE))
        message(sprintf("'xPierGenes' is being called to prioritise target genes (%s):", as.character(now)), appendLF=TRUE)
        message(sprintf("#######################################################", appendLF=TRUE))
    }
    
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
		message(sprintf("\t%d SNPs scored positively (including %d 'Lead' and %d 'LD');", nrow(df_SNP), sum(df_SNP$Flag=='Lead'), sum(df_SNP$Flag=='LD')), appendLF=TRUE)
		if(!is.null(df_nGenes)){
			message(sprintf("\t%d nearby genes within %d(bp) genomic distance window of %d SNPs", length(unique(df_nGenes$Gene)), distance.max, length(unique(df_nGenes$SNP))), appendLF=TRUE)
		}
		if(!is.null(df_eGenes)){
			message(sprintf("\t%d eQTL genes with expression modulated by %d SNPs", length(unique(df_eGenes$Gene)), length(unique(df_eGenes$SNP))), appendLF=TRUE)
		}
		if(!is.null(df_cGenes)){
			message(sprintf("\t%d HiC genes physically interacted with %d SNP", length(unique(df_cGenes$Gene)), length(unique(df_cGenes$SNP))), appendLF=TRUE)
		}
		message(sprintf("\t%d genes defined as seeds from %d SNPs", length(seeds.genes), ncol(G2S_score)), appendLF=TRUE)
		message(sprintf("\trandomly walk the network (%d nodes and %d edges) starting from %d seed genes/nodes (with %.2f restarting prob.)", vcount(pNode$g), ecount(pNode$g), length(seeds.genes), restart), appendLF=TRUE)
	}
    
    #####
    ## SNP
    df_SNP <- df_SNP[order(df_SNP$Flag,df_SNP$Score,df_SNP$SNP,decreasing=TRUE),]
    
    ## nGenes
    if(is.null(df_nGenes)){
    	nGenes <- NULL
    }else{
		nGenes <- df_nGenes
		ind <- match(nGenes$SNP, df_SNP$SNP)
		nGenes$SNP_Flag <- df_SNP$Flag[ind]
	}
    ## eGenes
    if(is.null(df_eGenes)){
    	eGenes <- NULL
    }else{
		eGenes <- xDefineEQTL(data=df_SNP$SNP, include.eQTL=include.eQTL, eQTL.customised=eQTL.customised, verbose=FALSE, RData.location=RData.location, guid=guid)
		ind <- match(eGenes$SNP, df_SNP$SNP)
		eGenes$SNP_Flag <- df_SNP$Flag[ind]
	}
    ## cGenes
    if(is.null(df_cGenes)){
    	cGenes <- NULL
    }else{
		cGenes <- xDefineHIC(data=df_SNP$SNP, entity="SNP", include.HiC=include.HiC, GR.SNP=GR.SNP, verbose=FALSE, RData.location=RData.location, guid=guid)
		cGenes <- cGenes$df
		ind <- match(cGenes$SNP, df_SNP$SNP)
		cGenes$SNP_Flag <- df_SNP$Flag[ind]
    }
    
    #####
    ## append
    pNode[['SNP']] <- df_SNP
    pNode[['Gene2SNP']] <- Gene2SNP
    pNode[['nGenes']] <- nGenes
    pNode[['eGenes']] <- eGenes
    pNode[['cGenes']] <- cGenes

	
    ####################################################################################
    endT <- Sys.time()
    if(verbose){
        message(paste(c("\nFinish at ",as.character(endT)), collapse=""), appendLF=TRUE)
    }
    
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=TRUE)
    
    invisible(pNode)
}
