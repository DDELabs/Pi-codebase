PierSNPsAdv <- function(data, include.LD=NA, LD.customised=NULL, LD.r2=0.8, significance.threshold=5e-5, score.cap=10, distance.max=2000, decay.kernel=c("slow","constant","linear","rapid"), decay.exponent=2, GR.SNP=c("dbSNP_GWAS","dbSNP_Common","dbSNP_Single"), GR.Gene=c("UCSC_knownGene","UCSC_knownCanonical"), include.TAD=c("none","GM12878","IMR90","MSC","TRO","H1","MES","NPC"), include.eQTL=NA, eQTL.customised=NULL, include.HiC=NA, cdf.function=c("empirical","exponential"), scoring.scheme=c("max","sum","sequential"), network=c("STRING_highest","STRING_high","STRING_medium","STRING_low","PCommonsUN_high","PCommonsUN_medium","PCommonsDN_high","PCommonsDN_medium","PCommonsDN_Reactome","PCommonsDN_KEGG","PCommonsDN_HumanCyc","PCommonsDN_PID","PCommonsDN_PANTHER","PCommonsDN_ReconX","PCommonsDN_TRANSFAC","PCommonsDN_PhosphoSite","PCommonsDN_CTD", "KEGG","KEGG_metabolism","KEGG_genetic","KEGG_environmental","KEGG_cellular","KEGG_organismal","KEGG_disease","REACTOME"), STRING.only=c(NA,"neighborhood_score","fusion_score","cooccurence_score","coexpression_score","experimental_score","database_score","textmining_score")[1], weighted=FALSE, network.customised=NULL, seeds.inclusive=TRUE, normalise=c("laplacian","row","column","none"), restart=0.7, normalise.affinity.matrix=c("none","quantile"), parallel=TRUE, multicores=NULL, verbose=TRUE, verbose.details=FALSE, RData.location="http://galahad.well.ox.ac.uk/bigdata", guid=NULL)
{

    startT <- Sys.time()
    if(verbose){
        message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=TRUE)
        message("", appendLF=TRUE)
    }
    ####################################################################################
    
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    decay.kernel <- match.arg(decay.kernel)
    cdf.function <- match.arg(cdf.function)
    scoring.scheme <- match.arg(scoring.scheme)
    network <- match.arg(network)
    normalise <- match.arg(normalise)
    normalise.affinity.matrix <- match.arg(normalise.affinity.matrix)
    
    ## force verbose.details to be FALSE if verbose is FALSE
    if(verbose==FALSE){
    	verbose.details <- FALSE
    }
    ####################################################################################
	if(verbose){
		now <- Sys.time()
		message(sprintf("Preparing the distance predictor (%s) ...", as.character(now)), appendLF=TRUE)
	}
	relative.importance <- c(1,0,0)
    pNode_distance <- PierSNPs(data=data, include.LD=include.LD, LD.customised=LD.customised, LD.r2=LD.r2, significance.threshold=significance.threshold, score.cap=score.cap, distance.max=distance.max, decay.kernel=decay.kernel, decay.exponent=decay.exponent, GR.SNP=GR.SNP, GR.Gene=GR.Gene, include.TAD=include.TAD, include.eQTL=NA, eQTL.customised=NULL, include.HiC=NA, cdf.function=cdf.function, relative.importance=relative.importance, scoring.scheme=scoring.scheme, network=network, weighted=weighted, network.customised=network.customised, seeds.inclusive=seeds.inclusive, normalise=normalise, restart=restart, normalise.affinity.matrix=normalise.affinity.matrix, parallel=parallel, multicores=multicores, verbose=verbose.details, RData.location=RData.location, guid=guid)
    ls_pNode_distance <- list(pNode_distance)
    names(ls_pNode_distance) <- paste('nGene_', distance.max, '_', decay.kernel, sep='')
    
    ####################################################################################
    
    ls_pNode_eQTL <- NULL
    
    include.eQTLs <- include.eQTL[!is.na(include.eQTL)]
    if(length(include.eQTLs)>0){
		names(include.eQTLs) <- include.eQTLs
		ls_pNode_eQTL <- lapply(include.eQTLs, function(x){
			if(verbose){
				now <- Sys.time()
				message(sprintf("Preparing the eQTL predictor '%s' (%s) ...", x, as.character(now)), appendLF=TRUE)
			}
			relative.importance <- c(0,1,0)
			pNode <- PierSNPs(data=data, include.LD=include.LD, LD.customised=LD.customised, LD.r2=LD.r2, significance.threshold=significance.threshold, score.cap=score.cap, distance.max=distance.max, decay.kernel=decay.kernel, decay.exponent=decay.exponent, GR.SNP=GR.SNP, GR.Gene=GR.Gene, include.eQTL=x, eQTL.customised=NULL, include.HiC=NA, cdf.function=cdf.function, relative.importance=relative.importance, scoring.scheme=scoring.scheme, network=network, weighted=weighted, network.customised=network.customised, seeds.inclusive=seeds.inclusive, normalise=normalise, restart=restart, normalise.affinity.matrix=normalise.affinity.matrix, parallel=parallel, multicores=multicores, verbose=verbose.details, RData.location=RData.location, guid=guid)
			if(verbose & is.null(pNode)){
				message(sprintf("\tNote: this predictor '%s' is NULL", x), appendLF=TRUE)
			}
			return(pNode)
		})
		names(ls_pNode_eQTL) <- paste('eGene_', names(ls_pNode_eQTL), sep='')
    }
    
    ################################
    ################################
    ls_pNode_eQTL_customised <- NULL
    df_SGS_customised <- NULL
    if(!is.null(eQTL.customised)){
    
		if(is.vector(eQTL.customised)){
			# assume a file
			df <- utils::read.delim(file=eQTL.customised, header=TRUE, row.names=NULL, stringsAsFactors=FALSE)
		}else if(is.matrix(eQTL.customised) | is.data.frame(eQTL.customised)){
			df <- eQTL.customised
		}
		
		if(!is.null(df)){
			colnames(df) <- c("SNP", "Gene", "Sig", "Context")
			SGS_customised <- df
			#SGS_customised <- cbind(df, Context=rep('Customised',nrow(df)))
			
			############################
			# remove Gene if NA
			# remove SNP if NA
			df_SGS_customised <- SGS_customised[!is.na(SGS_customised[,1]) & !is.na(SGS_customised[,2]),]
			############################
		}
    }
    if(!is.null(df_SGS_customised)){
		ls_df <- split(x=df_SGS_customised, f=df_SGS_customised$Context)
		ls_pNode_eQTL_customised <- lapply(1:length(ls_df), function(i){
			if(verbose){
				now <- Sys.time()
				message(sprintf("Preparing the customised eQTL predictor '%s' (%s) ...", names(ls_df)[i], as.character(now)), appendLF=TRUE)
			}
			relative.importance <- c(0,1,0)
			pNode <- PierSNPs(data=data, include.LD=include.LD, LD.customised=LD.customised, LD.r2=LD.r2, significance.threshold=significance.threshold, score.cap=score.cap, distance.max=distance.max, decay.kernel=decay.kernel, decay.exponent=decay.exponent, GR.SNP=GR.SNP, GR.Gene=GR.Gene, include.eQTL=NA, eQTL.customised=ls_df[[i]], include.HiC=NA, cdf.function=cdf.function, relative.importance=relative.importance, scoring.scheme=scoring.scheme, network=network, weighted=weighted, network.customised=network.customised, seeds.inclusive=seeds.inclusive, normalise=normalise, restart=restart, normalise.affinity.matrix=normalise.affinity.matrix, parallel=parallel, multicores=multicores, verbose=verbose.details, RData.location=RData.location, guid=guid)
			if(verbose & is.null(pNode)){
				message(sprintf("\tNote: this predictor '%s' is NULL", names(ls_df)[i]), appendLF=TRUE)
			}
			return(pNode)
		})
		names(ls_pNode_eQTL_customised) <- paste('eGene_', names(ls_df), sep='')
    }
    ls_pNode_eQTL <- c(ls_pNode_eQTL, ls_pNode_eQTL_customised)
    ################################
    ################################
    
    include.HiCs <- include.HiC[!is.na(include.HiC)]
    if(length(include.HiCs)>0){
		names(include.HiCs) <- include.HiCs
		ls_pNode_HiC <- lapply(include.HiCs, function(x){
			if(verbose){
				now <- Sys.time()
				message(sprintf("Preparing the HiC predictor '%s' (%s) ...", x, as.character(now)), appendLF=TRUE)
			}
			relative.importance <- c(0,0,1)
			pNode <- PierSNPs(data=data, include.LD=include.LD, LD.customised=LD.customised, LD.r2=LD.r2, significance.threshold=significance.threshold, score.cap=score.cap, distance.max=distance.max, decay.kernel=decay.kernel, decay.exponent=decay.exponent, GR.SNP=GR.SNP, GR.Gene=GR.Gene, include.eQTL=NA, eQTL.customised=NULL, include.HiC=x, cdf.function=cdf.function, relative.importance=relative.importance, scoring.scheme=scoring.scheme, network=network, weighted=weighted, network.customised=network.customised, seeds.inclusive=seeds.inclusive, normalise=normalise, restart=restart, normalise.affinity.matrix=normalise.affinity.matrix, parallel=parallel, multicores=multicores, verbose=verbose.details, RData.location=RData.location, guid=guid)
			if(verbose & is.null(pNode)){
				message(sprintf("\tNote: this predictor '%s' has NULL", x), appendLF=TRUE)
			}
			return(pNode)
		})
		names(ls_pNode_HiC) <- paste('cGene_', names(ls_pNode_HiC), sep='')
	}else{
		ls_pNode_HiC <- NULL
	}
    
    ##########################################################################################
    ## prioritisation equally
    #relative.importance <- c(1/3,1/3,1/3)
    #pNode_all <- xPierSNPs(data=data, include.LD=include.LD, LD.customised=LD.customised, LD.r2=LD.r2, significance.threshold=significance.threshold, score.cap=score.cap, distance.max=distance.max, decay.kernel=decay.kernel, decay.exponent=decay.exponent, GR.SNP=GR.SNP, GR.Gene=GR.Gene, include.eQTL=include.eQTLs, eQTL.customised=NULL, include.HiC=include.HiCs, cdf.function=cdf.function, relative.importance=relative.importance, scoring.scheme=scoring.scheme, network=network, weighted=weighted, network.customised=network.customised, seeds.inclusive=seeds.inclusive, normalise=normalise, restart=restart, normalise.affinity.matrix=normalise.affinity.matrix, parallel=parallel, multicores=multicores, verbose=verbose, RData.location=RData.location, guid=guid)
    ##########################################################################################
    ls_pNode <- c(ls_pNode_distance, ls_pNode_eQTL, ls_pNode_HiC)
    
    ####################################################################################
    endT <- Sys.time()
    if(verbose){
        message(paste(c("\nFinish at ",as.character(endT)), collapse=""), appendLF=TRUE)
    }
    
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=TRUE)
    
    invisible(ls_pNode)
}
