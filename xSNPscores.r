xSNPscores <- function(data, include.LD=NA, LD.customised=NULL, LD.r2=0.5, significance.threshold=5e-5, score.cap=10, verbose=TRUE, RData.location="http://galahad.well.ox.ac.uk/bigdata", guid=NULL)
{

    if(is.null(data)){
        stop("The input data must be not NULL.\n")
    }else{
    	
    	if(is(data,'DataFrame')){
    		#data <- S4Vectors::as.matrix(data)
    	}
    
		if (is.vector(data)){
			if(length(data)>1){
				# assume a vector
				if(is.null(names(data))){
					stop("The input data must have names with attached dbSNP ID.\n")
				}
			}else{
				# assume a file
				data <- utils::read.delim(file=data, header=FALSE, row.names=NULL, stringsAsFactors=FALSE)
			}
		}
		
		if (is.vector(data)){
			pval <- data
		}else if(is.matrix(data) | is.data.frame(data)){
			data <- as.matrix(data)
			data_list <- split(x=data[,2], f=as.character(data[,1]))
			res_list <- lapply(data_list, function(x){
				x <- as.numeric(x)
				x <- x[!is.na(x)]
				if(length(x)>0){
					min(x)
				}else{
					NULL
				}
			})
			pval <- unlist(res_list)
		}
		
		# force those zeros to be miminum of non-zeros
		#tmp <- as.numeric(format(.Machine)['double.xmin'])
		tmp <- min(pval[pval!=0])
		pval[pval < tmp] <- tmp
	}
	
	## replace '_' with ':'
	tmp <- names(pval)
	tmp <- gsub("_", ":", tmp, perl=TRUE)
	## replace 'imm:' with 'chr'
	names(pval) <- gsub("imm:", "chr", tmp, perl=TRUE)
	
	Lead_Sig <- data.frame(SNP=names(pval), Sig=pval, row.names=NULL, stringsAsFactors=FALSE)
	leads <- Lead_Sig[,1]
	sigs <- Lead_Sig[,2]

	if(verbose){
		now <- Sys.time()
		message(sprintf("A total of %d Lead SNPs are input", length(leads)), appendLF=TRUE)
	}

	###########################
	## include additional SNPs that are in LD with input SNPs
	if(LD.r2>=0.5 & LD.r2<=1){
		default.include.LD <- c("ACB","AFR","AMR","ASW","BEB","CDX","CEU","CHB","CHS","CLM","EAS","ESN","EUR","FIN","GBR","GIH","GWD","IBS","ITU","JPT","KHV","LWK","MSL","MXL","PEL","PJL","PUR","SAS","STU","TSI","YRI")
		ind <- match(default.include.LD, include.LD)
		include.LD <- default.include.LD[!is.na(ind)]
	}else{
		include.LD <- NULL
	}
	
	LLR <- NULL
	if(length(include.LD) > 0 & is.null(LD.customised)){
	
		if(verbose){
			now <- Sys.time()
			message(sprintf("Inclusion of LD SNPs is based on population (%s) with R2 >= %f", paste(include.LD, collapse=','), LD.r2), appendLF=TRUE)
		}
	
		GWAS_LD <- xRDataLoader('GWAS_LD', RData.location=RData.location, guid=guid, verbose=verbose)
		res_list <- lapply(include.LD, function(x){
			ind <- match(x, names(GWAS_LD))
			data_ld <- GWAS_LD[[ind]]
			#data_ld <- ''
			#eval(parse(text=paste("data_ld <- GWAS_LD$", x, sep="")))
			ind <- match(rownames(data_ld), leads)
			ind_lead <- which(!is.na(ind))
			
			if(length(ind_lead) >= 2){
				ind_ld <- which(Matrix::colSums(data_ld[ind_lead,]>=LD.r2)>0)
				sLL <- data_ld[ind_lead, ind_ld]
				summ <- summary(sLL)
				res <- data.frame(Lead=rownames(sLL)[summ$i], LD=colnames(sLL)[summ$j], R2=summ$x, stringsAsFactors=FALSE)
			}else if(length(ind_lead) == 1){
				ind_ld <- which(data_ld[ind_lead,]>=LD.r2)
				sLL <- data_ld[ind_lead, ind_ld]
				res <- data.frame(Lead=rep(rownames(data_ld)[ind_lead],length(sLL)), LD=colnames(data_ld)[ind_ld], R2=sLL, stringsAsFactors=FALSE)
			}else{
				NULL
			}
		})
		## get data frame (Lead LD R2)
		LLR <- do.call(rbind, res_list)
		
		###########################
		## also based on ImmunoBase
		if(1){
			ImmunoBase_LD <- xRDataLoader('ImmunoBase_LD', RData.location=RData.location, guid=guid, verbose=verbose)
			res_list <- lapply(include.LD, function(x){
				ind <- match(x, names(ImmunoBase_LD))
				data_ld <- ImmunoBase_LD[[ind]]
				#data_ld <- ''
				#eval(parse(text=paste("data_ld <- ImmunoBase_LD$", x, sep="")))
				ind <- match(rownames(data_ld), leads)
				ind_lead <- which(!is.na(ind))
				
				if(length(ind_lead) >= 2){
					ind_ld <- which(Matrix::colSums(data_ld[ind_lead,]>=LD.r2)>0)
					sLL <- data_ld[ind_lead, ind_ld]
					summ <- summary(sLL)
					res <- data.frame(Lead=rownames(sLL)[summ$i], LD=colnames(sLL)[summ$j], R2=summ$x, stringsAsFactors=FALSE)
				}else if(length(ind_lead) == 1){
					ind_ld <- which(data_ld[ind_lead,]>=LD.r2)
					sLL <- data_ld[ind_lead, ind_ld]
					res <- data.frame(Lead=rep(rownames(data_ld)[ind_lead],length(sLL)), LD=colnames(data_ld)[ind_ld], R2=sLL, stringsAsFactors=FALSE)
				}else{
					NULL
				}
				
			})
			## get data frame (Lead LD R2)
			LLR_tmp <- do.call(rbind, res_list)
			LLR <- rbind(LLR, LLR_tmp)
		}
		
	###########################	
	}else if(!is.null(LD.customised)){
		if (is.vector(LD.customised)){
			# assume a file
			LLR <- utils::read.delim(file=LD.customised, header=FALSE, row.names=NULL, stringsAsFactors=FALSE)
		}else if(is.matrix(LD.customised) | is.data.frame(LD.customised)){
			LLR <- LD.customised
		}
		
		if(!is.null(LLR)){
			#############
			ind <- match(LLR[,1], leads)
			LLR <- LLR[!is.na(ind),]
			#############
			flag <- LLR[,3]>=LD.r2
			if(sum(flag)>0){
				LLR <- LLR[LLR[,3]>=LD.r2,]
				colnames(LLR) <- c("Lead", "LD", "R2")
			
				if(verbose){
					now <- Sys.time()
					message(sprintf("Inclusion of LD SNPs is based on customised data (%d Lead SNPs and %d LD SNPs) with R2>=%f", length(unique(LLR[,1])), length(unique(LLR[,2])), LD.r2), appendLF=TRUE)
				}
			}else{
				LLR <- NULL
			}
		}
		
	}
	
	if(!is.null(LLR)){
		## get data frame (LD Sig)
		ld_list <- split(x=LLR[,-2], f=LLR[,2])
		res_list <- lapply(ld_list, function(x){
			ind <- match(x$Lead, leads)
			## power transformation of p-values X R2, then keep the min (the most significant)
			min(sigs[ind] ^ x$R2)
		})
		vec <- unlist(res_list)
		LD_Sig <- data.frame(SNP=names(vec), Sig=vec, row.names=NULL, stringsAsFactors=FALSE)

		## merge Lead and LD
		df <- base::rbind(Lead_Sig, LD_Sig)
		res_list <- split(x=df$Sig, f=df$SNP)
		vec <- unlist(lapply(res_list, min))
		SNP_Sig <- data.frame(SNP=names(vec), FDR=vec, row.names=NULL, stringsAsFactors=FALSE)
	}else{
		if(verbose){
			now <- Sys.time()
			message(sprintf("Do not include any LD SNPs"), appendLF=TRUE)
		}
	
		SNP_Sig <- Lead_Sig
	}
	###########################
	
	if(verbose){
		now <- Sys.time()
		message(sprintf("A total of %d Lead/LD SNPs are considered", nrow(SNP_Sig)), appendLF=TRUE)
	}
	
	pval <- as.numeric(SNP_Sig[,2])
	names(pval) <- SNP_Sig[,1]
	
	# transformed into scores according to log-likelihood ratio between the true positives and the false positivies
    ## also take into account the given threshold of the significant level
    ## SNPs with p-value below this are considered significant and thus scored positively
    ## Instead, SNPs with p-values above this are considered insigificant and thus scored negatively (zero-out)
	
	if(is.null(significance.threshold)){
        scores <- log10((1-pval)/pval)
        #scores <- log10(1/pval)
    }else{
		scores <- log10((1-pval)/pval) - log10((1-significance.threshold)/significance.threshold)
		#scores <- log10(1/pval) - log10(1/significance.threshold)
	}
    ## replace those infinite values with the next finite ones
    tmp_max <- max(scores[!is.infinite(scores)])
    tmp_min <- min(scores[!is.infinite(scores)])
    scores[scores>tmp_max] <- tmp_max
    scores[scores<tmp_min] <- tmp_min
	## zero-out SNPs with negative scores
	ind_remained <- which(scores>0)
	seeds.snps <- scores[ind_remained]
	pval <- pval[ind_remained]
    
	if(verbose){
		now <- Sys.time()
		message(sprintf("A total of %d Lead/LD SNPs are scored positively", sum(seeds.snps>0)), appendLF=TRUE)
	}
    
    #########
    flag <- rep('Lead', length(pval))
    ind <- match(names(pval), Lead_Sig$SNP)
    flag[is.na(ind)] <- 'LD'
    
    df_SNP <- data.frame(SNP=names(pval), Score=seeds.snps, Pval=pval, Flag=flag, row.names=NULL, stringsAsFactors=FALSE)
    
    ##############################
    ## cap the maximum score
    if(!is.null(score.cap)){
    	score.cap <- as.numeric(score.cap)
    	if(score.cap <= max(df_SNP$Score)){
    		df_SNP$Score[df_SNP$Score>=score.cap] <- score.cap
    		
			if(verbose){
				now <- Sys.time()
				message(sprintf("SNP score capped to the maximum score %d.", score.cap), appendLF=TRUE)
			}
    	}
    }
    ##############################
        
    df_SNP <- df_SNP[order(df_SNP$Flag,df_SNP$Score,-df_SNP$Pval,df_SNP$SNP,decreasing=TRUE),]
    #########


    
    invisible(df_SNP)
}
