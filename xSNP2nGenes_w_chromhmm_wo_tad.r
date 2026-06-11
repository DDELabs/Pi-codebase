xSNP2nGenes <- function(data, distance.max=200000, decay.kernel=c("rapid","slow","linear","constant"), decay.exponent=2, GR.SNP=c("dbSNP_GWAS","dbSNP_Common","dbSNP_Single"), GR.Gene=c("UCSC_knownGene","UCSC_knownCanonical"), include.TAD=c("none","GM12878","IMR90","MSC","TRO","H1","MES","NPC"), verbose=TRUE, RData.location="http://galahad.well.ox.ac.uk/bigdata", guid=NULL)
{
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    decay.kernel <- match.arg(decay.kernel)
    # include.TAD <- match.arg(include.TAD)
	
	## replace '_' with ':'
	data <- gsub("_", ":", data, perl=TRUE)
	## replace 'imm:' with 'chr'
	data <- gsub("imm:", "chr", data, perl=TRUE)
	
	data <- unique(data)
	
    ######################################################
    # Link to targets based on genomic distance
    ######################################################
    
    gr_SNP <- xSNPlocations(data=data, GR.SNP=GR.SNP, verbose=verbose, RData.location=RData.location, guid=guid)

  	#######################################################
  	
	if(verbose){
		now <- Sys.time()
		message(sprintf("Load positional information for Genes (%s) ...", as.character(now)), appendLF=TRUE)
	}
	if(is(GR.Gene,"GRanges")){
			gr_Gene <- GR.Gene
	}else{
		gr_Gene <- xRDataLoader(GR.Gene[1], verbose=verbose, RData.location=RData.location, guid=guid)
		if(is.null(gr_Gene)){
			GR.Gene <- "UCSC_knownGene"
			if(verbose){
				message(sprintf("Instead, %s will be used", GR.Gene), appendLF=TRUE)
			}
			gr_Gene <- xRDataLoader(GR.Gene, verbose=verbose, RData.location=RData.location, guid=guid)
		}
    }
    
	if(verbose){
		now <- Sys.time()
		message(sprintf("Define nearby genes (%s) ...", as.character(now)), appendLF=TRUE)
	}
    
	# genes: get all UCSC genes within defined distance window away from variants
	maxgap <- distance.max-1
	#minoverlap <- 1L # 1b overlaps
	minoverlap <- 0L
	subject <- gr_Gene
	query <- gr_SNP
	q2r <- as.matrix(as.data.frame(suppressWarnings(GenomicRanges::findOverlaps(query=query, subject=subject, maxgap=maxgap, minoverlap=minoverlap, type="any", select="all", ignore.strand=TRUE))))
	
	if(length(q2r) > 0){
	
		if(verbose){
			now <- Sys.time()
			message(sprintf("Calculate distance (%s) ...", as.character(now)), appendLF=TRUE)
		}
	
		if(1){
			### very quick
			x <- subject[q2r[,2],]
			y <- query[q2r[,1],]
			
			chrom_overlaps <- findOverlaps(y, n_chrom, type = "within")
			reg_values <- n_chrom$reg[subjectHits(chrom_overlaps)]
			
			dists <- GenomicRanges::distance(x, y, select="all", ignore.strand=TRUE)
			
			###
			df_y <- GenomicRanges::as.data.frame(y, row.names=NULL)
			df_x <- GenomicRanges::as.data.frame(x, row.names=NULL)
			
			if(0){
				## Gap defined as the regions from an SNP to the closest end of a Gene
				df_interval <- data.frame(seqnames=df_y$seqnames, start=df_y$start, end=df_y$end, stringsAsFactors=FALSE)
				ind <- df_y$start < df_x$start
				df_interval[ind,] <- data.frame(seqnames=df_y$seqnames[ind], start=df_y$start[ind], end=df_x$start[ind], stringsAsFactors=FALSE)
				ind <- df_y$start > df_x$end
				df_interval[ind,] <- data.frame(seqnames=df_y$seqnames[ind], start=df_x$end[ind], end=df_y$start[ind], stringsAsFactors=FALSE)
				
			}else{
				## Gap defined as the regions from an SNP to the TSS of a Gene
				
				df_interval <- data.frame(seqnames=df_y$seqnames, start=df_y$start, end=df_y$end, stringsAsFactors=FALSE)
				ind <- df_x$strand=='+'
				df_interval[ind,] <- data.frame(seqnames=df_y$seqnames[ind], start=df_y$start[ind], end=df_x$start[ind], stringsAsFactors=FALSE)
				ind <- df_x$strand=='-'
				df_interval[ind,] <- data.frame(seqnames=df_y$seqnames[ind], start=df_y$start[ind], end=df_x$end[ind], stringsAsFactors=FALSE)
				ind <- df_x$strand=='*'
				df_interval[ind,] <- data.frame(seqnames=df_y$seqnames[ind], start=df_y$start[ind], end=(df_x$start[ind]+df_x$end[ind])/2, stringsAsFactors=FALSE)
				
				## swap the location of start and end
				ind <- df_interval$start > df_interval$end
				df_interval[ind,] <- data.frame(seqnames=df_interval$seqnames[ind], start=df_interval$end[ind], end=df_interval$start[ind], stringsAsFactors=FALSE)
			}
			
			###########################################################
			########## very important!
			## if SNP within a gene body no restriction will apply (that is, SNP location)
			df_interval[dists==0, 'start'] <-  df_y$start[dists==0]
			df_interval[dists==0, 'end'] <-  df_y$end[dists==0]
			###########################################################			
			
			vec_interval <- paste0(df_interval$seqnames, ':', as.character(df_interval$start), '-', as.character(df_interval$end))
			###
			
			df_nGenes <- data.frame(Gene=names(x), SNP=names(y), Dist=dists, r=reg_values, Gap=vec_interval, stringsAsFactors=FALSE)
		}
	
		## weights according to distance away from SNPs
		if(distance.max==0){
			x <- df_nGenes$Dist
		}else{
			x <- df_nGenes$Dist / distance.max
		}
		if(decay.kernel == 'slow'){
			# y <- 1-(x)^decay.exponent
		  sigmoid <- function(x) {
		    return(1 / (1 + exp(-x)))
		  }
			y <- (1-(x)^decay.exponent) * sigmoid (10 * df_nGenes$r)
		}else if(decay.kernel == 'rapid'){
			# y <- (1-x)^decay.exponent
		  sigmoid <- function(x) {
		    return(1 / (1 + exp(-x)))
		  }
		  y <- ((1-x)^decay.exponent) * sigmoid (10 * df_nGenes$r)
		}else if(decay.kernel == 'linear'){
			# y <- 1-x
			sigmoid <- function(x) {
			  return(1 / (1 + exp(-x)))
			}
			y <- (1-x) * sigmoid (10 * df_nGenes$r)
		}else{
			# y <- 1
		  sigmoid <- function(x) {
		    return(1 / (1 + exp(-x)))
		  }
		  y <- (1) * sigmoid (10 * df_nGenes$r)
		}
		df_nGenes$Weight <- y
	
		if(verbose){
			now <- Sys.time()
			message(sprintf("%d Genes are defined as nearby genes within %d(bp) genomic distance window using '%s' decay kernel (%s)", length(unique(df_nGenes$Gene)), distance.max, decay.kernel, as.character(now)), appendLF=TRUE)
		}
		
		df_nGenes <- df_nGenes[,c('Gene','SNP','Dist','Weight','Gap')]
		df_nGenes <- df_nGenes[order(df_nGenes$Gene,df_nGenes$Dist,decreasing=FALSE),]
		
	}else{
		df_nGenes <- NULL
		
		if(verbose){
			now <- Sys.time()
			message(sprintf("No nearby genes are defined"), appendLF=TRUE)
		}
	}
	
	###########################
	## include TAD
	default.include.TAD <- c("GM12878","IMR90","MSC","TRO","H1","MES","MES")
	ind <- match(default.include.TAD, include.TAD)
	include.TAD <- default.include.TAD[!is.na(ind)]
	if(length(include.TAD) > 0){
		if(verbose){
			now <- Sys.time()
			message(sprintf("Inclusion of TAD boundary regions is based on '%s'", include.TAD), appendLF=TRUE)
		}
		df_nGenes$TAD <- rep('Excluded', nrow(df_nGenes))
		
		TAD <- xRDataLoader(paste0('TAD.',include.TAD), RData.location=RData.location, guid=guid, verbose=verbose)
		# TAD <- include.TAD
		iGR <- xGR(data=df_nGenes$Gap, format="chr:start-end", RData.location=RData.location, guid=guid)
		q2r <- as.matrix(as.data.frame(suppressWarnings(GenomicRanges::findOverlaps(query=iGR, subject=TAD, type="within", select="all", ignore.strand=TRUE))))
		q2r <- q2r[!duplicated(q2r[,1]), ]
		df_nGenes$TAD[q2r[,1]] <- GenomicRanges::mcols(TAD)[q2r[,2],]
		
		###########################################################
		########## very important!
		## if SNP within a gene body no restriction will apply (that is, SNP location)
		df_nGenes$TAD[df_nGenes$Dist==0] <-  df_nGenes$Gap[df_nGenes$Dist==0]
		###########################################################	
		
		if(verbose){
			now <- Sys.time()
			message(sprintf("\t%d out of %d SNP-nGene pairs are within the same TAD boundary regions", sum(df_nGenes$TAD!='Excluded'), length(iGR)), appendLF=TRUE)
			message(sprintf("\t%d out of %d genes are defined as nearby genes after considering TAD boundary regions", length(unique(df_nGenes[df_nGenes$TAD!='Excluded','Gene'])), length(unique(df_nGenes$Gene))), appendLF=TRUE)
		}
		
	}
	
	####################################
	# only keep those genes with GeneID
	####################################
	if(!is.null(df_nGenes)){
		ind <- xSymbol2GeneID(df_nGenes$Gene, details=FALSE, verbose=verbose, RData.location=RData.location, guid=guid)
		df_nGenes <- df_nGenes[!is.na(ind), ]
		if(nrow(df_nGenes)==0){
			df_nGenes <- NULL
		}
	}
	####################################
	
    invisible(df_nGenes)
}