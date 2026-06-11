#' Function to visualise ABF evidence using heatmap
#'
#' \code{xPierABFheatmap} is supposed to visualise ABF evidence using heatmap. It returns an object of class "ggplot".
#'
#' @param data an input vector containing gene symbols
#' @param xTarget a "dTarget" or "sTarget" object with the componet 'list_pNode' related to 'eGene' predictors. Alternatively, it can be a data frame with columns ('context','mode','probeID','Symbol','gene_cse','snps','snp_cse','A1','A2','b_GWAS','b_eQTL','b_ABF','p_GWAS','p_eQTL','pp_ABF','H4')
#' @param type the type of the heatmap. It can be "Gene" (gene-centric heatmap) or "Gene_SNP" (heatmap for the gene-snp pair)
#' @param colormap short name for the colormap. It can be one of "jet" (jet colormap), "bwr" (blue-white-red colormap), "gbr" (green-black-red colormap), "wyr" (white-yellow-red colormap), "br" (black-red colormap), "yr" (yellow-red colormap), "wb" (white-black colormap), and "rainbow" (rainbow colormap, that is, red-yellow-green-cyan-blue-magenta). Alternatively, any hyphen-separated HTML color names, e.g. "blue-black-yellow", "royalblue-white-sandybrown", "darkgreen-white-darkviolet". A list of standard color names can be found in \url{http://html-color-codes.info/color-names}
#' @param zlim the minimum and maximum z values for which colors should be plotted
#' @return an object of class "ggplot" appended with 'mat' (the matrix colored by 'b_ABF') and 'df' (a data frame with columns 'priority','code','context','mode','ProbeID','Symbol','gene_cse','snps','snp_cse','A1','A2','b_GWAS','b_eQTL','b_ABF','p_GWAS','p_eQTL','pp_ABF','direction_GWAS','direction_eQTL','direction_ABF').
#' @note none
#' @export
#' @seealso \code{\link{xSparseMatrix}}, \code{\link{xHeatmap}}, \code{\link{xSparseMatrix}}
#' @include xPierABFheatmap.r
#' @examples
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' \dontrun{
#' gp <- xPierABFheatmap(data, dTarget)
#' }

xPierABFheatmap <- function(data, xTarget, type=c('Gene','Gene_SNP'), colormap='steelblue-lightyellow-orange', zlim=c(-0.5,0.5))
{
    type <- match.arg(type)
    
    ## df_evidence
    if(is(xTarget,"sTarget") | is(xTarget,"dTarget")){
    
    	if(any(names(xTarget) %in% 'list_pNode')){
			list_pNode <- xTarget$list_pNode
			ind <- which(grepl('eGene_', names(list_pNode)))
			ls_df <- lapply(list_pNode[ind], function(x){
				x$evidence
			})
			df_evidence <- do.call(rbind, ls_df)
		
			## df_gene_priority: extract priority
			df_gene_priority <- data.frame(gene=xTarget$priority$name, priority=xTarget$priority$rating, stringsAsFactors=FALSE)
			ind <- match(df_gene_priority$gene, data)
			if(sum(!is.na(ind))>0){
				df_gene_priority <- df_gene_priority[!is.na(ind),]
			}
		
		}else{
			return(NULL)
		}
		
    }else if(is(xTarget,"data.frame")){
    
    	df_evidence <- xTarget[,c('context','mode','probeID','Symbol','gene_cse','snps','snp_cse','A1','A2','b_GWAS','b_eQTL','b_ABF','p_GWAS','p_eQTL','pp_ABF','H4')]
    	
		## df_gene_priority: always priority=1
		df_gene_priority <- data.frame(gene=data, priority=1, stringsAsFactors=FALSE)
    	
    }else{
    	return(NULL)
    	
    }
    
    ######################################################################################
    if(is.null(df_evidence)){
    	return(NULL)
    }
    ######################################################################################
    
	## df_evidence_priority
	ind <- match(df_evidence$Symbol, df_gene_priority$gene)
	df_evidence_priority <- df_evidence[!is.na(ind),]
	df_evidence_priority$priority <- df_gene_priority$priority[ind[!is.na(ind)]]

	if(nrow(df_evidence_priority)==0){
		return(NULL)
	}

	#######################
	# different effect alleles for the same SNP
	#######################
	## df_tmp
	df_tmp <- df_evidence_priority
	df_tmp$uid <- paste0(df_tmp$snps, '|', df_tmp$A1, '|', df_tmp$A2)
	## df_ttmp	
	df_ttmp <- unique(df_tmp[,c('snps','A1','A2')])
	df_ttmp$uid <- paste0(df_ttmp$snps, '|', df_ttmp$A1, '|', df_ttmp$A2)
	## vec
	ind <- which(duplicated(df_ttmp[,'snps']))
	vec <- df_ttmp[ind,'snps']
	names(vec) <- df_ttmp[ind,'uid']
	## wrap
	for(i in 1:length(vec)){
		x <- vec[i]
		ind <- match(df_tmp$uid, names(x))
		j <- which(!is.na(ind))
		tmp <- df_tmp[j,'A1']
		df_tmp[j,'A1'] <- df_tmp[j,'A2']
		df_tmp[j,'A2'] <- tmp
		df_tmp[j,'b_GWAS'] <- -1 * df_tmp[j,'b_GWAS']
		df_tmp[j,'b_eQTL'] <- -1 * df_tmp[j,'b_eQTL']
	}
	df_evidence_priority <- df_tmp[,-ncol(df_tmp)]
	
	######################################################################################
	
	context <- code <- code1 <- priority <- pp_ABF <- NULL
	
	## append 'direction_GWAS', 'direction_eQTL', 'direction_ABF' and 'code'
	df <- df_evidence_priority
	df$direction_GWAS <- ifelse(df$b_GWAS>0, 'RSK', 'PRT')
	df$direction_eQTL <- ifelse(df$b_eQTL>0, 'ASC', 'DES')
	df$direction_ABF <- ifelse(df$b_ABF>0, 'POS', 'NEG')
	df$code <- paste0(df$Symbol,' (',df$snps,') [',df$A1,'<',df$A2,'; ',df$direction_GWAS,'; ',df$direction_eQTL,' - ', df$mode,  ']')
	
	## keep the best pp_ABF per context and code1
	code1 <- NULL
	df$code1 <- paste0(df$Symbol,' (',df$snps,') [',df$A1,'<',df$A2,']')
	df <- df %>% dplyr::arrange(context,code1,-pp_ABF)
	ind <- which(!duplicated(df[,c("context","code1")]))
	df <- df[ind,]	
	
	## order by priority
	df_output <- df %>% dplyr::arrange(-priority)
	
	if(type=='Gene_SNP'){
		##################
		df <- df %>% dplyr::arrange(-priority,code)
		##################
	
		## columns alwawys sorted by pre-defined ones
		## rows always sorted by priority
		default.Context <- c("CD14","LPS2","LPS24","IFN","Bcell","NK","Neutrophil","CD4","CD8","Blood","Monocyte","shared_CD14","shared_LPS2","shared_LPS24","shared_IFN")
		ind <- match(default.Context, unique(df$context))
		columns <- default.Context[!is.na(ind)]
		mat <- as.matrix(xSparseMatrix(df[,c('code','context','b_ABF')], rows=unique(df$code), columns=columns, verbose=FALSE))
		mat[mat==0] <- NA
		gp <- xHeatmap(mat, reorder="none", colormap=colormap, ncolors=64, zlim=zlim, legend.title="Effect", barwidth=0.4, x.rotate=60, shape=19, size=2, x.text.size=7, y.text.size=6, legend.text.size=5, legend.title.size=7, na.color='transparent', barheight=max(3,min(5,nrow(mat))))
		## append 'mat'
		gp$mat <- mat
	
	}else if(type=='Gene'){
		##################
		## only the one with the highest pp_ABF if multiple found in the same context given a gene
		Symbol <- NULL
		df <- df %>% dplyr::arrange(context,Symbol,-pp_ABF)
		ind <- which(!duplicated(df[,c("context","Symbol")]))
		df <- df[ind,]
		df <- df %>% dplyr::arrange(-priority)
		##################
		
		## columns alwawys sorted by pre-defined ones
		## rows always sorted by priority
		default.Context <- c("CD14","LPS2","LPS24","IFN","Bcell","NK","Neutrophil","CD4","CD8","Blood","Monocyte","shared_CD14","shared_LPS2","shared_LPS24","shared_IFN")
		ind <- match(default.Context, unique(df$context))
		columns <- default.Context[!is.na(ind)]
		mat <- as.matrix(xSparseMatrix(df[,c('Symbol','context','b_ABF')], rows=unique(df$Symbol), columns=columns, verbose=FALSE))
		mat[mat==0] <- NA
		gp <- xHeatmap(mat, reorder="none", colormap=colormap, ncolors=64, zlim=zlim, legend.title="Effect", barwidth=0.4, x.rotate=60, shape=19, size=2, x.text.size=7, y.text.size=6, legend.text.size=5, legend.title.size=7, na.color='transparent', barheight=max(3,min(5,nrow(mat))))
		## append 'mat'
		gp$mat <- mat
		
	}
	
	###########################
	## append 'df'
	gp$df <- df_output[,c('priority','code','context','mode','ProbeID','Symbol','gene_cse','snps','snp_cse','A1','A2','b_GWAS','b_eQTL','b_ABF','p_GWAS','p_eQTL','pp_ABF','direction_GWAS','direction_eQTL','direction_ABF')]
	###########################
		
	invisible(gp)
}

