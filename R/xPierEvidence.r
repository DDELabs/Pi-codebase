#' Function to extract evidence from a list of pNode objects
#'
#' \code{xPierEvidence} is supposed to extract evidence from a list of pNode objects, in terms of seed genes under genetic influence.
#'
#' @param list_pNode a list of "pNode" objects or a "pNode" object
#' @param target.query which gene is in query. If NULL, all genes will be queried
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @return
#' a data frame of nPair X 5 containing Gene-SNP pair info per context, where the 6 columns are "Gene" (seed genes), "SNP" (dbSNP), "Score" (an SNP's genetic influential score on a seed gene), "Context" (predictors), "Flag" (indicative of Lead SNPs or LD SNPs), and "Pval" (the SNP p-value)
#' @note none
#' @export
#' @seealso \code{\link{xPierEvidence}}
#' @include xPierEvidence.r
#' @examples
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' \dontrun{
#' df_Gene2SNP <- xPierEvidence(ls_pNode)
#' }

xPierEvidence <- function(list_pNode, target.query=NULL, verbose=TRUE)
{
    
   	if(is(list_pNode,"pNode")){
		if(is.null(list_pNode$Gene2SNP)){
			return(NULL)
		}else{
			list_pNode <- list(list_pNode)
		}
	}else if(is(list_pNode,"list")){
		######
		## check pNode has a componenet called 'Gene2SNP'
		list_pNode <- lapply(list_pNode, function(pNode){
			if(is.null(pNode$Gene2SNP)){
				return(NULL)
			}else{
				pNode
			}
		})
		######
		## Remove null elements in a list
		list_pNode <- base::Filter(base::Negate(is.null), list_pNode)
		if(length(list_pNode)==0){
			return(NULL)
		}
	}else{
		stop("The function must apply to 'list' of 'pNode' objects or a 'pNode' object.\n")
	}
	
	## Combine into a data frame called 'df_predictor'
	list_names <- names(list_pNode)
	if(is.null(list_names)){
		list_names <- paste('Predictor', 1:length(list_pNode), sep=' ')
		names(list_pNode) <- list_names
	}
	ls_Gene2SNP <- lapply(1:length(list_pNode), function(i){
		pNode <- list_pNode[[i]]
		res <- cbind(pNode$Gene2SNP, Context=rep(list_names[i],nrow(pNode$Gene2SNP)))
	})
	df_Gene2SNP <- do.call(rbind, ls_Gene2SNP)
	
	df_SNP <- list_pNode[[1]]$SNP[,c("SNP","Flag","Pval")]
	ind <- match(df_Gene2SNP$SNP, df_SNP$SNP)
	df_Gene2SNP$Flag <- df_SNP$Flag[ind]
	df_Gene2SNP$Pval <- df_SNP$Pval[ind]
	
	if(!is.null(target.query)){
		ind <- match(df_Gene2SNP$Gene, target.query)
		ind <- which(!is.na(ind))
		if(length(ind)!=0){
			df_Gene2SNP <- df_Gene2SNP[ind,]
		}else{
			df_Gene2SNP <- NULL
		}
	}
	
    invisible(df_Gene2SNP)
}
