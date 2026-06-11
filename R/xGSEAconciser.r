#' Function to make GSEA results conciser by removing redundant terms
#'
#' \code{xGSEAconciser} is supposed to make GSEA results conciser by removing redundant terms. A redundant term (called 'B') is claimed if its overlapped part (A&B) with a more significant term (called 'A') meets both criteria: 1) |A&B| > 0.9*|B|; and 2) |A&B| > 0.5*|A|.
#'
#' @param eGSEA an object of class "eGSEA"
#' @param cutoff a cutoff vector used to remove redundant terms. By default, it has the first element 0.9 and the second element 0.5. It means, for a term (less significant; called 'B'), if there is a more significant term (called 'A'), their overlapped members cover at least 90% of the B's members, and also at least 50% of the A's members, then this term B will be defined as redundant and thus being removed
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to false for no display
#' @return
#' an object of class "eGSEA", after redundant terms being removed.
#' @note none
#' @export
#' @seealso \code{\link{xGSEAconciser}}
#' @include xGSEAconciser.r
#' @examples
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' \dontrun{
#' eGSEA_concise <- xGSEAconciser(eGSEA)
#' }

xGSEAconciser <- function(eGSEA, cutoff=c(0.9,0.5), verbose=TRUE) 
{
        
    if(is(eGSEA,"eGSEA")){
		
		cross <- eGSEA$cross
		
		df <- eGSEA$df_summary
		ind <- match(df$setID, rownames(cross))
		cross <- cross[ind, ind]
		
		nRedundant_1 <- matrix(0, nrow=ncol(cross), ncol=ncol(cross))
		nRedundant_2 <- matrix(0, nrow=ncol(cross), ncol=ncol(cross))
		for(j in seq(2, ncol(cross))){
			for(i in seq(1, j-1)){
				nRedundant_1[i,j] <- cross[i,j] >= cross[j,j]*cutoff[1]
				nRedundant_2[i,j] <- cross[i,j] >= cross[i,i]*cutoff[2]
			}
		}
		nRedundant <- apply(nRedundant_1 & nRedundant_2, 2, sum)
		names(nRedundant) <- colnames(cross)
		
		ind <- match(df$setID, names(nRedundant))
		nRedundant <- nRedundant[ind]
		
		## update eGSEA by removing redundant terms
		flag <- nRedundant == 0
		
        if(verbose){
            now <- Sys.time()
            message(sprintf("\tAmong %d terms, there are %d non-redundant terms", length(flag), sum(flag)), appendLF=TRUE)
        }
		
		eGSEA$df_summary <- eGSEA$df_summary[flag, ]
		eGSEA$leading <- eGSEA$leading[flag]
		eGSEA$full <- eGSEA$full[flag]
		eGSEA$cross <- eGSEA$cross[flag, flag]
		
		res <- eGSEA
		
	}else{
		res <- NULL
		message("The provided object can be used for removing redundant terms.\n")
	}
    
    res
}
