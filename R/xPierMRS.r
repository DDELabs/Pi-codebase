#' Function to calculate multi-trait rating score from a list of dTarget/sTarget objects
#'
#' \code{xPierMRS} is supposed to calculate multi-trait rating score (MRS) from a list of dTarget/sTarget objects.
#'
#' @param list_xTarget a list of "dTarget"/"sTarget" objects
#' @param cutoff.rank the rank cutoff. By default it is 150
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @return
#' a data frame containing columns 'Target', 'MRS', 'rating' (in the form of "rating.trait_names") and 'rank' (in the form of "rank.trait_names").
#' @note none
#' @export
#' @seealso \code{\link{xPierCross}}
#' @include xPierMRS.r
#' @examples
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' \dontrun{
#' df_MRS <- xPierMRS(ls_xTarget)
#' }

xPierMRS <- function(list_xTarget, cutoff.rank=150, verbose=TRUE)
{
    
   	if(is(list_xTarget,"list")){
		## Remove null elements in a list
		list_xTarget <- base::Filter(base::Negate(is.null), list_xTarget)
		if(length(list_xTarget)==0){
			return(NULL)
		}
	}else{
		stop("The function must apply to 'list' of 'dTarget' objects.\n")
	}
	
	## df_rank
	df_rank <- xPierCross(list_xTarget, displayBy="rank", combineBy='union', aggregateBy="none")
	## df_rating
	df_rating <- xPierCross(list_xTarget, displayBy="rating", combineBy='union', aggregateBy="none")
	
	## MRS
	ind <- which(apply(df_rank<=cutoff.rank, 1, sum)>0)
	df <- df_rank[ind,]
	df[df>cutoff.rank] <- NA
	ntop <- apply(df, 1, function(x) sum(!is.na(x)))
	mrank <- apply(df, 1, function(x) median(x,na.rm=TRUE))
	n <- ncol(df)
	mrs <- (ntop+1-(mrank-1)/cutoff.rank)/(n+1)
	mrs <- sort(mrs, decreasing=TRUE)
	
	## df_output
	ind <- match(names(mrs), rownames(df_rating))
	df_MRS <- data.frame(Target=names(mrs), MRS=mrs, rating=df_rating[ind,], rank=df_rank[ind,], stringsAsFactors=FALSE)

	if(verbose){
		message(sprintf("A total of %d genes used for MRS calculation (based on rank cutoff %d and %d traits", nrow(df_MRS), cutoff.rank, ncol(df_rank)), appendLF=TRUE)
		
	}
	
    invisible(df_MRS)
}
