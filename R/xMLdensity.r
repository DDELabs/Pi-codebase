#' Function to visualise machine learning results using density plot
#'
#' \code{xMLdensity} is supposed to visualise machine learning results using density plot. It returns an object of class "ggplot".
#'
#' @param xTarget an object of class "xTarget" or "dTarget" (with the component 'pPerf')
#' @param displayBy which targets will be used for displaying. It can be one of "GS" for gold standard targets, "GSN" for gold standard negatives, "GSP" for gold standard positives, "NEW" for putative/new targets (non-GS), "All" for all targets (by default)
#' @param x.scale how to transform the x scale. It can be "normal" for no transformation, and "sqrt" for square root transformation (by default)
#' @param font.family the font family for texts
#' @param signature logical to indicate whether the signature is assigned to the plot caption. By default, it sets TRUE showing which function is used to draw this graph
#' @return an object of class "ggplot"
#' @note none
#' @export
#' @seealso \code{\link{xColormap}}
#' @include xMLdensity.r
#' @examples
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' \dontrun{
#' gp <- xMLdensity(xTarget, displayBy="All")
#' gp
#' }

xMLdensity <- function(xTarget, displayBy=c("All","GS","GSN","GSP","NEW"), x.scale=c("sqrt","normal"), font.family="sans", signature=TRUE) 
{
    
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    displayBy <- match.arg(displayBy)
    x.scale <- match.arg(x.scale)
    
    if(is(xTarget,"dTarget")){
    	if(is.null(xTarget$pPerf)){
    		stop("The function must apply to a 'dTarget' object with the component 'pPerf'.\n")
    	}
    }else{
		if(!is(xTarget,"sTarget")){
			stop("The function must apply to a 'sTarget' object.\n")
		}
    }

	priority <- xTarget$priority
	df <- data.frame(GS=priority$GS, Score=priority$rating, stringsAsFactors=FALSE)
    
    #GS_level<-c("GSN","GSP","NEW")
    GS_level <- sort(unique(df$GS))
    df$GS <- factor(df$GS, levels=GS_level)
    color <- xColormap("ggplot2")(length(GS_level))
    
    if(displayBy == "GS"){
    	df <- df[df$GS!='NEW',]
    	color <- color[1:2]
    }else if(displayBy=='GSN'){
    	df <- df[df$GS=='GSN',]
    	color <- color[1]
    }else if(displayBy=='GSP'){
    	df <- df[df$GS=='GSP',]
    	color <- color[2]
    }else if(displayBy=='NEW'){
    	df <- df[df$GS=='NEW',]
    	color <- color[3]
    }
    
    GS <- Score <- ''
    
	gp <- ggplot(df, aes(Score, fill=GS, color=GS)) + geom_density(alpha=0.1,adjust=1)
	gp <- gp + scale_color_manual(values=color) + scale_fill_manual(values=color)
	gp <- gp + theme_bw() + theme(legend.position="right", legend.title=element_blank(), axis.title.y=element_text(size=14,color="black",face="bold"), axis.title.x=element_text(size=14,color="black",face="bold"), panel.background=element_rect(fill=rgb(0.95,0.95,0.95,1)))
	gp <- gp + xlab("Composite scores\n(5-star ratings quantifying separation between GSP and GSN)") + ylab("Density of target genes")
	
	## x scale
    if(x.scale=="sqrt"){
    	x <- NULL
    	gp <- gp + scale_x_continuous(trans=scales::sqrt_trans(), breaks=scales::trans_breaks("log10",function(x) 10^x, n=2), limits=c(0,5))
    }else{
    	gp <- gp + xlim(0,5)
    }

	## caption
    if(signature){
    	caption <- paste("Created by xMLdensity from Pi version", utils ::packageVersion("Pi"))
    	gp <- gp + labs(caption=caption) + theme(plot.caption=element_text(hjust=1,face='bold.italic',size=8,colour='#002147'))
    }

	## change font family to 'Arial'
	gp <- gp + theme(text=element_text(family=font.family))

	## put arrows on x- and y-axis
	gp <- gp + theme(axis.line.x=element_line(arrow=arrow(angle=30,length=unit(0.25,"cm"), type="open")), axis.line.y=element_line(arrow=arrow(angle=30,length=unit(0.25,"cm"), type="open")))
	
	invisible(gp)
}
