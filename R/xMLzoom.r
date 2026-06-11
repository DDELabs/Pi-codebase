#' Function to visualise machine learning results using zoom plot
#'
#' \code{xMLzoom} is supposed to visualise machine learning results using zoom plot. It returns an object of class "ggplot".
#'
#' @param xTarget an object of class "sTarget" or "dTarget" (with the component 'pPerf')
#' @param top the number of the top targets to be labelled/highlighted
#' @param top.label.type how to label the top targets. It can be "box" drawing a box around the labels , and "text" for the text only
#' @param top.label.size the highlight label size
#' @param top.label.query which top genes in query will be labelled. By default, it sets to NULL meaning all top genes will be displayed. If labels in query can not be found, then none will be displayed
#' @param point.shape an integer specifying point shapes. By default, it is 3 for cross. For details, please refere to \url{http://sape.inf.usi.ch/quick-reference/ggplot2/shape}
#' @param font.family the font family for texts
#' @param signature logical to indicate whether the signature is assigned to the plot caption. By default, it sets TRUE
#' @return an object of class "ggplot"
#' @note none
#' @export
#' @seealso \code{\link{xColormap}}
#' @include xMLzoom.r
#' @examples
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' \dontrun{
#' gp <- xMLzoom(sTarget)
#' gp
#' }

xMLzoom <- function(xTarget, top=20, top.label.type=c("box","text"), top.label.size=3, top.label.query=NULL, point.shape=3, font.family="sans", signature=TRUE)
{

    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    top.label.type <- match.arg(top.label.type)

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
	df <- data.frame(Symbol=rownames(priority), GS=priority$GS, Score=priority$rating, stringsAsFactors=FALSE)
    
    #GS_level<-c("GSN","GSP","NEW")
    GS_level <- sort(unique(df$GS))
    df$GS <- factor(df$GS, levels=GS_level)
    color <- xColormap("ggplot2")(length(GS_level))
    
	df_highlight <- df[1:top, ]
	top_cutoff <- base::min(df_highlight[,3])
    
    GS <- Score <- Symbol <- ''
    
	gp <- ggplot(df, aes(x=Score, y=GS, color=GS)) + geom_point(alpha=0.95,shape=point.shape)
	gp <- gp + scale_color_manual(values=color)
	gp <- gp + theme_bw() + theme(legend.position="none", legend.title=element_blank(), axis.title.y=element_blank(), axis.text.y=element_text(size=12,face='bold'), axis.title.x=element_text(size=14,color="black",face="bold"))
	gp <- gp + xlab("Pi rating\n(quantifying separation between GSP and GSN)")
    
    ## zoom
    gp <- gp + ggforce::facet_zoom(x=(Score>=top_cutoff), zoom.data=(Score>=top_cutoff), zoom.size=1.5, show.area=TRUE)
    
	########################################
	## ONLY restricted to genes in query
	if(!is.null(top.label.query)){
		top.label.query <- as.vector(t(top.label.query)) # just in case converting data.frame to vector
		ind <- match(df_highlight$Symbol, top.label.query)
		if(sum(!is.na(ind)) >= 1){
			df_highlight <- df_highlight[!is.na(ind),]
		}else{
			message(sprintf("Note: none of genes in query found at the top %d list", top), appendLF=TRUE)
			df_highlight <- NULL
		}
	}
	########################################
	if(!is.null(df_highlight)){
		if(top.label.type=="text"){
			gp <- gp + ggrepel::geom_text_repel(data=df_highlight, aes(x=Score,y=GS,label=Symbol), size=top.label.size, fontface='bold', point.padding=unit(0.2,"lines"), segment.color='grey50', segment.alpha=0.5, arrow=arrow(length=unit(0.01,'npc')))
		}else if(top.label.type=="box"){
			gp <- gp + ggrepel::geom_label_repel(data=df_highlight, aes(x=Score,y=GS,label=Symbol), size=top.label.size, fontface='bold', box.padding=unit(0.35,"lines"), point.padding=unit(0.35,"lines"), segment.color='grey50', segment.alpha=0.5, arrow=arrow(length=unit(0.01,'npc')))
		}
	}
	
	## caption
    if(signature){
    	caption <- paste("Created by xMLzoom from Pi version", utils ::packageVersion("Pi"))
    	gp <- gp + labs(caption=caption) + theme(plot.caption=element_text(hjust=1,face='bold.italic',size=8,colour='#002147'))
    }
    
	## change font family to 'Arial'
	gp <- gp + theme(text=element_text(family=font.family))
	
	## put arrows on x-axis
	gp <- gp + theme(axis.line.x=element_line(arrow=arrow(angle=30,length=unit(0.25,"cm"), type="open")))
    
    invisible(gp)
}


