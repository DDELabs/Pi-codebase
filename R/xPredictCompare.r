#' Function to compare prediction performance results
#'
#' \code{xPredictCompare} is supposed to compare prediction performance results. It returns an object of class "ggplot".
#'
#' @param list_pPerf a list of "pPerf" objects
#' @param displayBy which performance will be used for comparison. It can be "ROC" for ROC curve (by default), "PR" for PR curve
#' @param type the type of plot to draw. It can be "curve" for curve plot (by default), "bar" for bar plot
#' @param sort logical to indicate whether to sort methods according to performance. By default, it sets TRUE
#' @param detail logical to indicate whether to label methods along with performance. By default, it sets TRUE
#' @param facet logical to indicate whether to facet/wrap a 1d of panels into 2d. By default, it sets FALSE
#' @param font.family the font family for texts
#' @param signature a logical to indicate whether the signature is assigned to the plot caption. By default, it sets TRUE showing which function is used to draw this graph
#' @return an object of class "ggplot" or NULL (if all input pPerf objects are NULL)
#' @note none
#' @export
#' @seealso \code{\link{xPredictCompare}}
#' @include xPredictCompare.r
#' @examples
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' \dontrun{
#' bp <- xPredictCompare(ls_pPerf, displayBy="ROC")
#' print(bp)
#' ## modify legend position
#' bp + theme(legend.position=c(0.75,0.25))
#' }

xPredictCompare <- function(list_pPerf, displayBy=c("ROC","PR"), type=c("curve","bar"), sort=TRUE, detail=TRUE, facet=FALSE, font.family="sans", signature=TRUE)
{

    displayBy <- match.arg(displayBy)
    type <- match.arg(type)
    
   	if(is(list_pPerf,"ggplot") | is(list_pPerf,"pPerf")){
		list_pPerf <- list(list_pPerf)
	}else if(is(list_pPerf,"list")){
		## Remove null elements in a list
		list_pPerf <- base::Filter(base::Negate(is.null), list_pPerf)
		if(length(list_pPerf)==0){
			warnings("All pPerf objects are NULL!")
			return(NULL)
		}
	}else{
		stop("The function must apply to a 'list' object, or a 'pPerf'/'ggplot' object.\n")
	}
    
	## Combine into a data frame called 'df_PRS'
	list_names <- names(list_pPerf)
	if(is.null(list_names)){
		list_names <- paste('Method', 1:length(list_pPerf), sep=' ')
	}
	ls_prs <- lapply(1:length(list_pPerf), function(i){
		prs <- list_pPerf[[i]]$PRS
		fmax <- signif(list_pPerf[[i]]$Fmax, digits=3)
		auroc <- signif(list_pPerf[[i]]$AUROC, digits=3)
		method <- list_names[i]
		
		#label <- paste(method, ' (AUC=', auroc, '; Fmax=', fmax,')', sep='')
		if(displayBy=='ROC'){
			label <- paste(method, ' (AUC=', auroc, ')', sep='')
		}else if(displayBy=='PR'){
			label <- paste(method, ' (Fmax=', fmax,')', sep='')
		}
		
		cbind(prs, Method=rep(method,nrow(prs)), fmax=rep(fmax,nrow(prs)), auroc=rep(auroc,nrow(prs)), Label=rep(label,nrow(prs)), stringsAsFactors=FALSE)
	})
	df_PRS <- do.call(rbind, ls_prs)

	## Method factor
	df_PRS$Method <- factor(df_PRS$Method, levels=list_names)
	
	if(type=='curve'){
		## draw curves
		Recall <- Precision <- Specificity <- Method <- Label <- auroc <- NULL
		if(displayBy=='ROC'){
			## sort by: auroc
			if(sort){
				df_PRS <- df_PRS[with(df_PRS, order(-auroc)), ]
				## define levels
				if(detail){
					df_PRS$Label <- factor(df_PRS$Label, levels=unique(df_PRS$Label))
				}else{
					df_PRS$Method <- factor(df_PRS$Method, levels=unique(df_PRS$Method))
				}
			}
		
			## ggplot
			p <- ggplot(df_PRS, aes(x=1-Specificity,y=Recall))
		
			if(detail){
				p <- p + geom_line(aes(colour=factor(Label)))
			}else{
				p <- p + geom_line(aes(colour=factor(Method)))
			}
		
			p <- p + ylab("True Positive Rate = TP/(TP+FN)") + xlab("False Positive Rate = FP/(FP+TN)") + ylim(0,max(df_PRS$Recall)) + xlim(0,max(1-df_PRS$Specificity))
		
		}else if(displayBy=='PR'){
			## sort by: fmax
			fmax <- ''
			if(sort){
				df_PRS <- df_PRS[with(df_PRS, order(-fmax)), ]
				## define levels
				if(detail){
					df_PRS$Label <- factor(df_PRS$Label, levels=unique(df_PRS$Label))
				}else{
					df_PRS$Method <- factor(df_PRS$Method, levels=unique(df_PRS$Method))
				}
			}
			## ggplot
			p <- ggplot(df_PRS, aes(x=Recall,y=Precision)) 

			if(detail){
				p <- p + geom_line(aes(colour=factor(Label)))
			}else{
				p <- p + geom_line(aes(colour=factor(Method)))
			}

			p <- p + ylab("Precision = TP/(TP+FP)") + xlab("Recall = TP/(TP+FN)") + ylim(0,max(df_PRS$Precision)) + xlim(0,max(df_PRS$Recall))
		
		}
	
		p <- p + theme_bw() + theme(axis.title.y=element_text(size=12,color="black"), axis.title.x=element_text(size=12,color="black"))
	
		if(facet){
			if(detail){
				p <- p + facet_wrap(~Label)
			}else{
				p <- p + facet_wrap(~Method)
			}
		
			## strip
			p <- p + theme(strip.background=element_rect(fill="transparent",color="transparent"), strip.text=element_text(face="italic"))
		
			p <- p + theme(legend.position="none", legend.title=element_blank())
		}else{
			p <- p + theme(legend.title=element_blank(), legend.key=element_rect(colour="transparent"))
		
			#p + theme(legend.position=c(0.75,0.25))
		}
	
	}else if(type=='bar'){
		
		df <- unique(df_PRS[,c("Method","fmax","auroc","Label")])
		
		## draw bar
		Method <- Label <- auroc <- NULL
		if(displayBy=='ROC'){
			## sort by: auroc
			if(sort){
				df <- df[with(df, order(-auroc)), ]
				## define levels
				df$Method <- factor(df$Method, levels=unique(df$Method))
			}
		
			## ggplot
			p <- ggplot(df, aes(x=Method,y=auroc))
		
			if(detail){
				p <- p + geom_col(aes(fill=factor(Label)))
			}else{
				p <- p + geom_col(aes(fill=factor(Method)))
			}
			p <- p + ylab("AUC\n(a measure of ROC)")
			p <- p + geom_text(aes(label=auroc), hjust=1)
			
			#ylim_low <- ifelse(min(df$auroc)>0.5, 0.5, min(df$auroc))
			#p <- p + coord_cartesian(ylim=c(ylim_low,1))
			
		}else if(displayBy=='PR'){
			## sort by: fmax
			if(sort){
				df <- df[with(df, order(-fmax)), ]
				## define levels
				df$Method <- factor(df$Method, levels=unique(df$Method))
			}
		
			## ggplot
			p <- ggplot(df, aes(x=Method,y=fmax))
		
			if(detail){
				p <- p + geom_col(aes(fill=factor(Label)))
			}else{
				p <- p + geom_col(aes(fill=factor(Method)))
			}
			p <- p + ylab("F-max\n(a measure of Precision-Recall curve)")
			p <- p + geom_text(aes(label=fmax), hjust=1)
		}
		
		p <- p + theme_bw() + theme(legend.position="none",axis.title.y=element_blank(), axis.text.y=element_text(size=12,color="black"), axis.title.x=element_text(size=14,color="black")) + coord_flip()
	
	}
	
	if(signature){
		caption <- paste("Created by xPredictROCR from Pi version", utils ::packageVersion("Pi"))
		p <- p + labs(caption=caption) + theme(plot.caption=element_text(hjust=1,face='bold.italic',size=8,colour='#002147'))
	}
	
	p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	
	## change font family to 'Arial'
	p <- p + theme(text=element_text(family=font.family))
	
	## put arrows on both axes
	p <- p + theme(axis.line=element_line(arrow=arrow(angle=30,length=unit(0.25,"cm"), type="open")))

	
	invisible(p)
}
