#' Function to visualise machine learning results using dot plot
#'
#' \code{xMLdotplot} is supposed to visualise machine learning results using dot plot. It returns an object of class "ggplot".
#'
#' @param sTarget an object of class "sTarget"
#' @param displayBy which statistics will be used for displaying. It can be either statistics across folds ("importance2fold" for predictor importance, "roc2fold" for AUC in ROC, "fmax2fold" for F-max in Precision-Recall curve) or overall statistics ("importance_accurancy" for predictor importance measured by accuracy decrease, "importance_gini" for predictor importance measured by Gini decrease, "ROC" for AUC in ROC, "Fmax" for F-max in Precision-Recall curve)
#' @param font.family the font family for texts
#' @param ML.included logical to indicate whether to use ML results
#' @param signature logical to indicate whether the signature is assigned to the plot caption. By default, it sets TRUE showing which function is used to draw this graph
#' @return an object of class "ggplot"
#' @note none
#' @export
#' @seealso \code{\link{xMLdotplot}}
#' @include xMLdotplot.r
#' @examples
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' \dontrun{
#' gp <- xMLdotplot(sTarget, displayBy="importance_accurancy")
#' gp
#' }

xMLdotplot <- function(sTarget, displayBy=c("importance2fold","roc2fold","fmax2fold","importance_accurancy","importance_gini","ROC","Fmax"), ML.included=TRUE, font.family="sans", signature=TRUE) 
{
    
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    displayBy <- match.arg(displayBy)

    if(!is(sTarget,"sTarget")){
    	stop("The function must apply to a 'sTarget' object.\n")
    }
    
    nfold <- length(sTarget$model)
    
    if(displayBy == "importance2fold"){
    	df <- sTarget$importance2fold
    	xlab <- paste0("Decrease in accuracy across ", nfold, " folds\n(a measure of predictor importance)")
    }else if(displayBy=='roc2fold'){
    	df <- sTarget$roc2fold
    	xlab <- paste0("AUC across ", nfold, " folds\n(a measure of ROC)")
    }else if(displayBy=='fmax2fold'){
    	df <- sTarget$fmax2fold
    	xlab <- paste0("F-max across ", nfold, " folds\n(a measure of Precision-Recall curve)")
    }else if(displayBy=='importance_accurancy'){
    	df <- data.frame(Val=sTarget$importance[,1], stringsAsFactors=FALSE)
    	rownames(df) <- rownames(sTarget$importance)
    	xlab <- "Decrease in accuracy\n(a measure of predictor importance)"
    }else if(displayBy=='importance_gini'){
    	df <- data.frame(Val=sTarget$importance[,2], stringsAsFactors=FALSE)
    	rownames(df) <- rownames(sTarget$importance)
    	xlab <- "Decrease in gini\n(a measure of predictor importance)"
    }else if(displayBy=='ROC'){
    	#### replace with roc2fold for Supervised_randomforest
    	Val <- sTarget$performance[,1]
    	if(ML.included){
    		Val[1] <- sTarget$roc2fold[1,1]
    	}else{
    		Val <- Val[-1]
    	}
    	############
    	df <- data.frame(Val=Val, stringsAsFactors=FALSE)
    	#rownames(df) <- rownames(sTarget$performance)
    	xlab <- "AUC\n(a measure of ROC)"
    }else if(displayBy=='Fmax'){
    	#### replace with fmax2fold for Supervised_randomforest
    	Val <- sTarget$performance[,2]
    	Val[1] <- sTarget$fmax2fold[1,1]
    	############
    	df <- data.frame(Val=Val, stringsAsFactors=FALSE)
    	rownames(df) <- rownames(sTarget$performance)
    	xlab <- "F-max\n(a measure of Precision-Recall curve)"
    }
    
    median <- name <- max <- min <- Val <- ''
    
    ## extract info on 'Predictor' and 'Method'
	tmp <- rownames(df)
	tmp <- gsub('^Supervised_', 'Supervised\n(', tmp)
	tmp <- gsub('^Annotation_', 'Annotation\n(', tmp)
	tmp <- gsub('^nearbyGenes_', 'nearbyGenes\n(', tmp)
	tmp <- gsub('^eQTL_', 'eQTL\n(', tmp)
	tmp <- gsub('^HiC_', 'Hi-C\n(', tmp)
	tmp <- gsub('^Expression_', 'Expression\n(', tmp)
	tmp <- paste(tmp,')',sep='')
	df <- data.frame(name=tmp, df, stringsAsFactors=FALSE)
	Predictor <- gsub('\n.*', '', as.character(df$name), perl=TRUE)
	Method <- gsub('.*\n\\(|\\)$', '', as.character(df$name), perl=TRUE)
	df$Predictor <- Predictor
	df$Method <- Method
	
	if(displayBy == "importance2fold" | displayBy == "roc2fold" | displayBy == "fmax2fold"){
		## order by 'Predictor', 'median'
		df <- df[with(df,order(Predictor,-median)),]
		df$Predictor <- factor(df$Predictor, levels=unique(df$Predictor))
		df$Method <- factor(df$Method, levels=rev(unique(df$Method)))
		bp <- ggplot(df, aes(median, Method, colour=Predictor))
		bp <- bp + geom_point() + geom_errorbarh(aes(xmax=max, xmin=min, height=.2))
	}else{	
		## order by 'Predictor', 'Val'
		df <- df[with(df,order(Predictor,-Val)),]
		df$Predictor <- factor(df$Predictor, levels=unique(df$Predictor))
		df$Method <- factor(df$Method, levels=rev(unique(df$Method)))
		bp <- ggplot(df, aes(Val, Method, colour=Predictor))
		bp <- bp + geom_point()
	}
	
	bp <- bp + scale_color_manual(values=xColormap("ggplot2")(length(levels(df$Predictor))))
	bp <- bp  + theme_bw() + theme(legend.position="right", legend.title=element_blank(), axis.title.y=element_blank(), axis.text.y=element_text(size=10,color="black"), axis.title.x=element_text(size=14,color="black"), panel.background=element_rect(fill=rgb(0.95,0.95,0.95,1)))
	bp <- bp + xlab(xlab)
	
	## caption
    if(signature){
    	caption <- paste("Created by xMLdotplot from Pi version", utils ::packageVersion("Pi"))
    	bp <- bp + labs(caption=caption) + theme(plot.caption=element_text(hjust=1,face='bold.italic',size=8,colour='#002147'))
    }
	
	## change font family to 'Arial'
	bp <- bp + theme(text=element_text(family=font.family))
	
	## put arrows on x-axis
	bp <- bp + theme(axis.line.x=element_line(arrow=arrow(angle=30,length=unit(0.25,"cm"), type="open")))
	
	## x-axis position
    if(displayBy == "ROC" | displayBy == "roc2fold"){
    	bp <- bp + scale_x_continuous(position="top", limits=c(0.5,1))
    }else{
		bp <- bp + scale_x_continuous(position="top")
	}
	
	invisible(bp)
}
