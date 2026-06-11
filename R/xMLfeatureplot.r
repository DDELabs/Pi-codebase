#' Function to visualise/assess features used for machine learning
#'
#' \code{xMLfeatureplot} is supposed to visualise/assess features used for machine learning. Visualisation can be made using either boxplot or dot plot for AUC and F-max. It returns an object of class "ggplot" for AUC and F-max, and an object of class "trellis" for boxplot.
#'
#' @param df_predictor a data frame containing genes (in rows) and predictors (in columns), with their predictive scores inside it. This data frame must has gene symbols as row names
#' @param GSP a vector containing Gold Standard Positive (GSP)
#' @param GSN a vector containing Gold Standard Negative (GSN)
#' @param displayBy which statistics will be used for displaying. It can be either "boxplot" for features themselves, "ROC" for AUC in ROC, "Fmax" for F-max in Precision-Recall curve)
#' @param font.family the font family for texts
#' @param ... additional parameters. Please refer to 'lattice::bwplot' for the complete list.
#' @return an object of class "ggplot" for AUC and F-max, and an object of class "trellis" for boxplot
#' @note none
#' @export
#' @seealso \code{\link{xPredictROCR}}, \code{\link{xPredictCompare}}
#' @include xMLfeatureplot.r
#' @examples
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' \dontrun{
#' gp <- xMLfeatureplot(df_predictor, GSP, GSN, displayBy="ROC")
#' }

xMLfeatureplot <- function(df_predictor, GSP, GSN, displayBy=c("boxplot","ROC","Fmax"), font.family="sans", ...)
{

    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    displayBy <- match.arg(displayBy)

	## pre-process GSP and GSN
	gsp <- unique(GSP)
	gsn <- unique(GSN)
	gsn <- setdiff(gsn, gsp)
	gs_names <- union(gsp, gsn)
	gs_targets <- rep(0, length(gs_names))
	names(gs_targets) <- gs_names
	gs_targets[gsp] <- 1
	
	## predictors + class
	ind <- match(names(gs_targets), rownames(df_predictor))
	df_predictor_class <- as.data.frame(df_predictor[ind[!is.na(ind)],])
	class <- as.factor(gs_targets[!is.na(ind)])
	levels(class) <- c("GSN","GSP")
	df_predictor_class$class <- class
	
	if(displayBy=="boxplot"){

		tmp <- colnames(df_predictor_class)
		tmp <- gsub('^Supervised_', 'Supervised\n(', tmp)
		tmp <- gsub('^Annotation_', 'Annotation\n(', tmp)
		tmp <- gsub('^nearbyGenes_', 'nearbyGenes\n(', tmp)
		tmp <- gsub('^eQTL_', 'eQTL\n(', tmp)
		tmp <- gsub('^HiC_', 'Hi-C\n(', tmp)
		tmp <- gsub('^Expression_', 'Expression\n(', tmp)
		tmp <- paste(tmp,')',sep='')
		tmp <- gsub('\n', ' ', tmp)
		colnames(df_predictor_class) <- tmp
	
		strip.left.aligned <- function(which.given, which.panel, factor.levels, ...){
			lattice::panel.rect(0, 0, 1, 1, col="transparent", border=0)
			lattice::panel.text(x=0, y=0.5, pos=4, lab=factor.levels[which.panel[which.given]], cex=0.45, col="black", font=2, srt=5)
		}
	
		res <- caret::featurePlot(x=df_predictor_class[,-ncol(df_predictor_class)],
					y=df_predictor_class$class, 
					plot="box",
					## Pass in options to bwplot() 
					scales = list(x=list(relation="free",rot=0,cex=0.45), y=list(relation="free",log=TRUE,rot=0,cex=0.35)),
					#strip = lattice::strip.custom(bg="lightgrey"),
					strip = strip.left.aligned,
					par.settings = list(axis.line=list(col="grey")),
					labels = c("","Affinity score"),
					#layout = c(6,6),
					...
				)
    }else{
	
		######################
		## evaluation
		######################
		### do preparation
		df_pred <- df_predictor
		ls_predictors <- lapply(colnames(df_pred), function(x){
			data.frame(rownames(df_pred), df_pred[,x], stringsAsFactors=FALSE)
		})
		names(ls_predictors) <- colnames(df_pred)
		# do evaluation
		ls_pPerf <- lapply(ls_predictors, function(x){
			pPerf <- xPredictROCR(prediction=x, GSP=GSP, GSN=GSN, verbose=FALSE)
		})
		# do plotting
		bp <- xPredictCompare(ls_pPerf, displayBy=c("ROC","PR"))
		df <- unique(bp$data[,c('Method','auroc','fmax')])
		df_performance <- cbind(ROC=df$auroc, Fmax=df$fmax)
		rownames(df_performance) <- df$Method
		
		#############
		if(is.null(df_performance)){
			return(NULL)
		}
		#############
		
		if(displayBy=='ROC'){
			df <- data.frame(Val=df_performance[,1], stringsAsFactors=FALSE)
			rownames(df) <- rownames(df_performance)
			df$Val <- ifelse(df$Val>=0.5, df$Val, 1-df$Val)
			xlab <- "AUC\n(a measure of ROC)"
		}else if(displayBy=='Fmax'){
			df <- data.frame(Val=df_performance[,2], stringsAsFactors=FALSE)
			rownames(df) <- rownames(df_performance)
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
    
		## order by 'Predictor', 'Val'
		df <- df[with(df,order(Predictor,-Val)),]
		df$Predictor <- factor(df$Predictor, levels=unique(df$Predictor))
		df$Method <- factor(df$Method, levels=rev(unique(df$Method)))
		bp <- ggplot(df, aes(Val, Method, colour=Predictor))
		bp <- bp + geom_point()
	
		bp <- bp + scale_color_manual(values=xColormap("ggplot2")(length(levels(df$Predictor))))
		bp <- bp + theme_bw() + theme(legend.position="right", legend.title=element_blank(), axis.title.y=element_blank(), axis.text.y=element_text(size=10,color="black"), axis.title.x=element_text(size=14,color="black"), panel.background=element_rect(fill=rgb(0.95,0.95,0.95,1)))
		bp <- bp + xlab(xlab)
	
		## change font family to 'Arial'
		bp <- bp + theme(text=element_text(family=font.family))
	
		## put arrows on x-axis
		bp <- bp + theme(axis.line.x=element_line(arrow=arrow(angle=30,length=unit(0.25,"cm"), type="open")))
	
		## x-axis position
		if(displayBy == "ROC"){
			bp <- bp + scale_x_continuous(position="top", limits=c(0.5,1))
		}else{
			bp <- bp + scale_x_continuous(position="top")
		}
		
		res <- bp
    
    }
    
    invisible(res)
}


