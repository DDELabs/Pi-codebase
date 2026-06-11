#' Function to assess the prediction performance via ROC and Precision-Recall (PR) analysis
#'
#' \code{xPredictROCR} is supposed to assess the prediction performance via Receiver Operating Characteristic (ROC) and Precision-Recall (PR) analysis. It requires three inputs: 1) Gold Standard Positive (GSP) targets; 2) Gold Standard Negative (GSN) targets; 3) prediction containing predicted targets and predictive scores.
#'
#' @param prediction a data frame containing predictions along with predictive scores. It has two columns: 1st column for target, 2nd column for predictive scores (the higher the better). Alternatively, it can be an object of class "pNode" (or "sTarget" or "dTarget") from which a data frame is extracted
#' @param GSP a vector containing Gold Standard Positives (GSP)
#' @param GSN a vector containing Gold Standard Negatives (GSN)
#' @param rescale logical to indicate whether to linearly rescale predictive scores for GSP/GSN targets to the range [0,1]. By default, it sets to TRUE
#' @param plot the way to plot performance curve. It can be 'none' for no curve returned, 'ROC' for ROC curve, and 'PR' for PR curve. 
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to TRUE for display
#' @param font.family the font family for texts
#' @param signature a logical to indicate whether the signature is assigned to the plot caption. By default, it sets TRUE showing which function is used to draw this graph
#' @return 
#' If plot is 'none' (by default), an object of class "pPerf", a list with following components:
#' \itemize{
#'  \item{\code{PRS}: a data frame with 3 columns ('Precision', 'Recall' and 'Specificity')}
#'  \item{\code{AUROC}: a scalar value for ROC AUC}
#'  \item{\code{Fmax}: a scalar value for maximum F-measure}
#'  \item{\code{ROC_perf}: a ROCR performance-class object for ROC curve}
#'  \item{\code{PR_perf}: a ROCR performance-class object for PR curve}
#'  \item{\code{Pred_obj}: a ROCR prediction-class object (potentially used for calculating other performance measures)}
#' }
#' If plot is 'ROC' or 'PR', it will return a ggplot object after being appended with the same components as mentioned above.
#' If no GSP and/or GSN is predicted, it will return NULL
#' @note
#' AUC: the area under ROC
#' F-measure: the maximum of a harmonic mean between precision and recall along PR curve
#' @export
#' @seealso \code{\link{xPredictROCR}}
#' @include xPredictROCR.r
#' @examples
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' \dontrun{
#' pPerf <- xPredictROCR(prediction, GSP, GSN)
#' }

xPredictROCR <- function(prediction, GSP, GSN, rescale=TRUE, plot=c("none","ROC","PR"), verbose=TRUE, font.family="sans", signature=TRUE)
{

    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    plot <- match.arg(plot)

    if (is(prediction,"pNode")){
        prediction <- prediction$priority[,c("name","priority")]
    }else if(is(prediction,"sTarget")){
    	prediction <- prediction$priority[, c("name","rating")]
    }else if(is(prediction,"dTarget")){
    	prediction <- prediction$priority[, c("name","rating")]
    }

    res_ls <- split(x=as.numeric(prediction[,2]), f=as.character(prediction[,1]))
    pred <- unlist(lapply(res_ls, base::max))
    ####
    pred <- pred[!is.na(pred)]
    ####
    if(verbose){
        now <- Sys.time()
        message(sprintf("There are %d targets in predictions (%s).", length(pred), as.character(now)), appendLF=TRUE)
    }
	
	## GSP
    gsp <- as.character(unique(GSP))
	### GSP but only predicted
    ind <- match(gsp, names(pred))
    gsp_predicted <- gsp[!is.na(ind)]
    if(verbose){
        now <- Sys.time()
        message(sprintf("Of %d targets in GSP, %d also predicted for evaluation (%s).", length(gsp), length(gsp_predicted), as.character(now)), appendLF=TRUE)
    }
    
	## GSN
    gsn <- as.character(unique(GSN))
	### GSN but only predicted
    ind <- match(gsn, names(pred))
    gsn_predicted <- gsn[!is.na(ind)]
    if(verbose){
        now <- Sys.time()
        message(sprintf("Of %d targets in GSN, %d also predicted for evaluation (%s).", length(gsn), length(gsn_predicted), as.character(now)), appendLF=TRUE)
    }
    
    ########################################
    # NULL if no GSP and/or GSN is predicted
    ########################################
    if(length(gsp_predicted)==0 | length(gsn_predicted)==0){
    	warnings("No GSP and/or GSN is predicted!")
    	return(NULL)
    }
    ########################################
    
    ######################################
	## prepare input for ROCR	
	gsp_pred_label <- data.frame(pred=pred[gsp_predicted], label=rep(1,length(gsp_predicted)), stringsAsFactors=FALSE)
	gsn_pred_label <- data.frame(pred=pred[gsn_predicted], label=rep(0,length(gsn_predicted)), stringsAsFactors=FALSE)
	pred_label <- rbind(gsp_pred_label, gsn_pred_label)
	
	## whether prediction scores are rescaled to the range [0,1]
	if(rescale){
		x <- pred_label$pred
		x_scaled <- (x - min(x)) / (max(x) - min(x))
		pred_label$pred <- x_scaled
	}
	
	## ROCR
	if(length(suppressWarnings(tryCatch(pred_obj <- ROCR::prediction(predictions=pred_label$pred, labels=pred_label$label), error=function(e) e, warning=function(w) w)))==2){
		return(NULL)
	}
	#pred_obj <- ROCR::prediction(predictions=pred_label$pred, labels=pred_label$label)

	## acc: Accuracy=(TP+TN)/(P+N)
	res <- ROCR::performance(pred_obj, measure="acc")
	acc <- unlist(res@y.values)
	
	## auc: Area under the ROC curve
	res <- ROCR::performance(pred_obj, measure="auc")
	auroc <- unlist(res@y.values)
	
	## f: Precision-recall F measure
	res <- ROCR::performance(pred_obj, measure="f", x.measure="rec")
	fmax <- base::max(unlist(res@y.values), na.rm=TRUE)
    i <- which(unlist(res@y.values)==fmax)[1]
    
	## ROC curves
	perf_roc <- ROCR::performance(pred_obj, measure="tpr", x.measure="fpr")
    tpr <- unlist(perf_roc@y.values)
    fpr <- unlist(perf_roc@x.values)
    
    ## PR curves
	perf_pr <- ROCR::performance(pred_obj, measure="prec", x.measure="rec")
    prec <- unlist(perf_pr@y.values)
    rec <- unlist(perf_pr@x.values)
    
    ## TPR vs pcfall (Prediction-conditioned fallout: FP/(TP+FP))
	#perf_tprfdr <- ROCR::performance(pred_obj, measure="tpr", x.measure="pcfall")    
    #tpr <- unlist(perf_tprfdr@y.values)
    #fdr <- unlist(perf_tprfdr@x.values)
    ##############################################################################################
    
    if(verbose){
        message(sprintf("In summary, Area under ROC: %.3f, and PR F-max: %.3f.", auroc, fmax), appendLF=TRUE)
    }
	
	df_PRS <- data.frame(Precision=prec, Recall=rec, Specificity=1-fpr, Accurancy=acc)
    
    Recall <- Precision <- Specificity <- NULL
    if(plot=='none'){
    	pPerf <- list(
    		PRS=df_PRS,
    		AUROC=auroc,
    		Fmax=fmax,
    		ROC_perf=perf_roc,
    		PR_perf=perf_pr,
    		Pred_obj=pred_obj
    	)
    	class(pPerf) <- "pPerf"
    	
    	invisible(pPerf)
    	
    }else if(plot=='PR'){
		p <- ggplot(df_PRS, aes(x=Recall,y=Precision)) 
		p <- p + geom_line() + theme_bw() + ylab("Precision = TP/(TP+FP)") + xlab("Recall = TP/(TP+FN)") + ylim(0,max(df_PRS$Precision)) + xlim(0,max(df_PRS$Recall)) 
		
		## title
		#p <- p + ggtitle(paste0('PR curve (Fmax = ',signif(fmax,digits=3),')'))
		title <- 'PR curve'
		p <- p + labs(title=title) + theme(plot.title=element_text(hjust=0.5))
		## caption
		if(signature){
			caption <- paste("Created by xPredictROCR from Pi version", utils ::packageVersion("Pi"))
			p <- p + labs(caption=caption) + theme(plot.caption=element_text(hjust=1,face='bold.italic',size=8,colour='#002147'))
		}
		
		## Add the point with the maximum F
		p <- p + geom_point(data=df_PRS[i,], aes(x=Recall,y=Precision), colour="red")
		#p <- p + geom_text_repel(data=df_PRS[i,], family="Times New Roman", aes(x=Recall, y=Precision, label=paste0('Fmax = MAX{2*P*R/(P+R)} = ',signif(fmax,digits=3))))
		#p <- p + geom_text_repel(data=df_PRS[i,], family="Times New Roman", aes(x=Recall, y=Precision, label=paste0('Fmax = ',signif(fmax,digits=3))))
		p <- p + geom_text(data=df_PRS[i,], family="Times New Roman", aes(x=Recall, y=Precision, label=paste0('Fmax = ',signif(fmax,digits=3))))
		
		## change font family to 'Arial'
		p <- p + theme(text=element_text(family=font.family))
		
		## put arrows on both axes
		p <- p + theme(axis.line=element_line(arrow=arrow(angle=30,length=unit(0.25,"cm"), type="open")))
		
		p$PRS <- df_PRS
		p$AUROC <- auroc
		p$Fmax <- fmax
		p$ROC_perf <- perf_roc
		p$PR_perf <- perf_pr
		p$Pred_obj <- pred_obj
		
		invisible(p)
		
    }else if(plot=='ROC'){
		p <- ggplot(df_PRS, aes(x=1-Specificity,y=Recall)) 
		p <- p + geom_line() + theme_bw() + ylab("True Positive Rate = TP/(TP+FN)") + xlab("False Positive Rate = FP/(FP+TN)") + ylim(0,max(df_PRS$Recall)) + xlim(0,max(1-df_PRS$Specificity))
		
		## title
		title <- paste0('ROC (AUC = ',signif(auroc,digits=3),')')
		p <- p + labs(title=title) + theme(plot.title=element_text(hjust=0.5))
		## caption
		if(signature){
			caption <- paste("Created by xPredictROCR from Pi version", utils ::packageVersion("Pi"))
			p <- p + labs(caption=caption) + theme(plot.caption=element_text(hjust=1,face='bold.italic',size=8,colour='#002147'))
		}
		
		## change font family to 'Arial'
		p <- p + theme(text=element_text(family=font.family))
		
		## put arrows on both axes
		p <- p + theme(axis.line=element_line(arrow=arrow(angle=30,length=unit(0.25,"cm"), type="open")))
		
		p$PRS <- df_PRS
		p$AUROC <- auroc
		p$Fmax <- fmax
		p$ROC_perf <- perf_roc
		p$PR_perf <- perf_pr
		p$Pred_obj <- pred_obj
		
		invisible(p)
	}

}


