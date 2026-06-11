#' Function to integrate predictor matrix in a supervised manner via machine learning algorithm glmnet.
#'
#' \code{xMLglmnet} is supposed to integrate predictor matrix in a supervised manner via machine learning algorithm glmnet. It requires three inputs: 1) Gold Standard Positive (GSP) targets; 2) Gold Standard Negative (GSN) targets; 3) a predictor matrix containing genes in rows and predictors in columns, with their predictive scores inside it. It returns an object of class 'pTarget'.
#'
#' @param df_predictor a data frame containing genes (in rows) and predictors (in columns), with their predictive scores inside it. This data frame must has gene symbols as row names
#' @param GSP a vector containing Gold Standard Positive (GSP)
#' @param GSN a vector containing Gold Standard Negative (GSN)
#' @param family response family type. It can be one of "binomial" for two-class logistic model or "gaussian" for gaussian model
#' @param type.measure loss to use for cross-validation. It can be one of "auc" for two-class logistic model, "mse" for the deviation from the fitted mean to the response using gaussian model
#' @param nfold an integer specifying the number of folds for cross validataion
#' @param alphas a vector specifying a range of alphas. Alpha is an elasticnet mixing parameter, with 0<=alpha<=1. By default, seq(0,1,by-0.1)
#' @param standardize logical specifying whether to standardise the predictor. If yes (by default), the predictor standardised prior to fitting the model. The coefficients are always returned on the original scale
#' @param lower.limits vector of lower limits for each coefficient (by default, '-Inf'; all should be non-positive). A single value provided will apply to every coefficient
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to TRUE for display
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @param guid a valid (5-character) Global Unique IDentifier for an OSF project. See \code{\link{xRDataLoader}} for details
#' @param ... additional parameters. Please refer to 'glmnet::cv.glmnet' for the complete list.
#' @return 
#' an object of class "pTarget", a list with following components:
#' \itemize{
#'  \item{\code{model}: an object of class "cv.glmnet" as a best model}
#'  \item{\code{priority}: a data frame of nGene X 5 containing gene priority information, where nGene is the number of genes in the input data frame, and the 5 columns are "GS" (either 'GSP', or 'GSN', or 'NEW'), "name" (gene names), "rank" (ranks of the priority scores), "priority" (priority score; rescaled into the 5-star ratings), and "description" (gene description)}
#'  \item{\code{predictor}: a data frame, which is the same as the input data frame but inserting an additional column 'GS' in the first column}
#'  \item{\code{cvm2alpha}: a data frame of nAlpha X 2 containing mean cross-validated error, where nAlpha is the number of alpha and the two columns are "min" (lambda.min) and "1se" (lambda.1se)}
#'  \item{\code{nonzero2alpha}: a data frame of nAlpha X 2 containing the number of non-zero coefficients, where nAlpha is the number of alpha and the two columns are "min" (lambda.min) and "1se" (lambda.1se)}
#'  \item{\code{importance}: a data frame of nPredictor X 1 containing the predictor importance/coefficient info}
#'  \item{\code{performance}: a data frame of 1+nPredictor X 2 containing the supervised/predictor performance info predictor importance info, where nPredictor is the number of predictors, two columns are "ROC" (AUC values) and "Fmax" (F-max values)}
#'  \item{\code{gp}: a ggplot object for the ROC curve}
#'  \item{\code{call}: the call that produced this result}
#' }
#' @note none
#' @export
#' @seealso \code{\link{xPredictROCR}}, \code{\link{xPredictCompare}}, \code{\link{xSymbol2GeneID}}
#' @include xMLglmnet.r
#' @examples
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' \dontrun{
#' pTarget <- xMLglmnet(df_prediction, GSP, GSN)
#' }

xMLglmnet <- function(df_predictor, GSP, GSN, family=c("binomial","gaussian"), type.measure=c("auc","mse"), nfold=3, alphas=seq(0,1,0.1), standardize=TRUE, lower.limits=-Inf, verbose=TRUE, RData.location="http://galahad.well.ox.ac.uk/bigdata", guid=NULL, ...)
{
	
    startT <- Sys.time()
    if(verbose){
        message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=TRUE)
        message("", appendLF=TRUE)
    }
    ####################################################################################
	
    family <- match.arg(family)
    type.measure <- match.arg(type.measure)
	
	if(family=="gaussian"){
		type.measure <- "mse"
	}
	
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
	df_predictor_class$class <- class
	
    if(verbose){
        now <- Sys.time()
        message(sprintf("1. Gold standards (%d in GSP, %d in GSN) are used for supervised integration of %d predictors/features (%s).", sum(class==1), sum(class==0), ncol(df_predictor), as.character(now)), appendLF=TRUE)
    }
	
	yTrain <- as.matrix(as.numeric(df_predictor_class[, ncol(df_predictor_class)])-1, ncol=1)
	xTrain <- as.matrix(df_predictor_class[, -ncol(df_predictor_class)])
	
	#####################################
    if(verbose){
        now <- Sys.time()
        message(sprintf("2. Do glmnet (%d-fold cross validation) using '%s' family measured by '%s' (%s) ...", nfold, family, type.measure, as.character(now)), appendLF=TRUE)
    }
	
	## elastic net with different tuning alpha parameters
	alphas <- sort(unique(alphas))
	alphas <- alphas[alphas>=0 & alphas<=1]
	if(length(alphas)==0){
		alphas <- seq(0,1,0.1)
	}
	names(alphas) <- alphas
	ls_cvfit <- lapply(alphas, function(alpha){
		if(verbose){
			message(sprintf("\talpha=%s", alpha), appendLF=TRUE)
		}
		
		set.seed(alpha*100)
		#cvfit <- glmnet::cv.glmnet(xTrain, yTrain, nfold=nfold, family=family, type.measure=type.measure, alpha=alpha/100, parallel=FALSE, nlambda=100, standardize=standardize, lower.limits=lower.limits, ...)
		cvfit <- glmnet::cv.glmnet(xTrain, yTrain, nfold=nfold, family=family, type.measure=type.measure, alpha=alpha/100, parallel=FALSE, nlambda=100, standardize=standardize, lower.limits=lower.limits)
		#coef(cvfit, s=c(cvfit$lambda.min, cvfit$lambda.1se))
		
		## non-zero predictors
		nonzero.min <- predict(cvfit, newx=xTrain, s="lambda.min", type="nonzero")
		nonzero.1se <- predict(cvfit, newx=xTrain, s="lambda.1se", type="nonzero")
		if(is.null(nrow(nonzero.min))){
			n.nonzero.min <- 0
		}else{
			n.nonzero.min <- nrow(nonzero.min)
		}
		if(is.null(nrow(nonzero.1se))){
			n.nonzero.1se <- 0
		}else{
			n.nonzero.1se <- nrow(nonzero.1se)
		}
		vec_nonzero <- c(n.nonzero.min, n.nonzero.1se)
		names(vec_nonzero) <- c("min", "1se")
		
		## performance
		ind <- which(cvfit$lambda==cvfit$lambda.min)
		cvm.min <- cvfit$cvm[ind]
		ind <- which(cvfit$lambda==cvfit$lambda.1se)
		cvm.1se <- cvfit$cvm[ind]
		vec_cvm <- c(cvm.min, cvm.1se)
		names(vec_cvm) <- c("min", "1se")
		
		cvfit$vec_nonzero <- vec_nonzero
		cvfit$vec_cvm <- vec_cvm
		return(cvfit)
	})
	#####################################
	
	mat_nonzero <- do.call(rbind, lapply(ls_cvfit, function(x) x$vec_nonzero))
	mat_cvm <- do.call(rbind, lapply(ls_cvfit, function(x) x$vec_cvm))
	
	# ?predict.glmnet
	if(type.measure=="auc"){
		vec_ind <- which(mat_cvm==max(mat_cvm), arr.ind=TRUE)[1,]
	}else{
		vec_ind <- which(mat_cvm==min(mat_cvm), arr.ind=TRUE)[1,]
	}
	cvfit_best <- ls_cvfit[[vec_ind[1]]]
    if(verbose){
        message(sprintf("\tBest model: alpha=%s, lambda.%s, %d nonzero predictors.", rownames(mat_cvm)[vec_ind[1]], colnames(mat_cvm)[vec_ind[2]], mat_nonzero[vec_ind[1],vec_ind[2]]), appendLF=TRUE)
    }
	
	## ROC curve
	if(vec_ind[2]==1){
		s <- "lambda.min"
	}else{
		s <- "lambda.1se"
	}
	df_predictor_glmnet <- predict(cvfit_best, newx=xTrain, s=s, type="response")
	x <- data.frame(name=rownames(df_predictor_glmnet), score=df_predictor_glmnet[,1], stringsAsFactors=FALSE)
	pPerf <- xPredictROCR(prediction=x, GSP=GSP, GSN=GSN, verbose=FALSE)
	ls_pPerf <- list(pPerf)
	names(ls_pPerf) <- paste0('alpha=', rownames(mat_cvm)[vec_ind[1]], '; lambda.', colnames(mat_cvm)[vec_ind[2]], '; n=', mat_nonzero[vec_ind[1],vec_ind[2]])
	gp <- xPredictCompare(ls_pPerf, displayBy="ROC", facet=TRUE, signature=FALSE)
	
	#####################################
    if(verbose){
        now <- Sys.time()
        message(sprintf("3. Do prediction for full set (%s).", as.character(now)), appendLF=TRUE)
    }
	
	df_full <- predict(cvfit_best, newx=df_predictor, s=s, type="response")
	vec_full <- df_full[,1]
	names(vec_full) <- rownames(df_full)
	vec_full <- sort(vec_full, decreasing=TRUE)
	
	## get rank
	vec_rank <- rank(-1*vec_full, ties.method="min")
	
	## priority: being rescaled into the [0,5] range
	priority <- vec_full
	vec_priority <- 5 * (priority - min(priority))/(max(priority) - min(priority))
	
	#########################################
	## output
	### df_priority
	output_gs <- rep('NEW', length(vec_priority))
	names(output_gs) <- names(vec_priority)
	ind <- match(names(vec_priority), names(gs_targets))
	output_gs[!is.na(ind)] <- gs_targets[ind[!is.na(ind)]]
	output_gs[output_gs=='0'] <- 'GSN'
	output_gs[output_gs=='1'] <- 'GSP'
	df_priority <- data.frame(GS=output_gs, name=names(vec_priority), rank=vec_rank, priority=vec_priority, stringsAsFactors=FALSE)
	### add description
	df_priority$description <- xSymbol2GeneID(df_priority$name, details=TRUE, RData.location=RData.location, guid=guid)$description
	###
	
	### df_predictor_gs
	ind <- match(names(vec_priority), rownames(df_predictor))
	output_df_predictor <- df_predictor[ind,]
	df_predictor_gs <- data.frame(GS=output_gs, name=names(vec_priority), output_df_predictor, stringsAsFactors=FALSE)
	
####################################################################################
	
	######################
	## overall evaluation
	######################
	### do preparation
	df_predictor_overall <- cbind(Supervised_glmnet=df_priority$priority, df_predictor_gs[,-c(1,2)])
	rownames(df_predictor_overall) <- rownames(df_priority)
	df_pred <- df_predictor_overall
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
	df_evaluation <- cbind(ROC=df$auroc, Fmax=df$fmax)
	rownames(df_evaluation) <- df$Method
	#####################
	#####################
	
    pTarget <- list(model = cvfit_best,
    				priority = df_priority,
    				predictor = df_predictor_gs,
    				cvm2alpha = mat_cvm,
    				nonzero2alpha = mat_nonzero,
    				importance = as.matrix(coef(cvfit_best)[-1,],ncol=1),
    				performance = df_evaluation,
    				gp = gp,
                  	Call     = match.call()
                 )
    class(pTarget) <- "pTarget"
    
    ####################################################################################
    endT <- Sys.time()
    if(verbose){
        message(paste(c("\nFinish at ",as.character(endT)), collapse=""), appendLF=TRUE)
    }
    
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=TRUE)
    
    invisible(pTarget)
}


