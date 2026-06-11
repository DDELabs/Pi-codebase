#' Function to integrate predictor matrix in a supervised manner via machine learning algorithm random forest.
#'
#' \code{xMLrandomforest} is supposed to integrate predictor matrix in a supervised manner via machine learning algorithm random forest. It requires three inputs: 1) Gold Standard Positive (GSP) targets; 2) Gold Standard Negative (GSN) targets; 3) a predictor matrix containing genes in rows and predictors in columns, with their predictive scores inside it. It returns an object of class 'sTarget'.
#'
#' @param list_pNode a list of "pNode" objects or a "pNode" object
#' @param df_predictor a data frame containing genes (in rows) and predictors (in columns), with their predictive scores inside it. This data frame must has gene symbols as row names
#' @param GSP a vector containing Gold Standard Positive (GSP)
#' @param GSN a vector containing Gold Standard Negative (GSN)
#' @param nfold an integer specifying the number of folds for cross validataion. Per fold creates balanced splits of the data preserving the overall distribution for each class (GSP and GSN), therefore generating balanced cross-vallidation train sets and testing sets. By default, it is 3 meaning 3-fold cross validation
#' @param nrepeat an integer specifying the number of repeats for cross validataion. By default, it is 10 indicating the cross-validation repeated 10 times
#' @param seed an integer specifying the seed
#' @param mtry an integer specifying the number of predictors randomly sampled as candidates at each split. If NULL, it will be tuned by `randomForest::tuneRF`, with starting value as sqrt(p) where p is the number of predictors. The minimum value is 3
#' @param ntree an integer specifying the number of trees to grow. By default, it sets to 2000
#' @param fold.aggregateBy the aggregate method used to aggregate results from k-fold cross validataion. It can be either "orderStatistic" for the method based on the order statistics of p-values, or "fishers" for Fisher's method, "Ztransform" for Z-transform method, "logistic" for the logistic method. Without loss of generality, the Z-transform method does well in problems where evidence against the combined null is spread widely (equal footings) or when the total evidence is weak; Fisher's method does best in problems where the evidence is concentrated in a relatively small fraction of the individual tests or when the evidence is at least moderately strong; the logistic method provides a compromise between these two. Notably, the aggregate methods 'Ztransform' and 'logistic' are preferred here
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to TRUE for display
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @param guid a valid (5-character) Global Unique IDentifier for an OSF project. See \code{\link{xRDataLoader}} for details
#' @param ... additional parameters. Please refer to 'randomForest::randomForest' for the complete list.
#' @return 
#' an object of class "sTarget", a list with following components:
#' \itemize{
#'  \item{\code{model}: a list of models, results from per-fold train set}
#'  \item{\code{priority}: a data frame of nGene X 5 containing gene priority information, where nGene is the number of genes in the input data frame, and the 5 columns are "GS" (either 'GSP', or 'GSN', or 'NEW'), "name" (gene names), "rank" (priority rank), "rating" (the 5-star priority score/rating), and "description" (gene description)}
#'  \item{\code{predictor}: a data frame, which is the same as the input data frame but inserting an additional column 'GS' in the first column}
#'  \item{\code{pred2fold}: a list of data frame, results from per-fold test set}
#'  \item{\code{prob2fold}: a data frame of nGene X 2+nfold containing the probability of being GSP, where nGene is the number of genes in the input data frame, nfold is the number of folds for cross validataion, and the first two columns are "GS" (either 'GSP', or 'GSN', or 'NEW'), "name" (gene names), and the rest columns storing the per-fold probability of being GSP}
#'  \item{\code{importance2fold}: a data frame of nPredictor X 4+nfold containing the predictor importance info per fold, where nPredictor is the number of predictors, nfold is the number of folds for cross validataion, and the first 4 columns are "median" (the median of the importance across folds), "mad" (the median of absolute deviation of the importance across folds), "min" (the minimum of the importance across folds), "max" (the maximum of the importance across folds), and the rest columns storing the per-fold importance}
#'  \item{\code{roc2fold}: a data frame of 1+nPredictor X 4+nfold containing the supervised/predictor ROC info (AUC values), where nPredictor is the number of predictors, nfold is the number of folds for cross validataion, and the first 4 columns are "median" (the median of the AUC values across folds), "mad" (the median of absolute deviation of the AUC values across folds), "min" (the minimum of the AUC values across folds), "max" (the maximum of the AUC values across folds), and the rest columns storing the per-fold AUC values}
#'  \item{\code{fmax2fold}: a data frame of 1+nPredictor X 4+nfold containing the supervised/predictor PR info (F-max values), where nPredictor is the number of predictors, nfold is the number of folds for cross validataion, and the first 4 columns are "median" (the median of the F-max values across folds), "mad" (the median of absolute deviation of the F-max values across folds), "min" (the minimum of the F-max values across folds), "max" (the maximum of the F-max values across folds), and the rest columns storing the per-fold F-max values}
#'  \item{\code{importance}: a data frame of nPredictor X 2 containing the predictor importance info, where nPredictor is the number of predictors, two columns for two types ("MeanDecreaseAccuracy" and "MeanDecreaseGini") of predictor importance measures. "MeanDecreaseAccuracy" sees how worse the model performs without each predictor (a high decrease in accuracy would be expected for very informative predictors), while "MeanDecreaseGini" measures how pure the nodes are at the end of the tree (a high score means the predictor was important if each predictor is taken out)}
#'  \item{\code{performance}: a data frame of 1+nPredictor X 2 containing the supervised/predictor performance info predictor performance info, where nPredictor is the number of predictors, two columns are "ROC" (AUC values) and "Fmax" (F-max values)}
#'  \item{\code{evidence}: an object of the class "eTarget", a list with following components "evidence" and "metag"}
#'  \item{\code{list_pNode}: a list of "pNode" objects}
#' }
#' @note none
#' @export
#' @seealso \code{\link{xPierMatrix}}, \code{\link{xSparseMatrix}}, \code{\link{xPredictROCR}}, \code{\link{xPredictCompare}}, \code{\link{xSymbol2GeneID}}
#' @include xMLrandomforest.r
#' @examples
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' \dontrun{
#' sTarget <- xMLrandomforest(df_prediction, GSP, GSN)
#' }

xMLrandomforest <- function(list_pNode=NULL, df_predictor=NULL, GSP, GSN, nfold=3, nrepeat=10, seed=825, mtry=NULL, ntree=1000, fold.aggregateBy=c("logistic","Ztransform","fishers","orderStatistic"), verbose=TRUE, RData.location="http://galahad.well.ox.ac.uk/bigdata", guid=NULL, ...)
{
	
    startT <- Sys.time()
    if(verbose){
        message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=TRUE)
        message("", appendLF=TRUE)
    }
    ####################################################################################
	
    fold.aggregateBy <- match.arg(fold.aggregateBy)
	
	if(is(list_pNode,"list") & is.null(df_predictor)){
		df_predictor <- xPierMatrix(list_pNode, displayBy="score", combineBy="union", aggregateBy="none", RData.location=RData.location, guid=guid)
		
		###########
		## evidence
		eTarget <- xPierMatrix(list_pNode, displayBy="evidence", combineBy="union", aggregateBy="none", verbose=FALSE, RData.location=RData.location, guid=guid)
		###########
				
	}else if(!is.null(df_predictor)){
		df_predictor <- df_predictor
		
		###########
		## evidence
		eTarget <- NULL
		###########
		
	}else{
		stop("The function must apply to 'list' of 'pNode' objects or a 'data.frame'.\n")
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
	
	#####################################
	
	nfold <- as.integer(nfold)
	nrepeat <- as.integer(nrepeat)
    if(verbose){
        now <- Sys.time()
        if(nfold==1){
			message(sprintf("2. GS matrix of %d rows/genes X %d columns (predictors+class) are used as train set (%s) ...", nrow(df_predictor_class), ncol(df_predictor_class), as.character(now)), appendLF=TRUE)
        }else{
        	message(sprintf("2. GS matrix of %d rows/genes X %d columns (predictors+class) are split (by rows/genes) into %d non-redundant sets: each fold '%d/%d' are used as train set and the remaining '1/%d' as test set. These spits are repeated over %d times (%s) ...", nrow(df_predictor_class), ncol(df_predictor_class), nfold, nfold-1, nfold, nfold, nrepeat, as.character(now)), appendLF=TRUE)
        }
    }
    
    ## create non-redundant sets
    # create balanced splits of the data
	# preserve the overall class distribution
	# generate balanced cross-validation groupings from a set of data
	if(TRUE){
		if(!is.null(seed)) set.seed(seed)
		index_sets <- caret::createMultiFolds(y=df_predictor_class$class, k=nfold, times=nrepeat)
		
	}else{
		index_sets <- lapply(1:nrepeat, function(i){
			if(!is.null(seed)) set.seed(seed+i)
			res_ls <- caret::createFolds(y=df_predictor_class$class, k=nfold, list=TRUE, returnTrain=TRUE)
			#length(base::Reduce(union, index_sets))
			names(res_ls) <- paste0("F",1:nfold)
			names(res_ls) <- paste0("R",i,names(res_ls))
			res_ls
		})
		index_sets <- unlist(index_sets, recursive=FALSE)
	}
	
	#####################################
	
    if(verbose){
        now <- Sys.time()
        if(nfold==1){
			message(sprintf("3. Do random forest: '1/%d' as train set (%s) ...", nfold, as.character(now)), appendLF=TRUE)
		}else{
			message(sprintf("3. Do random forest: %d-fold cross-validation (%s) ...", nfold, as.character(now)), appendLF=TRUE)
		}
    }
	
	if(nfold==1 & nrepeat==1){
		ls_model <- lapply(1:length(index_sets), function(i){
			trainset <- df_predictor_class[index_sets[[i]],]
			
			if(verbose){
				message(sprintf("\tFold %d: %d GSP + %d GSN", i, table(trainset$class)[2], table(trainset$class)[1]), appendLF=TRUE)
			}
			
			if(is.null(mtry)){
				mtry <- as.integer(sqrt(ncol(trainset)))
				suppressMessages(df_mtry <- randomForest::tuneRF(x=trainset[,-ncol(trainset)], y=trainset[, ncol(trainset)], ntreeTry=ntree, mtryStart=mtry, stepFactor=2, trace=FALSE, plot=FALSE))
				ind <- which(df_mtry[,2] == min(df_mtry[,2]))[1]
				mtry <- df_mtry[ind,1]
				if(mtry<3){
					mtry <- 3
				}
			}
			set.seed(i)
			#suppressMessages(rf.model <- randomForest::randomForest(class ~ ., data=trainset, importance=TRUE, ntree=ntree, mtry=mtry, ...))
			suppressMessages(rf.model <- randomForest::randomForest(class ~ ., data=trainset, importance=TRUE, ntree=ntree, mtry=mtry))
		})
		
	}else{
		ls_model <- lapply(1:length(index_sets), function(i){
			trainindex <- index_sets[[i]]
			trainset <- df_predictor_class[trainindex,]
		
			if(verbose){
				message(sprintf("\tRepeatFold %d: %d GSP + %d GSN", i, table(trainset$class)[2], table(trainset$class)[1]), appendLF=TRUE)
			}
		
			if(is.null(mtry)){
				mtry <- as.integer(sqrt(ncol(trainset)))
				suppressMessages(df_mtry <- randomForest::tuneRF(x=trainset[,-ncol(trainset)], y=trainset[, ncol(trainset)], ntreeTry=ntree, mtryStart=mtry, stepFactor=2, trace=FALSE, plot=FALSE))
				ind <- which(df_mtry[,2] == min(df_mtry[,2]))[1]
				mtry <- df_mtry[ind,1]
				if(mtry<3){
					mtry <- 3
				}
			}

			set.seed(i)
			#suppressMessages(rf.model <- randomForest::randomForest(class ~ ., data=trainset, importance=TRUE, ntree=ntree, mtry=mtry, ...))
			suppressMessages(rf.model <- randomForest::randomForest(class ~ ., data=trainset, importance=TRUE, ntree=ntree, mtry=mtry))
		})
		names(ls_model) <- names(index_sets)
	
	}
	
	#####################################
	
    if(verbose){
        now <- Sys.time()
        message(sprintf("4. Extract the predictor/feature importance matrix of %d rows/predictors X %d columns/repeats*folds (%s).", ncol(df_predictor_class)-1, nfold*nrepeat, as.character(now)), appendLF=TRUE)
    }
	
	######################
	## importance per fold
	######################
	ls_res <- lapply(1:length(ls_model), function(i){
		rf.model <- ls_model[[i]]
		rf.model.importance <- randomForest::importance(rf.model)[,3]
		df <- data.frame(predictor=names(rf.model.importance), model=rep(names(ls_model)[i],length(rf.model.importance)), importance=rf.model.importance, stringsAsFactors=FALSE)
	})
	df_res <- do.call(rbind, ls_res)
	df_res <- as.matrix(xSparseMatrix(df_res, verbose=FALSE))
	ind <- match(colnames(df_predictor), rownames(df_res))
	if(nfold==1 & nrepeat==1){
		df_res <- as.matrix(df_res[ind,], ncol=nfold)
		colnames(df_res) <- 'Fold1'
	}else{
		df_res <- df_res[ind,]
	}
	vec_median <- apply(df_res, 1, stats::median)
	vec_mad <- apply(df_res, 1, stats::mad)
	vec_min <- apply(df_res, 1, base::min)
	vec_max <- apply(df_res, 1, base::max)
	#####
	res_x <- apply(df_res, 1, stats::t.test)
	vec_mean <- unlist(lapply(res_x, function(x) x$estimate))
	vec_conf_lower <- unlist(lapply(res_x, function(x) x$conf.int[1]))
	vec_conf_upper <- unlist(lapply(res_x, function(x) x$conf.int[2]))
	#####
    df_importance <- cbind(median=vec_median, mad=vec_mad, min=vec_min, max=vec_max, mean=vec_mean, conf_lower=vec_conf_lower, conf_upper=vec_conf_upper, df_res)

	######################
	## overall importance
	######################
	trainset <- df_predictor_class
	if(is.null(mtry)){
		mtry <- as.integer(sqrt(ncol(trainset)))
		suppressMessages(df_mtry <- randomForest::tuneRF(x=trainset[,-ncol(trainset)], y=trainset[, ncol(trainset)], ntreeTry=ntree, mtryStart=mtry, stepFactor=2, trace=FALSE, plot=FALSE))
		ind <- which(df_mtry[,2] == min(df_mtry[,2]))[1]
		mtry <- df_mtry[ind,1]
		if(mtry<3){
			mtry <- 3
		}
	}
	set.seed(1)
	suppressMessages(rf.model.overall <- randomForest::randomForest(class ~ ., data=trainset, importance=TRUE, ntree=ntree, mtry=mtry))
	rf.model.overall.importance <- randomForest::importance(rf.model.overall, type=NULL, class=NULL, scale=TRUE)[,3:4]
	#randomForest::varImpPlot(rf.model.overall)
	
	#####################################

    if(verbose){
        now <- Sys.time()
        message(sprintf("5. Performance evaluation using test sets (%s).", as.character(now)), appendLF=TRUE)
        message(sprintf("Extract the ROC matrix of %d rows (Supervised + predictors) X %d columns/repeats*folds (%s).", ncol(df_predictor_class), nfold*nrepeat, as.character(now)), appendLF=TRUE)
    }
	
	######################
	## evaluation per fold
	######################
	lsls_predictors <- lapply(1:length(ls_model), function(i){
		rf.model <- ls_model[[i]]
		## prediction for testset: ?predict.randomForest
		trainindex <- index_sets[[i]]
		testset <- df_predictor_class[-trainindex,]
		vec_predict_test <- predict(rf.model, newdata=testset[,-ncol(testset)], type='prob')[,2]
		### do preparation
		ind <- match(rownames(testset), rownames(df_predictor))
		df_predictor_test <- cbind(Supervised_randomforest=as.numeric(vec_predict_test), df_predictor[ind,])
		rownames(df_predictor_test) <- rownames(df_predictor[ind,])
		df_pred <- df_predictor_test
		ls_predictors <- lapply(colnames(df_pred), function(x){
			data.frame(rownames(df_pred), df_pred[,x], stringsAsFactors=FALSE)
		})
		names(ls_predictors) <- colnames(df_pred)
		return(ls_predictors)
	})
	names(lsls_predictors) <- names(ls_model)
	ls_res <- lapply(1:length(lsls_predictors), function(i){
		ls_predictors <- lsls_predictors[[i]]
		# do evaluation
		ls_pPerf <- lapply(ls_predictors, function(x){
			pPerf <- xPredictROCR(prediction=x, GSP=GSP, GSN=GSN, verbose=FALSE)
		})
		# do plotting
		bp <- xPredictCompare(ls_pPerf, displayBy=c("ROC","PR"))
		
		if(is.null(bp)){
			df <- NULL
		}else{
			df <- unique(bp$data[,c('Method','auroc','fmax')])
			df <- data.frame(predictor=df$Method, model=rep(names(ls_model)[i],nrow(df)), ROC=df$auroc, Fmax=df$fmax, stringsAsFactors=FALSE)
		}
		
		return(df)
	})
	df_res <- do.call(rbind, ls_res)
	if(!is.null(df_res)){
		## df_ROC
		df_res <- as.matrix(xSparseMatrix(df_res[,-4], verbose=FALSE))
		ind <- match(c("Supervised_randomforest",colnames(df_predictor)), rownames(df_res))
		ind <- ind[!is.na(ind)]
		if(nfold==1 & nrepeat==1){
			df_res <- as.matrix(df_res[ind,], ncol=nfold)
			colnames(df_res) <- 'Fold1'
		}else{
			df_res <- df_res[ind,]
		}
		vec_median <- apply(df_res, 1, stats::median)
		vec_mad <- apply(df_res, 1, stats::mad)
		vec_min <- apply(df_res, 1, base::min)
		vec_max <- apply(df_res, 1, base::max)
		#####
		res_x <- apply(df_res, 1, stats::t.test)
		vec_mean <- unlist(lapply(res_x, function(x) x$estimate))
		vec_conf_lower <- unlist(lapply(res_x, function(x) x$conf.int[1]))
		vec_conf_upper <- unlist(lapply(res_x, function(x) x$conf.int[2]))
		#####
		df_ROC <- cbind(median=vec_median, mad=vec_mad, min=vec_min, max=vec_max, mean=vec_mean, conf_lower=vec_conf_lower, conf_upper=vec_conf_upper, df_res)
		
		## df_Fmax
		df_res <- do.call(rbind, ls_res)
		df_res <- as.matrix(xSparseMatrix(df_res[,-3], verbose=FALSE))
		ind <- match(c("Supervised_randomforest",colnames(df_predictor)), rownames(df_res))
		if(nfold==1 & nrepeat==1){
			df_res <- as.matrix(df_res[ind,], ncol=nfold)
			colnames(df_res) <- 'Fold1'
		}else{
			df_res <- df_res[ind,]
		}
		vec_median <- apply(df_res, 1, stats::median)
		vec_mad <- apply(df_res, 1, stats::mad)
		vec_min <- apply(df_res, 1, base::min)
		vec_max <- apply(df_res, 1, base::max)
		#####
		if(0){
		res_x <- apply(df_res, 1, stats::t.test)
		vec_mean <- unlist(lapply(res_x, function(x) x$estimate))
		vec_conf_lower <- unlist(lapply(res_x, function(x) x$conf.int[1]))
		vec_conf_upper <- unlist(lapply(res_x, function(x) x$conf.int[2]))
		#####
		df_Fmax <- cbind(median=vec_median, mad=vec_mad, min=vec_min, max=vec_max, mean=vec_mean, conf_lower=vec_conf_lower, conf_upper=vec_conf_upper, df_res)
		}else{
		df_Fmax <- cbind(median=vec_median, mad=vec_mad, min=vec_min, max=vec_max, df_res)
		}
    
    
    }else{
    	df_ROC <- NULL
    	df_Fmax <- NULL
    }
    
	#####################################
    if(verbose){
        now <- Sys.time()
        message(sprintf("6. Do prediction for fullset (%s).", as.character(now)), appendLF=TRUE)
        message(sprintf("Extract the full prediction matrix of %d rows/genes X %d columns/repeats*folds, aggregated via '%s' (%s) ...", nrow(df_predictor_class), nfold*nrepeat, fold.aggregateBy, as.character(now)), appendLF=TRUE)
    }
	
	ls_full <- lapply(1:length(ls_model), function(i){
		rf.model <- ls_model[[i]]
		## prediction for fullset: ?predict.randomForest
		vec_predict_full <- predict(rf.model, newdata=df_predictor, type='prob')[,2]
		# output
		df <- data.frame(genes=names(vec_predict_full), model=rep(names(ls_model)[i],length(vec_predict_full)), prob=vec_predict_full, stringsAsFactors=FALSE)
	})
	df_full <- do.call(rbind, ls_full)
	df_full <- as.matrix(xSparseMatrix(df_full, verbose=FALSE))

	## Convert into p-values by computing an empirical cumulative distribution function
    ls_pval <- lapply(1:ncol(df_full), function(j){
    	x <- df_full[,j]
    	my.CDF <- stats::ecdf(x)
    	pval <- 1 - my.CDF(x)
    })
	df_pval <- do.call(cbind, ls_pval)
	rownames(df_pval) <- rownames(df_full)
	
	## aggregate p values
	df_ap <- dnet::dPvalAggregate(pmatrix=df_pval, method=fold.aggregateBy)
	df_ap <- sort(df_ap, decreasing=FALSE)

	## get rank
	df_rank <- rank(df_ap, ties.method="min")
	#df_rank <- match(df_rank, sort(unique(df_rank)))

	######
	df_ap[df_ap==0] <- min(df_ap[df_ap!=0])
	######

	## adjp
	df_adjp <- stats::p.adjust(df_ap, method="BH")
	
	## rating: first log10-transformed ap and then being rescaled into the [0,5] range
	rating <- -log10(df_ap)
	####
	rating <- sqrt(rating)
	####
	rating <- 5 * (rating - min(rating))/(max(rating) - min(rating))
	
	#########################################
	## output
	### df_priority
	output_gs <- rep('NEW', length(df_ap))
	names(output_gs) <- names(df_ap)
	ind <- match(names(df_ap), names(gs_targets))
	output_gs[!is.na(ind)] <- gs_targets[ind[!is.na(ind)]]
	output_gs[output_gs=='0'] <- 'GSN'
	output_gs[output_gs=='1'] <- 'GSP'
	df_priority <- data.frame(GS=output_gs, name=names(df_ap), rank=df_rank, pvalue=df_ap, fdr=df_adjp, rating=rating, stringsAsFactors=FALSE)
	### add description
	df_priority$description <- xSymbol2GeneID(df_priority$name, details=TRUE, RData.location=RData.location, guid=guid)$description
	###
	
	### df_predictor_gs df_full_gs
	ind <- match(names(df_ap), rownames(df_predictor))
	if(nfold==1 & nrepeat==1){
		output_df_full <- as.matrix(df_full[ind,], ncol=nfold)
		colnames(output_df_full) <- 'Fold1'
	}else{
		output_df_full <- df_full[ind,]
	}
	output_df_predictor <- df_predictor[ind,]
	df_predictor_gs <- data.frame(GS=output_gs, name=names(df_ap), output_df_predictor, stringsAsFactors=FALSE)
	df_full_gs <- data.frame(GS=output_gs, name=names(df_ap), output_df_full, stringsAsFactors=FALSE)
	
	### pred2fold
	pred2fold <- lapply(1:length(lsls_predictors), function(i){
		x <- lsls_predictors[[i]][[1]]
		colnames(x) <- c("name","prob")
		return(x)
	})
	names(pred2fold) <- names(lsls_predictors)
    ####################################################################################

	######################
	## overall evaluation
	######################
	### do preparation
	df_predictor_overall <- cbind(Supervised_randomforest=df_priority$rating, df_predictor_gs[,-c(1,2)])
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
	
	if(!is.null(eTarget)){
		ind <- match(df_priority$name, rownames(eTarget$evidence))
		eTarget$evidence <- eTarget$evidence[ind,]
	}
	
	priority <- data.frame(df_priority[,c("GS","name","rank","rating","description")], stringsAsFactors=FALSE)
	
    sTarget <- list(ls_model = ls_model,
    				priority = priority,
    				predictor = df_predictor_gs,
    				pred2fold = pred2fold,
    				prob2fold = df_full_gs,
    				importance2fold = df_importance,
    				roc2fold = df_ROC,
    				fmax2fold = df_Fmax,
    				importance = rf.model.overall.importance,
    				performance = df_evaluation,
    				evidence = eTarget,
    				list_pNode  = list_pNode
                 )
    class(sTarget) <- "sTarget"
    
    ####################################################################################
    endT <- Sys.time()
    if(verbose){
        message(paste(c("\nFinish at ",as.character(endT)), collapse=""), appendLF=TRUE)
    }
    
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=TRUE)
    
    invisible(sTarget)
}


