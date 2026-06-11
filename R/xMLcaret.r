#' Function to integrate predictor matrix in a supervised manner via machine learning algorithms using caret.
#'
#' \code{xMLcaret} is supposed to integrate predictor matrix in a supervised manner via machine learning algorithms using caret. The caret package streamlines model building and performance evaluation. It requires three inputs: 1) Gold Standard Positive (GSP) targets; 2) Gold Standard Negative (GSN) targets; 3) a predictor matrix containing genes in rows and predictors in columns, with their predictive scores inside it. It returns an object of class 'sTarget'.
#'
#' @param list_pNode a list of "pNode" objects or a "pNode" object
#' @param df_predictor a data frame containing genes (in rows) and predictors (in columns), with their predictive scores inside it. This data frame must has gene symbols as row names
#' @param GSP a vector containing Gold Standard Positive (GSP)
#' @param GSN a vector containing Gold Standard Negative (GSN)
#' @param method machine learning method. It can be one of "gbm" for Gradient Boosting Machine (GBM), "svmRadial" for Support Vector Machines with Radial Basis Function Kernel (SVM), "rda" for Regularized Discriminant Analysis (RDA), "knn" for k-nearest neighbor (KNN), "pls" for Partial Least Squares (PLS), "nnet" for Neural Network (NNET), "rf" for Random Forest (RF), "myrf" for customised Random Forest (RF), "cforest" for Conditional Inference Random Forest, "glmnet" for glmnet, "glm" for Generalized Linear Model (GLM), "bayesglm" for Bayesian Generalized Linear Model (BGLM), "LogitBoost" for Boosted Logistic Regression (BLR), "xgbLinear" for eXtreme Gradient Boosting as linear booster (XGBL), "xgbTree" for eXtreme Gradient Boosting as tree booster (XGBT)
#' @param nfold an integer specifying the number of folds for cross validataion. Per fold creates balanced splits of the data preserving the overall distribution for each class (GSP and GSN), therefore generating balanced cross-vallidation train sets and testing sets. By default, it is 3 meaning 3-fold cross validation
#' @param nrepeat an integer specifying the number of repeats for cross validataion. By default, it is 10 indicating the cross-validation repeated 10 times
#' @param seed an integer specifying the seed
#' @param aggregateBy the aggregate method used to aggregate results from repeated cross validataion. It can be either "none" for no aggregration (meaning the best model based on all data used for cross validation is used), or "orderStatistic" for the method based on the order statistics of p-values, or "fishers" for Fisher's method, "Ztransform" for Z-transform method, "logistic" for the logistic method. Without loss of generality, the Z-transform method does well in problems where evidence against the combined null is spread widely (equal footings) or when the total evidence is weak; Fisher's method does best in problems where the evidence is concentrated in a relatively small fraction of the individual tests or when the evidence is at least moderately strong; the logistic method provides a compromise between these two. Notably, the aggregate methods 'Ztransform' and 'logistic' are preferred here
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to TRUE for display
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @param guid a valid (5-character) Global Unique IDentifier for an OSF project. See \code{\link{xRDataLoader}} for details
#' @return 
#' an object of class "sTarget", a list with following components:
#' \itemize{
#'  \item{\code{model}: an object of class "train" as a best model}
#'  \item{\code{ls_model}: a list of best models from repeated cross-validation}
#'  \item{\code{priority}: a data frame of n X 5 containing gene priority information, where n is the number of genes in the input data frame, and the 5 columns are "GS" (either 'GSP', or 'GSN', or 'Putative'), "name" (gene names), "rank" (priority rank), "rating" (5-star priority score/rating), and "description" (gene description)}
#'  \item{\code{predictor}: a data frame, which is the same as the input data frame but inserting two additional columns ('GS' in the first column, 'name' in the second column)}
#'  \item{\code{performance}: a data frame of 1+nPredictor X 2 containing the supervised/predictor performance info, where nPredictor is the number of predictors, two columns are "ROC" (AUC values) and "Fmax" (F-max values)}
#'  \item{\code{performance_cv}: a data frame of nfold*nrepeat X 2 containing the repeated cross-validation performance, where two columns are "ROC" (AUC values) and "Fmax" (F-max values)}
#'  \item{\code{importance}: a data frame of nPredictor X 1 containing the predictor importance info}
#'  \item{\code{gp}: a ggplot object for the ROC curve}
#'  \item{\code{gp_cv}: a ggplot object for the ROC curves from repeated cross-validation}
#'  \item{\code{evidence}: an object of the class "eTarget", a list with following components "evidence" and "metag"}
#'  \item{\code{list_pNode}: a list of "pNode" objects}
#' }
#' @note It will depend on whether a package "caret" and its suggested packages have been installed. It can be installed via: \code{BiocManager::install(c("caret","e1071","gbm","kernlab","klaR","pls","nnet","randomForest","party","glmnet","arm","caTools","xgboost"))}.
#' @export
#' @seealso \code{\link{xPierMatrix}}, \code{\link{xPredictROCR}}, \code{\link{xPredictCompare}}, \code{\link{xSparseMatrix}}, \code{\link{xSymbol2GeneID}}
#' @include xMLcaret.r
#' @examples
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' \dontrun{
#' sTarget <- xMLcaret(df_prediction, GSP, GSN, method="myrf")
#' }

xMLcaret <- function(list_pNode=NULL, df_predictor=NULL, GSP, GSN, method=c("gbm","svmRadial","rda","knn","pls","nnet","rf","myrf","cforest","glmnet","glm","bayesglm","LogitBoost","xgbLinear","xgbTree"), nfold=3, nrepeat=10, seed=825, aggregateBy=c("none","logistic","Ztransform","fishers","orderStatistic"), verbose=TRUE, RData.location="http://galahad.well.ox.ac.uk/bigdata", guid=NULL)
{
	
    startT <- Sys.time()
    if(verbose){
        message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=TRUE)
        message("", appendLF=TRUE)
    }
    ####################################################################################
	
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    method <- match.arg(method)
    aggregateBy <- match.arg(aggregateBy)
	
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
	levels(class) <- c("GSN","GSP")
	df_predictor_class$class <- class
	
    if(verbose){
        now <- Sys.time()
        message(sprintf("1. Gold standards (%d in GSP, %d in GSN) are used for supervised integration of %d predictors/features (%s).", sum(class=="GSP"), sum(class=="GSN"), ncol(df_predictor), as.character(now)), appendLF=TRUE)
    }
	
    if(verbose){
        now <- Sys.time()
        message(sprintf("2. Optimise tuning parameters of machine learning algorithm '%s' using %d-fold cross validation repeated %d times (%s) ....", method, nfold, nrepeat, as.character(now)), appendLF=TRUE)
    }
	
	fitControl <- caret::trainControl(method=c("repeatedcv","cv","oob")[1], number=nfold, repeats=nrepeat, classProbs=TRUE, summaryFunction=caret::twoClassSummary, allowParallel=FALSE)
	fitControl_withoutParameters <- caret::trainControl(method="none", classProbs=TRUE, allowParallel=FALSE)
	
    #http://topepo.github.io/caret/available-models.html
    #http://topepo.github.io/caret/train-models-by-tag.html
	if(method=="gbm"){
		# Gradient Boosting Machine (GBM) model: 
		# specify tuning parameter grid
		#library(gbm)
		#caret::getModelInfo("gbm")$gbm$grid
		grid_gbm <-  base::expand.grid(
						n.trees = (1:10)*50, 
						interaction.depth = c(1:10),
						shrinkage = 0.1,
						n.minobsinnode = 10
						)
		if(!is.null(seed)) set.seed(seed)
		fit_gbm <- caret::train(class ~ ., 
								data = df_predictor_class, 
								method = c("gbm"), 
								trControl = fitControl,
								tuneGrid = grid_gbm,
								#tuneLength = 10,
								verbose = FALSE,
								metric = "ROC"
								)
		fit_target <- fit_gbm
		df_fit <- fit_target$results[,c(4,2,5,6,7)]
		
		# get prob function
		func_prob <- caret::getModelInfo("gbm")$gbm$prob
		
		# get a list of models from repeated cross-validation
		ls_model <- lapply(fit_target$control$index, function(k){
			x <- fit_target$trainingData[k, -1]
			y <- fit_target$trainingData[k, 1]
			fit <- caret::train(x=x, y=y, 
								method = c("gbm"), 
								trControl = fitControl_withoutParameters,
								tuneGrid = fit_target$bestTune,
								verbose = FALSE,
								metric = "ROC"
								)
			fit$finalModel
		})
		
	}else if(method=="svmRadial"){
		# Support Vector Machines with Radial Basis Function Kernel
		# specify tuning parameter grid
		#library(kernlab)
		#caret::getModelInfo("svmRadial")$svmRadial$grid
		grid_svm <-  base::expand.grid(
						C = 2^(-2:8), 
						sigma = c(1:10)*10
						)
		if(!is.null(seed)) set.seed(seed)
		fit_svm <- caret::train(class ~ ., 
								data = df_predictor_class, 
								method = "svmRadial", 
								trControl = fitControl,
								#tuneLength = 10,
								tuneGrid = grid_svm,
								metric = "ROC"
								)
		fit_target <- fit_svm
		df_fit <- fit_target$results[,1:5]
		
		# get prob function
		func_prob <- caret::getModelInfo("svmRadial")$svmRadial$prob
		
		# get a list of models from repeated cross-validation
		ls_model <- lapply(fit_target$control$index, function(k){
			x <- fit_target$trainingData[k, -1]
			y <- fit_target$trainingData[k, 1]
			fit <- caret::train(x=x, y=y, 
								method = c("svmRadial"), 
								trControl = fitControl_withoutParameters,
								tuneGrid = fit_target$bestTune,
								metric = "ROC"
								)
			fit$finalModel
		})
		
	}else if(method=="rda"){
		# Regularized Discriminant Analysis
		# specify tuning parameter grid
		#library(klaR)
		#caret::getModelInfo("rda")$rda$grid
		grid_rda <-  base::expand.grid(
						gamma = seq(0,1,0.1), 
						lambda = seq(0,1,0.1)
						)
		if(!is.null(seed)) set.seed(seed)
		fit_rda <- caret::train(class ~ ., 
								data = df_predictor_class, 
								method = "rda", 
								trControl = fitControl,
								#tuneLength = 10,
								tuneGrid = grid_rda,
								metric = "ROC"
								)
		fit_target <- fit_rda 
		df_fit <- fit_target$results[,c(2,1,3,4,5)]
		
		# get prob function
		func_prob <- caret::getModelInfo("rda")$rda$prob
		#func_fit <- caret::getModelInfo("rda")$rda$fit
		
		# get a list of models from repeated cross-validation
		ls_model <- lapply(fit_target$control$index, function(k){
			x <- fit_target$trainingData[k, -1]
			y <- fit_target$trainingData[k, 1]
			fit <- caret::train(x=x, y=y, 
								method = c("rda"), 
								trControl = fitControl_withoutParameters,
								tuneGrid = fit_target$bestTune,
								metric = "ROC"
								)
			fit$finalModel
		})
		
	}else if(method=="knn"){   			
		# KNN
		# specify tuning parameter grid
		#caret::getModelInfo("knn")$knn$grid
		grid_knn <-  base::expand.grid(
						k = seq(1,20,1)
						)
		if(!is.null(seed)) set.seed(seed)
		fit_knn <- caret::train(class ~ ., 
								data = df_predictor_class, 
								method = "knn", 
								trControl = fitControl,
								#tuneLength = 10,
								tuneGrid = grid_knn,
								metric = "ROC"
								)
		fit_target <- fit_knn    
		df_fit <- fit_target$results[,c(1,2,3,4)]
		
		# get prob function
		func_prob <- caret::getModelInfo("knn")$knn$prob
		
		# get a list of models from repeated cross-validation
		ls_model <- lapply(fit_target$control$index, function(k){
			x <- fit_target$trainingData[k, -1]
			y <- fit_target$trainingData[k, 1]
			fit <- caret::train(x=x, y=y, 
								method = c("knn"), 
								trControl = fitControl_withoutParameters,
								tuneGrid = fit_target$bestTune,
								metric = "ROC"
								)
			fit$finalModel
		})
		
    }else if(method=="pls"){
		# Partial Least Squares
		# specify tuning parameter grid
		#library(pls)
		#caret::getModelInfo("pls")$pls$grid
		grid_pls <-  base::expand.grid(
						ncomp = seq(1,20,1)
						)
		if(!is.null(seed)) set.seed(seed)
		fit_pls <- caret::train(class ~ ., 
								data = df_predictor_class, 
								method = "pls", 
								trControl = fitControl,
								#tuneLength = 10,
								tuneGrid = grid_pls,
								metric = "ROC"
								)
		fit_target <- fit_pls
		df_fit <- fit_target$results[,c(1,2,3,4)]
		
		# get prob function
		func_prob <- caret::getModelInfo("pls")$pls$prob
		
		# get a list of models from repeated cross-validation
		ls_model <- lapply(fit_target$control$index, function(k){
			x <- fit_target$trainingData[k, -1]
			y <- fit_target$trainingData[k, 1]
			fit <- caret::train(x=x, y=y, 
								method = c("pls"), 
								trControl = fitControl_withoutParameters,
								tuneGrid = fit_target$bestTune,
								metric = "ROC"
								)
			fit$finalModel
		})
		
    }else if(method=="nnet"){
		# Neural Network
		# specify tuning parameter grid
		#library(nnet)
		#caret::getModelInfo("nnet")$nnet$grid
		len <- 10
		grid_nnet <-  base::expand.grid(
						size = ((1:len)*2) - 1,
						decay = c(0, 10^seq(-1,-4,length=len-1))
						)
		if(!is.null(seed)) set.seed(seed)
		suppressMessages(fit_nnet <- caret::train(class ~ ., 
								data = df_predictor_class, 
								method = "nnet", 
								trControl = fitControl,
								tuneLength = 10,
								#tuneGrid = grid_nnet,
								trace = FALSE,
								metric = "ROC"
								)
							)
		fit_nnet$results$decay <- signif(fit_nnet$results$decay, digits=2)
		fit_target <- fit_nnet
		df_fit <- fit_target$results[,c(1,2,3,4,5)]

		# get prob function
		func_prob <- caret::getModelInfo("nnet")$nnet$prob
		
		# get a list of models from repeated cross-validation
		ls_model <- lapply(fit_target$control$index, function(k){
			x <- fit_target$trainingData[k, -1]
			y <- fit_target$trainingData[k, 1]
			suppressMessages(fit <- caret::train(x=x, y=y, 
								method = c("nnet"), 
								trControl = fitControl_withoutParameters,
								tuneGrid = fit_target$bestTune,
								trace = FALSE,
								metric = "ROC"
								)
							)
			fit$finalModel
		})
    
    }else if(method=="rf"){		    
		# Random Forest
		# specify tuning parameter grid
		#library(randomForest)
		#caret::getModelInfo("rf")$rf$grid
		grid_rf <-  base::expand.grid(
						mtry = seq(2,10,1)
						)
		if(!is.null(seed)) set.seed(seed)
		fit_rf <- caret::train(class ~ ., 
								data = df_predictor_class, 
								method = "rf", 
								trControl = fitControl,
								#tuneLength = 10,
								tuneGrid = grid_rf,
								ntree = 1000,
								importance = TRUE,
								metric = "ROC"
								)
		fit_target <- fit_rf   
		df_fit <- fit_target$results[,c(1,2,3,4)]
    
		# get prob function
		func_prob <- caret::getModelInfo("rf")$rf$prob
		
		# get a list of models from repeated cross-validation
		ls_model <- lapply(fit_target$control$index, function(k){
			x <- fit_target$trainingData[k, -1]
			y <- fit_target$trainingData[k, 1]
			fit <- caret::train(x=x, y=y, 
								method = c("rf"), 
								trControl = fitControl_withoutParameters,
								tuneGrid = fit_target$bestTune,
								metric = "ROC"
								)
			fit$finalModel
		})

    
    }else if(method=="myrf"){
    
		##############################
		my_rf <- list(label = "my_RF",
               library = "randomForest",
               type = "Classification",
               ## Tune over both parameters at the same time
               parameters = data.frame(parameter = c('ntree', 'mtry'),
                                       class = c("numeric", 'numeric'),
                                       label = c('#Number of Trees', 
                                                 '#Randomly Selected Predictors')),
               grid = function(x, y, len = NULL, search = "grid") {
                	if(search == "grid") {
                   		grid <- expand.grid(ntree = seq(100, 100*len, by = 100),
                                       mtry = caret::var_seq(p=ncol(x),len=len))
                   	}
               },
               loop = NULL,
               fit = function(x, y, wts, param, lev, last, classProbs, ...) { 
            		randomForest::randomForest(x, y, mtry = param$mtry, ntree = param$ntree, ...)
               },
               predict = function(modelFit, newdata, submodels = NULL) {  
               		if(!is.null(newdata)) predict(modelFit, newdata) else predict(modelFit)
               },
               prob = function(modelFit, newdata, submodels = NULL){
               		if(!is.null(newdata)) predict(modelFit, newdata, type = "prob") else predict(modelFit, type = "prob")
               },
               varImp = function(object, ...){
                    varImp <- randomForest::importance(object, ...)
                    if(object$type == "regression")
                      varImp <- data.frame(Overall = varImp[,"%IncMSE"])
                    else {
                      retainNames <- levels(object$y)
                      if(all(retainNames %in% colnames(varImp))) {
                        varImp <- varImp[, retainNames]
                      } else {
                        varImp <- data.frame(Overall = varImp[,1])
                      }
                    }
                    
                    out <- as.data.frame(varImp)
                    if(dim(out)[2] == 2) {
                      tmp <- apply(out, 1, mean)
                      out[,1] <- out[,2] <- tmp  
                    }
                    out
                },
               predictors = function(x, ...) {
                    ## After doing some testing, it looks like randomForest
                    ## will only try to split on plain main effects (instead
                    ## of interactions or terms like I(x^2).
                    varIndex <- as.numeric(names(table(x$forest$bestvar)))
                    varIndex <- varIndex[varIndex > 0]
                    varsUsed <- names(x$forest$ncat)[varIndex]
                    varsUsed
                  },
               levels = function(x) x$classes,
               sort = function(x) x[order(x[,1]),]
            )
		##############################
    	
		# my customised Random Forest
		# specify tuning parameter grid
		#library(randomForest)
		#caret::getModelInfo("rf")$rf$grid
		grid_myrf <-  base::expand.grid(
						ntree = seq(400,1600,200),
						mtry = seq(3,10,1)
						)
		if(!is.null(seed)) set.seed(seed)
		fit_myrf <- caret::train(class ~ ., 
								data = df_predictor_class, 
								method = my_rf, 
								trControl = fitControl,
								#tuneLength = 10,
								#ntree = 1000,
								tuneGrid = grid_myrf,
								importance = TRUE,
								metric = "ROC"
								)
		fit_target <- fit_myrf
		df_fit <- fit_target$results[,c(1,2,3,4,5)]
		#res <- xMLparameters(df_fit, nD="2D", contour=FALSE)
		#xMLparameters(df_fit, nD="3D", theta.3D=40, phi.3D=20)
		
		# get prob function
		func_prob <- my_rf$prob
		
		# get a list of models from repeated cross-validation
		ls_model <- lapply(fit_target$control$index, function(k){
			x <- fit_target$trainingData[k, -1]
			y <- fit_target$trainingData[k, 1]
			fit <- caret::train(x=x, y=y, 
								method = my_rf, 
								trControl = fitControl_withoutParameters,
								tuneGrid = fit_target$bestTune,
								metric = "ROC"
								)
			fit$finalModel
		})
		
    }else if(method=="cforest"){
		# Conditional Inference Random Forest
		# specify tuning parameter grid
		#library(party)
		#caret::getModelInfo("cforest")$cforest$grid
		grid_crf <-  base::expand.grid(
						mtry = seq(2,30,2)
						)
		if(!is.null(seed)) set.seed(seed)
		fit_crf <- caret::train(class ~ ., 
								data = df_predictor_class, 
								method = "cforest", 
								trControl = fitControl,
								#tuneLength = 10,
								tuneGrid = grid_crf,
								metric = "ROC"
								)
		fit_target <- fit_crf    
		df_fit <- fit_target$results[,c(1,2,3,4)]
    
		# get prob function
		func_prob <- caret::getModelInfo("cforest")$cforest$prob
		
		# get a list of models from repeated cross-validation
		ls_model <- lapply(fit_target$control$index, function(k){
			x <- fit_target$trainingData[k, -1]
			y <- fit_target$trainingData[k, 1]
			fit <- caret::train(x=x, y=y, 
								method = c("cforest"), 
								trControl = fitControl_withoutParameters,
								tuneGrid = fit_target$bestTune,
								metric = "ROC"
								)
			fit$finalModel
		})
    
    }else if(method=="glmnet"){
		# glmnet
		# specify tuning parameter grid
		#library(glmnet)
		#caret::getModelInfo("glmnet")$glmnet$grid
		len <- 5
		grid_glmnet <-  base::expand.grid(
						alpha = seq(0,1,length=len+1),
						lambda = 10^(-5:0)
						)
		if(!is.null(seed)) set.seed(seed)
		fit_glmnet <- caret::train(class ~ ., 
								data = df_predictor_class, 
								method = "glmnet", 
								trControl = fitControl,
								tuneLength = 10,
								#tuneGrid = grid_glmnet,
								metric = "ROC"
								)
		fit_glmnet$results$lambda <- signif(fit_glmnet$results$lambda, digits=2)
		fit_target <- fit_glmnet   
		df_fit <- fit_target$results[,1:5]
		
		# get prob function
		func_prob <- caret::getModelInfo("glmnet")$glmnet$prob

		# get a list of models from repeated cross-validation
		ls_model <- lapply(fit_target$control$index, function(k){
			x <- fit_target$trainingData[k, -1]
			y <- fit_target$trainingData[k, 1]
			fit <- caret::train(x=x, y=y, 
								method = c("glmnet"), 
								trControl = fitControl_withoutParameters,
								tuneGrid = fit_target$bestTune,
								metric = "ROC"
								)
			fit$finalModel
		})
		
    }else if(method=="glm"){
		# Generalized Linear Model (GLM)
		#caret::getModelInfo("glm")$glm$grid
		if(!is.null(seed)) set.seed(seed)
		fit_glm <- caret::train(class ~ ., 
								data = df_predictor_class, 
								method = c("glm"), 
								family = "binomial",
								trControl = fitControl,
								tuneLength = 10,
								metric = "ROC"
								)
		fit_target <- fit_glm
		df_fit <- fit_target$results[,1:3]
		
		# get prob function
		func_prob <- caret::getModelInfo("glm")$glm$prob

		# get a list of models from repeated cross-validation
		ls_model <- lapply(fit_target$control$index, function(k){
			x <- fit_target$trainingData[k, -1]
			y <- fit_target$trainingData[k, 1]
			fit <- caret::train(x=x, y=y, 
								method = c("glm"), 
								family = "binomial",
								trControl = fitControl_withoutParameters,
								tuneGrid = fit_target$bestTune,
								metric = "ROC"
								)
			fit$finalModel
		})
		
    }else if(method=="bayesglm"){
		# Bayesian Generalized Linear Model (BGLM)
		#library(arm)
		#caret::getModelInfo("bayesglm")$bayesglm$grid
		if(!is.null(seed)) set.seed(seed)
		fit_bglm <- caret::train(class ~ ., 
								data = df_predictor_class, 
								method = c("bayesglm"), 
								family = "binomial",
								trControl = fitControl,
								tuneLength = 10,
								metric = "ROC"
								)
		fit_target <- fit_bglm
		df_fit <- fit_target$results[,1:3]
		
		# get prob function
		func_prob <- caret::getModelInfo("bayesglm")$bayesglm$prob

		# get a list of models from repeated cross-validation
		ls_model <- lapply(fit_target$control$index, function(k){
			x <- fit_target$trainingData[k, -1]
			y <- fit_target$trainingData[k, 1]
			fit <- caret::train(x=x, y=y, 
								method = c("bayesglm"),
								family = "binomial",
								trControl = fitControl_withoutParameters,
								tuneGrid = fit_target$bestTune,
								metric = "ROC"
								)
			fit$finalModel
		})
		
    }else if(method=="LogitBoost"){
		# Boosted Logistic Regression (BLR)
		# specify tuning parameter grid
		#library(caTools)
		#caret::getModelInfo("LogitBoost")$LogitBoost$grid
		len <- 30
		grid_blr <-  base::expand.grid(
						nIter = 1 + (1:len)*10
						)
		if(!is.null(seed)) set.seed(seed)
		fit_blr <- caret::train(class ~ ., 
								data = df_predictor_class, 
								method = c("LogitBoost"),
								trControl = fitControl,
								tuneGrid = grid_blr,
								#tuneLength = 30,
								metric = "ROC"
								)
		fit_target <- fit_blr
		df_fit <- fit_target$results[,1:3]
		
		# get prob function
		func_prob <- caret::getModelInfo("LogitBoost")$LogitBoost$prob

		# get a list of models from repeated cross-validation
		ls_model <- lapply(fit_target$control$index, function(k){
			x <- fit_target$trainingData[k, -1]
			y <- fit_target$trainingData[k, 1]
			fit <- caret::train(x=x, y=y, 
								method = c("LogitBoost"), 
								trControl = fitControl_withoutParameters,
								tuneGrid = fit_target$bestTune,
								metric = "ROC"
								)
			fit$finalModel
		})
		
    }else if(method=="xgbLinear"){
		# eXtreme Gradient Boosting linear (XGBlinear)
		# specify tuning parameter grid
		#library(xgboost)
		#caret::getModelInfo("xgbLinear")$xgbLinear$grid
		len <- 10
		grid_xgbl <- base::expand.grid(
						  	lambda = 0,
							alpha = 0,
							nrounds = floor((1:len)*50),
							#alpha = c(0, 10^seq(-1,-4,length=len-1)),
							#nrounds = floor((1:len)*50),
							eta = 0.3
						  )
		if(!is.null(seed)) set.seed(seed)
		fit_xgbl <- caret::train(class ~ ., 
								data = df_predictor_class, 
								method = c("xgbLinear"),
								trControl = fitControl,
								tuneGrid = grid_xgbl,
								#tuneLength = 4,
								metric = "ROC"
								)
		fit_target <- fit_xgbl
		df_fit <- fit_target$results[,c(2:3,5:7)]
		
		# get prob function
		func_prob <- caret::getModelInfo("xgbLinear")$xgbLinear$prob

		# get a list of models from repeated cross-validation
		ls_model <- lapply(fit_target$control$index, function(k){
			x <- fit_target$trainingData[k, -1]
			y <- fit_target$trainingData[k, 1]
			fit <- caret::train(x=x, y=y, 
								method = c("xgbLinear"), 
								trControl = fitControl_withoutParameters,
								tuneGrid = fit_target$bestTune,
								metric = "ROC"
								)
			fit$finalModel
		})
		
    }else if(method=="xgbTree"){
		# eXtreme Gradient Boosting tree (XGBtree)
		# specify tuning parameter grid
		#library(xgboost)
		#caret::getModelInfo("xgbTree")$xgbTree$grid
		len <- 10
		grid_xgbt <- base::expand.grid(
							max_depth = seq(1, len),
                        	nrounds = floor((1:len) * 50),
                        	eta = 0.3,
                        	gamma = 0,
                        	colsample_bytree = 1,
                        	min_child_weight = 1,
                        	subsample = 0.75
						)
		if(!is.null(seed)) set.seed(seed)
		fit_xgbt <- caret::train(class ~ ., 
								data = df_predictor_class, 
								method = c("xgbTree"),
								trControl = fitControl,
								tuneGrid = grid_xgbt,
								#tuneLength = 5,
								metric = "ROC"
								)
		fit_target <- fit_xgbt
		df_fit <- fit_target$results[,c(2,7,8:10)]
		#res <- xMLparameters(df_fit, nD="2D", contour=FALSE)
		#xMLparameters(df_fit, nD="3D", theta.3D=40, phi.3D=25)
		
		# get prob function
		func_prob <- caret::getModelInfo("xgbTree")$xgbTree$prob

		# get a list of models from repeated cross-validation
		ls_model <- lapply(fit_target$control$index, function(k){
			x <- fit_target$trainingData[k, -1]
			y <- fit_target$trainingData[k, 1]
			fit <- caret::train(x=x, y=y, 
								method = c("xgbTree"), 
								trControl = fitControl_withoutParameters,
								tuneGrid = fit_target$bestTune,
								metric = "ROC"
								)
			fit$finalModel
		})
		
    }
	
	#####################################
    if(verbose){
        now <- Sys.time()
        message(sprintf("3. Performance evaluation using test sets (%s).", as.character(now)), appendLF=TRUE)
        message(sprintf("Extract the performance matrix of %d rows/repeats*folds X 2 (AUC and F-max) (%s).", nfold*nrepeat, as.character(now)), appendLF=TRUE)
    }
	
	######################
	## evaluation per fold*repeat
	######################
	ls_predictors <- lapply(1:length(ls_model), function(i){
		## prediction for testset
		newdata <- fit_target$trainingData[fit_target$control$indexOut[[i]], -1]
		df_predict_test <- func_prob(ls_model[[i]], newdata=newdata)
		rownames(df_predict_test) <- rownames(newdata)
		### do preparation
		data.frame(name=rownames(df_predict_test), score=as.numeric(df_predict_test[,"GSP"]), stringsAsFactors=FALSE)
	})
	names(ls_predictors) <- names(ls_model)
	# do evaluation
	ls_pPerf <- lapply(ls_predictors, function(x){
		pPerf <- xPredictROCR(prediction=x, GSP=GSP, GSN=GSN, verbose=FALSE)
	})
	# do plotting
	gp_cv <- xPredictCompare(ls_pPerf, displayBy=c("ROC","PR"), facet=FALSE, signature=FALSE)
	# extracting
	if(is.null(gp_cv)){
		df_ROC_Fmax <- NULL
	}else{
		df <- unique(gp_cv$data[,c('Method','auroc','fmax')])
		df_ROC_Fmax <- data.frame(predictor=df$Method, ROC=df$auroc, Fmax=df$fmax, stringsAsFactors=FALSE)
	}
    
	######################
	## evaluation using all training data
	######################
	newdata <- fit_target$trainingData[, -1]
	df_predict_test <- func_prob(fit_target$finalModel, newdata=newdata)
	rownames(df_predict_test) <- rownames(newdata)
	x <- data.frame(name=rownames(df_predict_test), score=as.numeric(df_predict_test[,"GSP"]), stringsAsFactors=FALSE)
	pPerf <- xPredictROCR(prediction=x, GSP=GSP, GSN=GSN, verbose=FALSE)
	ls_pPerf <- list(pPerf)
	gp <- xPredictCompare(ls_pPerf, displayBy="ROC", facet=TRUE, signature=FALSE)
	#####################################
	
	#####################################
    if(verbose){
        now <- Sys.time()
        message(sprintf("4. Do prediction for full set (%s).", as.character(now)), appendLF=TRUE)
        message(sprintf("Extract the full prediction matrix of %d rows/genes X %d columns/repeats*folds, aggregated via '%s' (%s) ...", nrow(df_predictor_class), nfold*nrepeat, aggregateBy, as.character(now)), appendLF=TRUE)
    }
	
	if(aggregateBy=="none"){
		## prediction for fullset
		df_full <- func_prob(fit_target$finalModel, newdata=df_predictor)
		vec_full <- df_full[,"GSP"]
		names(vec_full) <- rownames(df_predictor)
		vec_full <- sort(vec_full, decreasing=TRUE)
	
		## get rank
		vec_rank <- rank(-1*vec_full, ties.method="min")
	
		## priority: being rescaled into the [0,5] range
		priority <- vec_full
		vec_priority <- 5 * (priority - min(priority))/(max(priority) - min(priority))
		
	}else{

		ls_full <- lapply(1:length(ls_model), function(i){
			## prediction for fullset
			df_predict_full <- func_prob(ls_model[[i]], newdata=df_predictor)
			# output
			df <- data.frame(genes=rownames(df_predictor), model=rep(names(ls_model)[i],nrow(df_predict_full)), prob=df_predict_full[,"GSP"], stringsAsFactors=FALSE)
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
		vec_ap <- dnet::dPvalAggregate(pmatrix=df_pval, method=aggregateBy)
		vec_ap <- sort(vec_ap, decreasing=FALSE)

		## get rank
		vec_rank <- rank(vec_ap, ties.method="min")

		######
		vec_ap[vec_ap==0] <- min(vec_ap[vec_ap!=0])
		######

		## adjp
		vec_adjp <- stats::p.adjust(vec_ap, method="BH")
	
		## rating: first log10-transformed ap and then being rescaled into the [0,5] range
		rating <- -log10(vec_ap)
		####
		rating <- sqrt(rating)
		####
		vec_priority <- 5 * (rating - min(rating))/(max(rating) - min(rating))

	}
	
	#########################################
	## output
	### df_priority
	output_gs <- rep('NEW', length(vec_priority))
	names(output_gs) <- names(vec_priority)
	ind <- match(names(vec_priority), names(gs_targets))
	output_gs[!is.na(ind)] <- gs_targets[ind[!is.na(ind)]]
	output_gs[output_gs=='0'] <- 'GSN'
	output_gs[output_gs=='1'] <- 'GSP'
	df_priority <- data.frame(GS=output_gs, name=names(vec_priority), rank=vec_rank, rating=vec_priority, stringsAsFactors=FALSE)
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
	df_predictor_overall <- cbind(Supervised=df_priority$rating, df_predictor_gs[,-c(1,2)])
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
	
	######################
	## overall importance
	######################
	overall.importance <- caret::varImp(fit_target)$importance
	#####################
	
	if(!is.null(eTarget)){
		ind <- match(df_priority$name, rownames(eTarget$evidence))
		eTarget$evidence <- eTarget$evidence[ind,]
	}
	
    sTarget <- list(model = fit_target, 
    				ls_model = ls_model,
    				priority = df_priority,
    				predictor = df_predictor_gs,
    				performance = df_evaluation,
    				gp = gp,
    				performance_cv = df_ROC_Fmax,
    				gp_cv = gp_cv,
    				importance = overall.importance,
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


