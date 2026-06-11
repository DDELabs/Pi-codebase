#' Function to visualise cross-validation performance against tuning parameters
#'
#' \code{xMLparameters} is supposed to visualise cross-validation performance against tuning parameters. 
#'
#' @param data an object of the class "train" or "train.formula" (resulting from caret::train) used for 1D or 2D visualisation. Alternatively, it can be a data frame used for 2D or 3D visualisation
#' @param nD an integer specifying the dimension of the visualisation. It can be one of '1D', '2D' and '3D' and 'auto' (if input data is a "train" object)
#' @param contour logical to indicate whether coutour plot should be also included
#' @param main a title for the plot
#' @param xlab a title for the x axis
#' @param ylab a title for the y axis
#' @param zlab a title for the z axis
#' @param clab a title for the colorbar
#' @param nlevels the number of levels to partition the input matrix values. The same level has the same color mapped to
#' @param colormap short name for the colormap. It can be one of "jet" (jet colormap), "bwr" (blue-white-red colormap), "gbr" (green-black-red colormap), "wyr" (white-yellow-red colormap), "br" (black-red colormap), "yr" (yellow-red colormap), "wb" (white-black colormap), and "rainbow" (rainbow colormap, that is, red-yellow-green-cyan-blue-magenta). Alternatively, any hyphen-separated HTML color names, e.g. "blue-black-yellow", "royalblue-white-sandybrown", "darkgreen-white-darkviolet". A list of standard color names can be found in \url{http://html-color-codes.info/color-names}
#' @param highlight logical whether to highlight the point with the maximum value
#' @param x.label.cex the x-axis label size
#' @param x.label.srt the x-axis label angle (in degree from horizontal)
#' @param theta.3D the azimuthal direction. By default, it is 40
#' @param phi.3D the colatitude direction. By default, it is 20
#' @return 
#' invisible
#' @note none
#' @export
#' @seealso \code{\link{xSparseMatrix}}
#' @include xMLparameters.r
#' @examples
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' \dontrun{
#' xMLparameters(df_fit, nD="2D")
#' xMLparameters(df_fit, nD="3D", theta.3D=40, phi.3D=60)
#' }

xMLparameters <-function(data, nD=c("auto","1D","2D","3D"), contour=TRUE, main='Repeated cross-validation', xlab=NA, ylab=NA, zlab=NA, clab="AUC (repeated CV)", nlevels=50, colormap=c("lightblue-lightyellow-darkorange-darkred","bwr","jet","gbr","wyr","br","yr","rainbow","wb"), highlight=TRUE, x.label.cex=0.8, x.label.srt=30, 
theta.3D=40, phi.3D=25)
{

    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    nD <- match.arg(nD)

    if(is(data,"train")){
    	if(nD=="auto"){
			if(data$method %in% c("gbm","svmRadial","rda","nnet","glmnet","xgbTree","custom")){
				nD <- "2D"
			}else if(data$method %in% c("knn","pls","LogitBoost","xgbLinear")){
				nD <- "1D"
			}
    	}
    	
    	# ?caret::plot.train
    	if(nD=='2D'){
			res <- plot(data, metric="ROC", plotType="level", main=main, xlab=xlab, ylab=ylab, scales=list(x=list(rot=x.label.srt,cex=x.label.cex)), region=TRUE, pretty=TRUE, col.regions=xColormap(colormap)(nlevels), colorkey=list(tck=0.8))
			return(res)
			
        }else if(nD=='1D'){
        	res <- ggplot(data) + theme_bw() + labs(title=main) + theme(plot.title=element_text(hjust=0.5,face='bold'))
        	return(res)
        	
        }else if(nD=='3D'){
    		stop("The function only supports 1D and 2D for a 'train' or 'train.formula' object.\n")
    	}else{
    		return(NULL)
    	}
        
    }else if(is(data,"data.frame") & nD!="auto"){
		df_fit <- data
		## get parameter name
		if(is.na(xlab)){
			xlab <- colnames(df_fit)[1]
		}
		if(is.na(ylab)){
			ylab <- colnames(df_fit)[2]
		}
		if(is.na(zlab)){
			zlab <- colnames(df_fit)[3]
		}
		
    }else{
    	stop("The function must apply to a 'list' object, or a 'train'/'train.formula' object.\n")
    }

	if(nD=="2D"){
		## prepare data frame (x, y, z)
		x <- as.factor(df_fit[,1])
		y <- as.factor(df_fit[,2])
		z <- df_fit[,3]
		df_xyz <- data.frame(x=x, y=y, z=z)
		
		if(contour){
			res <- lattice::contourplot(z~x*y, df_xyz, xlab=xlab, ylab=ylab, main=main, scales=list(x=list(rot=x.label.srt,cex=x.label.cex)), region=TRUE, pretty=TRUE, col.regions=xColormap(colormap)(nlevels), colorkey=list(tck=0.8))
			
		}else{
			res <- lattice::levelplot(z~x*y, df_xyz, xlab=xlab, ylab=ylab, main=main, scales=list(x=list(rot=x.label.srt,cex=x.label.cex)), region=TRUE, pretty=TRUE, col.regions=xColormap(colormap)(nlevels), colorkey=list(tck=0.8))
			
		}
		
	}else if(nD=="3D"){
		
		## prepare data colvar matrix
		data_colvar <- as.matrix(xSparseMatrix(df_fit[,c(1:3)], verbose=FALSE))
	
		## prepare data matrix
		data <- as.matrix(xSparseMatrix(df_fit[,1:3], verbose=FALSE))
		rows <- as.character(unique(df_fit[,1]))
		columns <- as.character(unique(df_fit[,2]))
		ind_rows <- match(rows, rownames(data))
		ind_columns <- match(columns, colnames(data))
		data <- data[ind_rows, ind_columns]
		data_colvar <- data_colvar[ind_rows, ind_columns]
		
		## x- and y-axis
		x_at <- 1:nrow(data)
		y_at <- 1:ncol(data)
		#x_at <- as.numeric(rownames(data))
		#y_at <- as.numeric(colnames(data))
		
		# highlight
		hightlight_ind <- which(data_colvar == max(data_colvar), arr.ind=TRUE)
		hightlight_x <- x_at[hightlight_ind[, 1]]
		hightlight_y <- y_at[hightlight_ind[, 2]]
		hightlight_zval <- round(data[hightlight_ind], 3)
		hightlight_xval <- rownames(data)[hightlight_x]
		hightlight_yval <- colnames(data)[hightlight_y]
		
		if(highlight){
			plot <- FALSE
		}else{
			plot <- TRUE
		}
		
		## 3D plot
		zlim <- c(floor(min(data)*100)/100, ceiling(max(data)*100)/100)
		zlim <- c(1.5*zlim[1]-0.5*zlim[2], zlim[2])
		#zlim <- c(0.5, zlim[2])
		plot3D::persp3D(z=data, x=x_at, y=y_at, colvar=data_colvar, axes=TRUE, zlim=zlim, cex.axis=0.8, cex.lab=1.2, ticktype=c("simple","detailed")[1], col=xColormap(colormap)(nlevels), colkey=list(side=4,length=0.3,width=0.6,shift=0.2,cex.axis=0.6,cex.clab=0.8,side.clab=2,tck=-0.3), xlab=paste0("\n",xlab), ylab=paste0("\n",ylab), zlab=paste0("\n",zlab), clab=c("","",clab), bty="b", facets=TRUE, curtain=FALSE, plot=plot, lighting=FALSE, lphi=90, theta=theta.3D, phi=phi.3D, contour=list(col="grey",lwd=0.8))
		if(highlight){
			plot3D::scatter3D(x=hightlight_x, y=hightlight_y, z=zlim[1], type="h", colkey=FALSE, pch=19, cex=0.6, alpha=0.5, col="black", add=TRUE, plot=FALSE)
			label <- paste0(zlab,"=",hightlight_zval,"\n(",xlab,"=",hightlight_xval,")\n(",ylab,"=",hightlight_yval,")")
			plot3D::text3D(x=hightlight_x+0.5, y=hightlight_y, z=zlim[1], label=label, colkey=FALSE, cex=0.8, col="red", srt=0, add=TRUE, plot=TRUE)
		}
		
	}

    invisible()
}

    
