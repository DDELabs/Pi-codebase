#' Function to visualise a numeric matrix as a contour plot
#'
#' \code{xContour} is supposed to visualise a numeric matrix as a contour plot. 
#'
#' @param data a numeric matrix for the contour plot
#' @param main an overall title for the plot
#' @param xlab a title for the x axis. If specified, it will override 'names(dimnames(data))[1]'
#' @param ylab a title for the y axis. If specified, it will override 'names(dimnames(data))[2]'
#' @param key a title for the key plot (on the right)
#' @param nlevels the number of levels to partition the input matrix values. The same level has the same color mapped to
#' @param colormap short name for the colormap. It can be one of "jet" (jet colormap), "bwr" (blue-white-red colormap), "gbr" (green-black-red colormap), "wyr" (white-yellow-red colormap), "br" (black-red colormap), "yr" (yellow-red colormap), "wb" (white-black colormap), and "rainbow" (rainbow colormap, that is, red-yellow-green-cyan-blue-magenta). Alternatively, any hyphen-separated HTML color names, e.g. "blue-black-yellow", "royalblue-white-sandybrown", "darkgreen-white-darkviolet". A list of standard color names can be found in \url{http://html-color-codes.info/color-names}
#' @param highlight how to highlight the point. It can be 'none' for no highlight (by default), 'min' for highlighting the point with the minimum value of the matrix, and 'max' for highlighting the point with the maximum value of the matrix
#' @param highlight.col the highlight colors
#' @param x.label.cex the x-axis label size
#' @param x.label.srt the x-axis label angle (in degree from horizontal)
#' @param signature a logical to indicate whether the signature is assigned to the plot caption. By default, it sets FALSE
#' @param ... additional graphic parameters. For most parameters, please refer to \url{http://stat.ethz.ch/R-manual/R-devel/library/graphics/html/filled.contour.html}
#' @return 
#' invisible
#' @note none
#' @export
#' @seealso \code{\link{xContour}}
#' @include xContour.r
#' @examples
#' x <- y <- seq(-4*pi, 4*pi, len=10)
#' r <- sqrt(outer(x^2, y^2, "+"))
#' data <- cos(r^2)*exp(-r/(2*pi))
#' xContour(data)
#' #xContour(data, signature=TRUE)

xContour <-function(data, main='', xlab='', ylab='', key='', nlevels=50, colormap=c("darkblue-lightblue-lightyellow-darkorange","bwr","jet","gbr","wyr","br","yr","rainbow","wb"), highlight=c('none','min','max'), highlight.col="white", x.label.cex=0.95, x.label.srt=30, signature=FALSE, ...)
{

    highlight <- match.arg(highlight)
    
    palette.name <- supraHex::visColormap(colormap=colormap)
    
    if(is.null(rownames(data))){
    	rownames(data) <- as.integer(1:nrow(data))
    }
    if(is.null(colnames(data))){
    	colnames(data) <- as.integer(1:ncol(data))
    }
    
	x_at <- (1:dim(data)[1]-1)/(dim(data)[1]-1)
	y_at <- (1:dim(data)[2]-1)/(dim(data)[2]-1)
    
    if(xlab==''){
    	xlab <- names(dimnames(data))[1]
    }
    if(ylab==''){
    	ylab <- names(dimnames(data))[2]
    }
    
	#################################################################################
	.my.filled.contour <- function (x = seq(0, 1, length.out = nrow(z)), y = seq(0, 1, 
		length.out = ncol(z)), z, xlim = range(x, finite = TRUE), 
		ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE), 
		levels = pretty(zlim, nlevels), nlevels = 20, color.palette = cm.colors, 
		col = color.palette(length(levels) - 1), plot.title, plot.axes, 
		key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1, 
		axes = TRUE, frame.plot = axes, ...) 
	{
		if (missing(z)) {
			if (!missing(x)) {
				if (is.list(x)) {
					z <- x$z
					y <- x$y
					x <- x$x
				}
				else {
					z <- x
					x <- seq.int(0, 1, length.out = nrow(z))
				}
			}
			else stop("no 'z' matrix specified")
		}
		else if (is.list(x)) {
			y <- x$y
			x <- x$x
		}
		if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
			stop("increasing 'x' and 'y' values expected")
		mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
		on.exit(par(par.orig))
		#w <- (3 + mar.orig[2L]) * par("csi") * 2.54
		w <- (2.5 + mar.orig[2L]) * par("csi") * 2.54
		layout(matrix(c(2, 1), ncol = 2L), widths = c(1, lcm(w)))
		par(las = las)
		mar <- mar.orig
		mar[4L] <- mar[2L]
		mar[2L] <- 1
		par(mar = mar)
		plot.new()
		plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = "i", 
			yaxs = "i")
		#rect(0, levels[-length(levels)], 1, levels[-1L], col = col)
		rect(0, levels[-length(levels)], 1, levels[-1L], col = col, border = NA)
		if (missing(key.axes)) {
			if (axes) 
				axis(4)
		}
		else key.axes
		box()
		if (!missing(key.title)) 
			key.title
		mar <- mar.orig
		mar[4L] <- 1
		mar[2L] <- mar[2L]*1.5
		par(mar = mar)
		plot.new()
		plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)
		.filled.contour(x, y, z, levels, col)
		if (missing(plot.axes)) {
			if (axes) {
				title(main = "", xlab = "", ylab = "")
				Axis(x, side = 1)
				Axis(y, side = 2)
			}
		}
		else plot.axes
		if (frame.plot) 
			box()
		if (missing(plot.title)) 
			title(...)
		else plot.title
		invisible()
	}
	#################################################################################

	.my.filled.contour(data,
		plot.title=title(main=main, xlab=xlab, ylab=ylab),
		key.title=title(key),
		nlevels=nlevels,
		color.palette=palette.name,
		plot.axes={
			## x-axis
			#graphics::axis(side=1, at=x_at, labels=rownames(data), las=2)
			graphics::axis(side=1, at=x_at, labels=FALSE)
			graphics::text(x=x_at, y=par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]), labels=rownames(data), cex=x.label.cex, srt=x.label.srt, adj=0.5, xpd=TRUE)
			## y-axis
			graphics::axis(side=2, at=y_at, labels=colnames(data), las=2)
			## highlight
			if(highlight=='min'){
				hightlight_ind <- which(data==min(data), arr.ind=TRUE)
			}else if(highlight=='max'){
				hightlight_ind <- which(data==max(data), arr.ind=TRUE)
			}
			if(highlight!='none'){
				hightlight_x <- x_at[hightlight_ind[,1]]
				hightlight_y <- y_at[hightlight_ind[,2]]
				hightlight_z <- round(data[hightlight_ind],3)
				graphics::points(hightlight_x, hightlight_y, cex=0.8, pch=23, bg=highlight.col)
				graphics::text(hightlight_x, hightlight_y, hightlight_z, col=highlight.col, adj=c(-0.1,0))
			}
		},
		key.axes={
			graphics::axis(4)
		}
	)
    
    if(signature){
    	caption <- paste("Created by xContour from Pi version", utils::packageVersion("Pi"))
    	graphics::mtext(caption, side=1, line=4, adj=1, cex=.66, font=3)
    }
    
    invisible()
}

    
