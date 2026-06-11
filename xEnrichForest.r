xEnrichForest <- function(eTerm, top_num=10, FDR.cutoff=0.05, CI.one=TRUE, colormap="ggplot2.top", ncolors=64, zlim=NULL, barwidth=0.5, barheight=NULL, wrap.width=NULL, font.family="sans", drop=FALSE, sortBy=c("or","adjp","fdr","none"))
{
    
    sortBy <- match.arg(sortBy)
    
    if(is.null(eTerm)){
        warnings("There is no enrichment in the 'eTerm' object.\n")
        return(NULL)
    }
    
    if(is(eTerm,'eTerm')){
		## when 'auto', will keep the significant terms
		df <- xEnrichViewer(eTerm, top_num="all")
		
		############
		if(!CI.one){
			ind <- which(df$CIl>1 | df$CIu<1)
			df <- df[ind,]
		}
		############
		
		if(top_num=='auto'){
			top_num <- sum(df$adjp<FDR.cutoff)
			if(top_num <= 1){
				warnings("There is no enrichment in the 'eTerm' object.\n")
				top_num <- 10
				message(sprintf("0 or 1 term found; instead the top %d terms sorted by '%s' are shown!", top_num, sortBy), appendLF=TRUE)
			}
		}
		df <- xEnrichViewer(eTerm, top_num=top_num, sortBy=sortBy)
		df$group <- 'group'
		df$ontology <- 'ontology'
		
	}else if(is(eTerm,'ls_eTerm') | is(eTerm,'data.frame') | is(eTerm,'tbl')){
	
		if(is(eTerm,'ls_eTerm')){
			## when 'auto', will keep the significant terms
			df <- eTerm$df
			
		}else if(is(eTerm,'data.frame') | is(eTerm,'tbl')){
			eTerm <- as.data.frame(eTerm)
			
			if(all(c('group','ontology','name','adjp','or','CIl','CIu') %in% colnames(eTerm))){
				df <- eTerm[,c('group','ontology','name','adjp','or','CIl','CIu')]
			
			}else if(all(c('group','name','adjp','or','CIl','CIu') %in% colnames(eTerm))){
				df <- eTerm[,c('group','name','adjp','or','CIl','CIu')]
				df$ontology <- 'ontology'
			
			}else if(all(c('ontology','name','adjp','or','CIl','CIu') %in% colnames(eTerm))){
				df <- eTerm[,c('ontology','name','adjp','or','CIl','CIu')]
				df$group <- 'group'
			
			}else if(all(c('name','adjp','or','CIl','CIu') %in% colnames(eTerm))){
				df <- eTerm[,c('name','adjp','or','CIl','CIu')]
				df$group <- 'group'
				df$ontology <- 'ontology'
			
			}else{
				warnings("The input data.frame does not contain required columns: c('group','ontology','name','adjp','or','CIl','CIu').\n")
				return(NULL)
			}
			
		}
		
		## columns are ordered as indicated by inputs
		df$group <- factor(df$group, levels=unique(df$group))
		
		############
		if(!CI.one){
			ind <- which(df$CIl>1 | df$CIu<1)
			df <- df[ind,]
		}
		############
		
		or <- group <- ontology <- rank <- adjp <- NULL
		df <- df %>% dplyr::arrange(-or)
		if(top_num=='auto'){
			df <- subset(df, df$adjp<FDR.cutoff)
		}else{
			top_num <- as.integer(top_num)
			df <- as.data.frame(df %>% dplyr::group_by(group,ontology) %>% dplyr::group_by(rank=rank(-or),add=TRUE) %>% dplyr::filter(rank<=top_num & adjp<FDR.cutoff))
		}
		
	}

	##########################
	##########################
	if(nrow(df)==0){
		return(NULL)
	}
	##########################
	##########################

	## text wrap
	if(!is.null(wrap.width)){
		width <- as.integer(wrap.width)
		res_list <- lapply(df$name, function(x){
			x <- gsub('_', ' ', x)
			y <- strwrap(x, width=width)
			if(length(y)>1){
				paste0(y[1], '...')
			}else{
				y
			}
		})
		df$name <- unlist(res_list)
	}
	
	name <- fdr <- or <- CIl <- CIu <- NULL
	group <- ontology <- NULL
	
	df$fdr <- -log10(df$adjp)
	if(is.null(zlim)){
		tmp <- df$fdr
		zlim <- c(floor(min(tmp)), ceiling(max(tmp[!is.infinite(tmp)])))
	}
	df$fdr[df$fdr<=zlim[1]] <- zlim[1]
	df$fdr[df$fdr>=zlim[2]] <- zlim[2]
	
	## order by 'or', 'adjp'
	if(is(eTerm,'eTerm') & sortBy!='or'){
		df <- df[rev(1:nrow(df)),]
	}else{
		df <- df[with(df,order(group, ontology, or, -fdr)),]
	}
	df$name <- factor(df$name, levels=unique(df$name))

	###########################################
	bp <- ggplot(df, aes(x=name, y=log2(or), ymin=log2(CIl), ymax=log2(CIu), color=fdr))
	bp <- bp + geom_pointrange(size=0.3) + ylab(expression(log[2]("odds ratio")))
	bp <- bp + geom_hline(yintercept=0, color='black', linetype='dashed') + coord_flip()
	bp <- bp + theme_bw() + theme(legend.position="right",legend.direction = "horizontal",axis.title.y=element_blank(), axis.text.y=element_text(size=6,color="black"), axis.title.x=element_text(size=6,color="black"), axis.text.x=element_text(size=5.5,color="black"))
	bp <- bp + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	bp <- bp + scale_colour_gradientn(colors=xColormap(colormap)(ncolors), limits=zlim, breaks = c(0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28), guide=guide_colorbar(title=expression(-log[10]("FDR")),title.position="top",title.theme = element_text(size = 6),label.theme = element_text(size = 5.5),barwidth=unit(1.75
,"cm"),barheight=unit(0.15, "cm"),draw.ulim=FALSE,draw.llim=FALSE))
	
	## change font family to 'Arial'
	bp <- bp + theme(text=element_text(family=font.family))
	
	## put arrows on x-axis
	#bp <- bp + theme(axis.line.x=element_line(arrow=arrow(angle=30,length=unit(0.25,"cm"), type="open")))
	
	## x-axis (actually y-axis) position
	#bp <- bp + scale_y_continuous(position="bottom")
	
	# facet_grid: partitions a plot into a matrix of panels
	## group (columns), ontology (rows)
	ngroup <- length(unique(df$group))
	nonto <- length(unique(df$ontology))
	if(ngroup!=1 | nonto!=1){
		scales <- "free_y"
		space <- "free_y"
		
		if(ngroup==1){
			bp <- bp + facet_grid(ontology~., scales=scales, space=space, drop=drop)
		}else if(nonto==1){
			bp <- bp + facet_grid(.~group, scales=scales, space=space, drop=drop)
		}else{
			bp <- bp + facet_grid(ontology~group, scales=scales, space=space, drop=drop)
		}
		
		## strip
		bp <- bp + theme(strip.background=element_rect(fill="transparent",color="transparent"), strip.text=element_text(size=8,face="bold.italic"))
	}
	
	invisible(bp)
}

