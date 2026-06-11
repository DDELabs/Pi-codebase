xVisEvidenceAdv <- function(xTarget, g=NA, nodes=NULL, node.info=c("smart","none"), neighbor.order=1, neighbor.seed=TRUE, neighbor.top=NULL, largest.comp=TRUE, node.label.size=2, node.label.color='black', node.label.alpha=0.9, node.label.padding=0.5, node.label.arrow=0, node.label.force=0.1, node.shape=19, node.color.title='Pi\nrating', colormap='white-yellow-red', ncolors=64, zlim=c(0,5), node.size.range=5, title='', edge.color="orange", edge.color.alpha=0.5, edge.curve=0, edge.arrow.gap=0.025, pie.radius=NULL, pie.color='black', pie.color.alpha=1, pie.thick=0.1,...)
{

    node.info <- match.arg(node.info)
	
	subg <- xVisEvidence(xTarget=xTarget, g=g, nodes=nodes, node.info=node.info, neighbor.order=neighbor.order, neighbor.seed=neighbor.seed, neighbor.top=neighbor.top, largest.comp=largest.comp, show=FALSE)
	
	if(!is(subg,'igraph')){
		return(NULL)
	}
	if(ecount(subg)<=1){
		return(NULL)
	}
	
	## layout
	if(any(is.null(V(subg)$xcoord), is.null(V(subg)$ycoord))){
		#glayout <- igraph::layout_as_tree(subg,root=dnet::dDAGroot(subg),circular=TRUE,flip.y=TRUE)
		if(0){
			glayout <- igraph::layout_with_kk(subg)
			V(subg)$xcoord <- glayout[,1]
			V(subg)$ycoord <- glayout[,2]
		}else{
			subg <- subg %>% xLayout("graphlayouts.layout_with_stress")
		}
	}
	
	#################
	gp <- xGGnetwork(g=subg, node.label='vertex.label', node.label.size=node.label.size, node.label.color=node.label.color, node.label.alpha=node.label.alpha, node.label.padding=node.label.padding, node.label.arrow=node.label.arrow, node.label.force=node.label.force, node.shape=node.shape, node.xcoord='xcoord', node.ycoord='ycoord', node.color='priority', node.color.title=node.color.title, colormap=colormap, ncolors=ncolors, zlim=zlim, node.size.range=node.size.range, title=title, edge.color=edge.color, edge.color.alpha=edge.color.alpha, edge.curve=edge.curve, edge.arrow.gap=edge.arrow.gap,...)
	    #################
	

	#df <- gp$data_nodes
	df <- gp$data
	print(df)
	
	if(1){
		## previously
		columns <- c('dGene','pGene','fGene','nGene','eGene','cGene')
		df_sub <- as.data.frame(matrix(0, nrow=nrow(df), ncol=length(columns)+2))
		colnames(df_sub) <- c('x','y',columns)
		ind <- match(colnames(df_sub), colnames(df))
		df_sub[,!is.na(ind)] <- df[,ind[!is.na(ind)]]
		ind <- which(apply(df_sub[,c(-1,-2)],1,sum)!=0)
		if(length(ind)>0){
			df_sub <- df_sub[ind, ]
			df_sub_1 <- df_sub[,c(1,2)]
			df_sub_2 <- df_sub[,c(-1,-2)]
			df_sub_2[df_sub_2>=1] <- 1
			df_sub <- cbind(df_sub_1, df_sub_2)

			gp <- xPieplot(df_sub, columns, colormap='ggplot2', pie.radius=pie.radius, pie.color=pie.color, pie.color.alpha=pie.color.alpha, pie.thick=pie.thick, legend.title='Seed gene', gp=gp)
		}

	}else{
		#################
		## now (20211115)
		#################
		x <- y <- ycoord <- vertex.label <- . <- NULL

		df_sub_1 <- df %>% dplyr::select(x,y)
		df_sub_2 <- df %>% dplyr::select(ycoord:vertex.label) %>% dplyr::select(-1,-length(.))
		df_sub_2[df_sub_2>=1] <- 1
		ind <- which(apply(df_sub_2,1,sum)!=0)
		df_sub <- cbind(df_sub_1, df_sub_2) %>% dplyr::slice(ind)

		columns <- colnames(df_sub_2)
		gp <- xPieplot(df_sub, columns, colormap='ggplot2', pie.radius=pie.radius, pie.color=pie.color, pie.color.alpha=pie.color.alpha, pie.thick=pie.thick, legend.title='Seed gene', gp=gp)
	}
	
    invisible(gp)
}
