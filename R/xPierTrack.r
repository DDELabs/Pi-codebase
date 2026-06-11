#' Function to visualise a prioritised gene using track plot
#'
#' \code{xPierTrack} is supposed to visualise a prioritised gene using track plot. Priority for the gene in query is displayed on the data track and nearby genes on the annotation track. Genomic locations on the X-axis are indicated on the X-axis, and the gene in query is highlighted. If SNPs are also provided, SNP annotation track will be also displayed at the bottom.
#'
#' @param pNode an object of class "pNode" (or "sTarget" or "dTarget")
#' @param priority.top the number of the top targets used for track plot. By default, it is NULL meaning all targets are used
#' @param target.query which gene in query will be visualised. If NULL, the target gene with the top priority will be displayed
#' @param window the maximum distance defining nearby genes around the target gene in query. By default it is 1e6
#' @param nearby the maximum number defining nearby genes around the target gene in query. By default it is NULL. If not NULL, it will overwrite the parameter 'window'
#' @param query.highlight logical to indicate whether the gene in query will be highlighted
#' @param track.ideogram logical to indicate whether ideogram track is shown. By default, it is TRUE
#' @param track.genomeaxis logical to indicate whether genome axis track is shown. By default, it is TRUE
#' @param name.datatrack the name for the data track. By default, it is "Priority index"
#' @param name.annotrack the name for the annotation track. By default, it is "Genes". If NULL, the title for annotation track will be hided
#' @param GR.Gene the genomic regions of genes. By default, it is 'UCSC_knownGene', that is, UCSC known genes (together with genomic locations) based on human genome assembly hg19. It can be 'UCSC_knownCanonical', that is, UCSC known canonical genes (together with genomic locations) based on human genome assembly hg19. Alternatively, the user can specify the customised input. To do so, first save your RData file (containing an GR object) into your local computer, and make sure the GR object content names refer to Gene Symbols. Then, tell "GR.Gene" with your RData file name (with or without extension), plus specify your file RData path in "RData.location"
#' @param SNPs a input vector containing SNPs. SNPs should be provided as dbSNP ID (ie starting with rs). Alternatively, they can be in the format of 'chrN:xxx', where N is either 1-22 or X, xxx is genomic positional number; for example, 'chr16:28525386'. By default, it is NULL meaning the SNP annotation track will be not displayed
#' @param max.num.SNPs the maximum number (50 by default) of SNPs to be shown. If NULL, no such restriction. Also this parameter only works when the SNP annotation track is enabled 
#' @param GR.SNP the genomic regions of SNPs. By default, it is 'dbSNP_GWAS', that is, SNPs from dbSNP (version 146) restricted to GWAS SNPs and their LD SNPs (hg19). It can be 'dbSNP_Common', that is, Common SNPs from dbSNP (version 146) plus GWAS SNPs and their LD SNPs (hg19). Alternatively, the user can specify the customised input. To do so, first save your RData file (containing an GR object) into your local computer, and make sure the GR object content names refer to dbSNP IDs. Then, tell "GR.SNP" with your RData file name (with or without extension), plus specify your file RData path in "RData.location". Note: you can also load your customised GR object directly
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to false for no display
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @param guid a valid (5-character) Global Unique IDentifier for an OSF project. See \code{\link{xRDataLoader}} for details
#' @param ... additional graphic parameters. For example, the parameter "add" allows the plot added to an existing plotting canvas without re-initialising. See \url{http://www.rdocumentation.org/packages/Gviz/topics/plotTracks} for the complete list.
#' @return a list of GenomeGraph tracks, each one augmented by the computed image map coordinates in the 'imageMap' slot, along with the additional 'ImageMap' object 'titles' containing information about the title panels.
#' @note none
#' @export
#' @seealso \code{\link{xRDataLoader}}, \code{\link{xSNPlocations}}, \code{\link{xGR}}
#' @include xPierTrack.r
#' @examples
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' \dontrun{
#' # a) provide the SNPs with the significance info
#' ## get lead SNPs reported in AS GWAS and their significance info (p-values)
#' #data.file <- "http://galahad.well.ox.ac.uk/bigdata/AS.txt"
#' #AS <- read.delim(data.file, header=TRUE, stringsAsFactors=FALSE)
#' ImmunoBase <- xRDataLoader(RData.customised='ImmunoBase', RData.location=RData.location)
#' gr <- ImmunoBase$AS$variants
#' AS <- as.data.frame(GenomicRanges::mcols(gr)[, c('Variant','Pvalue')])
#'
#' # b) perform priority analysis
#' pNode <- xPierSNPs(data=AS, include.eQTL="JKng_mono", include.HiC='Monocytes', network="PCommonsUN_medium", restart=0.7, RData.location=RData.location)
#'
#' # c) track plot
#' library(Gviz)
#' #pdf(file="Gene_tracks.pdf", height=4, width=10, compress=TRUE)
#' xPierTrack(pNode, RData.location=RData.location)
#' #dev.off()
#' xPierTrack(pNode, priority.top=1000, nearby=20, RData.location=RData.location)
#' }

xPierTrack <- function(pNode, priority.top=NULL, target.query=NULL, window=1e6, nearby=NULL, query.highlight=TRUE, track.ideogram=TRUE, track.genomeaxis=TRUE, name.datatrack="5-star rating\n(Priority index)", name.annotrack="Targets", GR.Gene=c("UCSC_knownGene","UCSC_knownCanonical"), SNPs=NULL, max.num.SNPs=50, GR.SNP=c("dbSNP_GWAS","dbSNP_Common","dbSNP_Single"), verbose=TRUE, RData.location="http://galahad.well.ox.ac.uk/bigdata", guid=NULL, ...)
{

    if(is(pNode,"pNode")){
        df_priority <- pNode$priority[, c("seed","weight","priority")]
    }else if(is(pNode,"sTarget") | is(pNode,"dTarget")){
    	df_priority <- pNode$priority[, c("name","rank","rating")]
    	df_priority$priority <- df_priority$rating
    }else{
    	stop("The function must apply to a 'pNode' or 'sTarget' or 'dTarget' object.\n")
    }

   ###############
	## priority top
	if(!is.null(priority.top)){
		priority.top <- as.integer(priority.top)
		if(priority.top > nrow(df_priority)){
			priority.top <- nrow(df_priority)
		}else if(priority.top <= 1){
			priority.top <- nrow(df_priority)
		}
	}else{
		priority.top <- nrow(df_priority)
	}
    df_priority <- df_priority[1:priority.top,]
    ###############
    
	if(verbose){
		now <- Sys.time()
		message(sprintf("Load positional information for Genes (%s) ...", as.character(now)), appendLF=TRUE)
	}
    gr_Gene <- xRDataLoader(RData.customised=GR.Gene[1], verbose=verbose, RData.location=RData.location, guid=guid)
    if(is.null(gr_Gene)){
    	GR.Gene <- "UCSC_knownGene"
		if(verbose){
			message(sprintf("Instead, %s will be used", GR.Gene), appendLF=TRUE)
		}
    	gr_Gene <- xRDataLoader(RData.customised=GR.Gene, verbose=verbose, RData.location=RData.location, guid=guid)
    }
    
    ## ONLY restricted to genes with genomic locations
	ind <- match(rownames(df_priority), names(gr_Gene))
	p_gr <- gr_Gene[ind[!is.na(ind)],]
	p_matrix <- df_priority[!is.na(ind),]
	
	## append genomic locations to GR object
	gr <- p_gr
	GenomicRanges::mcols(gr) <- cbind(GenomicRanges::mcols(gr), p_matrix)
	
	########################################
	if(!is.null(target.query)){
		ind <- match(target.query, names(gr))
		if(length(ind) == 0 | is.na(ind)){
			#warning(sprintf("\tNo found! %s will be used instead of %s", names(gr[1]), target.query), appendLF=TRUE)
			#target.query <- names(gr[1])
			warning(sprintf("\tNo found for queried %s", target.query), appendLF=TRUE)
			return(NULL)
		}
	}else{
		target.query <- names(gr[1])
	}
	########################################

	## for ideogram track
	if(track.ideogram){
		chr <- as.character(unique(GenomicRanges::seqnames(gr[target.query])))
		if(length(suppressWarnings(tryCatch(itrack <- Gviz::IdeogramTrack(genome="hg19", chromosome=chr, showBandId=TRUE, cex.bands=0.8), error=function(e) e, warning=function(w) w)))==2){
			itrack <- NULL
		}
	}else{
		itrack <- NULL
	}
	## for genome axis track
	if(track.genomeaxis){
		gtrack <- Gviz::GenomeAxisTrack()
	}else{
		gtrack <- NULL
	}
	
	####################
	if(is.null(nearby)){
		## a window (eg 1e6) of upstream and downstream from the gene in query
		q2r <- as.matrix(as.data.frame(suppressWarnings(GenomicRanges::findOverlaps(query=gr, subject=gr[target.query], maxgap=window-1, minoverlap=0L, type="any", select="all", ignore.strand=TRUE))))
		gr_sub <- gr[q2r[,1],]
		
	}else{
		dists <- as.matrix(as.data.frame(suppressWarnings(GenomicRanges::distance(x=gr, y=gr[target.query], select="all", ignore.strand=TRUE))))
		x <- sort(dists)
		gr_sub <- gr[which(dists <= x[min(nearby+1,length(x))])]
	}
	
	## for annotation and data track
	gr_sub$group <- gr_sub$Symbol
	gr_sub$query <- factor(ifelse(gr_sub$Symbol==target.query, "q", "notq"))
	x <- NULL
	if(is.null(name.annotrack)){
		atrack <- Gviz::AnnotationTrack(gr_sub, name=name.annotrack, feature=gr_sub$query, shape="box", col="transparent", fill="grey", background.title="transparent")
	}else{
		atrack <- Gviz::AnnotationTrack(gr_sub, name=name.annotrack, feature=gr_sub$query, shape="box", col="transparent", fill="grey", background.title="darkblue")
	}
	dtrack <- Gviz::DataTrack(gr_sub, data="priority", baseline=0, ylim=c(min(gr$priority),max(gr$priority)), ytick=0:max(gr$priority), type="h", name=name.datatrack, strand="+", background.title="transparent", background.panel="#FFFEDB", cex.axis=0.7, col.title="darkblue",col.axis="grey",col.baseline="transparent",col.line="darkblue",lwd=2,transformation=function(x) x)
	#Gviz::availableDisplayPars(dtrack)
	
	## for SNP annotation track
	gr_SNP <- xSNPlocations(data=unique(SNPs), GR.SNP=GR.SNP, verbose=verbose, RData.location=RData.location, guid=guid)
	if(!is.null(gr_SNP)){
		## for genes
		df_gr_sub <- GenomicRanges::as.data.frame(gr_sub, row.names=NULL)
		df <- data.frame(chr=as.character(unique(df_gr_sub[,1])), start=min(df_gr_sub[,2:3]), end=max(df_gr_sub[,2:3]))
		gr_tmp <- xGR(df, format="data.frame", RData.location=RData.location, guid=guid)
		## find snps
		q2r <- as.matrix(as.data.frame(suppressWarnings(GenomicRanges::findOverlaps(query=gr_SNP, subject=gr_tmp, maxgap=-1L, minoverlap=0L, type="any", select="all", ignore.strand=TRUE))))
		gr_SNP_sub <- gr_SNP[q2r[,1],]
		## track
		if(length(gr_SNP_sub)==0){
			atrack_SNP <- NULL
		}else{
			if(!is.null(max.num.SNPs)){
				if(length(gr_SNP_sub)>max.num.SNPs){
					gr_SNP_sub <- gr_SNP_sub[1:max.num.SNPs]
				}
			}
		
			gr_SNP_sub$group <- names(gr_SNP_sub)
			atrack_SNP <- Gviz::AnnotationTrack(gr_SNP_sub, name="", shape="box", col="transparent", background.title="transparent")
		}
	}else{
		atrack_SNP <- NULL
	}

	#########
	## plot tracks
	ls_tracks <- list(itrack,gtrack,dtrack,atrack,atrack_SNP)
	ls_tracks <- base::Filter(base::Negate(is.null), ls_tracks)
	if(query.highlight){
		res <- Gviz::plotTracks(ls_tracks, groupAnnotation="group", q="salmon", notq="grey", cex.title=0.8, col=NULL, extend.left=0.05, extend.right=0.05, ...)
	}else{
		res <- Gviz::plotTracks(ls_tracks, groupAnnotation="group", cex.title=1, col=NULL, extend.left=0.05, extend.right=0.05, ...)
	}
	
    
    invisible(res)
}


