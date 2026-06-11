#' Function to visualise a list of prioritised genes using advanced track plot
#'
#' \code{xPierTrackAdv} is supposed to visualise prioritised genes using advanced track plot. Internally, it calls the function 'xPierTrack' per gene.
#'
#' @param pNode an object of class "pNode" (or "sTarget" or "dTarget")
#' @param priority.top the number of the top targets used for track plot. By default, it is NULL meaning all targets are used
#' @param targets.query which genes in query will be visualised. If NULL, the target gene with the top priority will be displayed
#' @param window the maximum distance defining nearby genes around the target gene in query. By default it is 1e6
#' @param nearby the maximum number defining nearby genes around the target gene in query. By default it is NULL. If not NULL, it will overwrite the parameter 'window'
#' @param query.highlight logical to indicate whether the gene in query will be highlighted
#' @param track.ideogram logical to indicate whether ideogram track is shown. By default, it is TRUE
#' @param track.genomeaxis logical to indicate whether genome axis track is shown. By default, it is TRUE
#' @param name.datatrack the name for the data track. By default, it is "Priority index"
#' @param name.annotrack the name for the annotation track. By default, it is "Target genes"
#' @param GR.Gene the genomic regions of genes. By default, it is 'UCSC_knownGene', that is, UCSC known genes (together with genomic locations) based on human genome assembly hg19. It can be 'UCSC_knownCanonical', that is, UCSC known canonical genes (together with genomic locations) based on human genome assembly hg19. Alternatively, the user can specify the customised input. To do so, first save your RData file (containing an GR object) into your local computer, and make sure the GR object content names refer to Gene Symbols. Then, tell "GR.Gene" with your RData file name (with or without extension), plus specify your file RData path in "RData.location"
#' @param SNPs a input vector containing SNPs. SNPs should be provided as dbSNP ID (ie starting with rs). Alternatively, they can be in the format of 'chrN:xxx', where N is either 1-22 or X, xxx is genomic positional number; for example, 'chr16:28525386'. By default, it is NLL meaning the SNP annotation track will be not displayed
#' @param GR.SNP the genomic regions of SNPs. By default, it is 'dbSNP_GWAS', that is, SNPs from dbSNP (version 146) restricted to GWAS SNPs and their LD SNPs (hg19). It can be 'dbSNP_Common', that is, Common SNPs from dbSNP (version 146) plus GWAS SNPs and their LD SNPs (hg19). Alternatively, the user can specify the customised input. To do so, first save your RData file (containing an GR object) into your local computer, and make sure the GR object content names refer to dbSNP IDs. Then, tell "GR.SNP" with your RData file name (with or without extension), plus specify your file RData path in "RData.location". Note: you can also load your customised GR object directly
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to false for no display
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @param guid a valid (5-character) Global Unique IDentifier for an OSF project. See \code{\link{xRDataLoader}} for details
#' @param ... additional graphic parameters. For example, the parameter "strip" allows the panel title is hided (FALSE), shown (TRUE) or without the background (lattice::strip.custom(bg="transparent")); the parameter "layout" allows specification of the layout (the first element for the columns and the second element for the rows). See \url{http://www.rdocumentation.org/packages/lattice/topics/xyplot} for the complete list.
#' @return an object of class "trellis"
#' @note none
#' @export
#' @seealso \code{\link{xRDataLoader}}
#' @include xPierTrackAdv.r
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
#' xPierTrackAdv(pNode, RData.location=RData.location)
#' #dev.off()
#' xPierTrackAdv(pNode, priority.top=1000, nearby=20, RData.location=RData.location)
#' }

xPierTrackAdv <- function(pNode, priority.top=NULL, targets.query=NULL, window=1e6, nearby=NULL, query.highlight=TRUE, track.ideogram=TRUE, track.genomeaxis=TRUE, name.datatrack="Priority index", name.annotrack="Genes", GR.Gene=c("UCSC_knownGene","UCSC_knownCanonical"), SNPs=NULL, GR.SNP=c("dbSNP_GWAS","dbSNP_Common","dbSNP_Single"), verbose=TRUE, RData.location="http://galahad.well.ox.ac.uk/bigdata", guid=NULL, ...)
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
    gr_Gene <- xRDataLoader(GR.Gene[1], verbose=verbose, RData.location=RData.location, guid=guid)
    if(is.null(gr_Gene)){
    	GR.Gene <- "UCSC_knownGene"
		if(verbose){
			message(sprintf("Instead, %s will be used", GR.Gene), appendLF=TRUE)
		}
    	gr_Gene <- xRDataLoader(GR.Gene, verbose=verbose, RData.location=RData.location, guid=guid)
    }
    
    ## ONLY restricted to genes with genomic locations
	ind <- match(rownames(df_priority), names(gr_Gene))
	p_gr <- gr_Gene[ind[!is.na(ind)],]
	p_matrix <- df_priority[!is.na(ind),]
	
	## append genomic locations to GR object
	gr <- p_gr
	GenomicRanges::mcols(gr) <- cbind(GenomicRanges::mcols(gr), p_matrix)
	
	########################################
	if(is.null(targets.query)){
		ind <- 1
	}else{
		ind <- match(targets.query, names(gr))
		ind <- ind[!is.na(ind)]
		if(length(ind) == 0){
			ind <- 1
		}
	}
	targets.query.found <- names(gr)[ind]
	targets.chr <- as.character(GenomicRanges::seqnames(gr)[ind])
	targets.gene_chr <- paste0(targets.query.found, ' (', targets.chr, ')')
	########################################


	panel.function <- function(x){
		xPierTrack(pNode=pNode, priority.top=priority.top, target.query=x, window=window, nearby=nearby, query.highlight=query.highlight, track.ideogram=track.ideogram, track.genomeaxis=track.genomeaxis, name.datatrack=name.datatrack, name.annotrack=name.annotrack, GR.Gene=GR.Gene, SNPs=SNPs, GR.SNP=GR.SNP, verbose=FALSE, RData.location=RData.location, guid=guid, add=TRUE)
	}
	
	#df_genes <- data.frame(gene=c("CREB1","IL6","NOS3","TNF"))
	#strip <- lattice::strip.custom(bg="transparent")
	#layout <- c(1,4) # 1 column and 4 rows
	
	gene_chr <- gene <- NULL
	df_genes <- data.frame(gene_chr=targets.gene_chr, gene=targets.query.found, chr=targets.chr)
	df_genes$gene_chr <- factor(df_genes$gene_chr, levels=df_genes$gene_chr)
	res <- lattice::xyplot(1~gene | gene_chr, data=df_genes, panel=panel.function, scales=list(draw=FALSE), xlab=NULL, ylab=NULL, ...)
    
    return(res)
}
