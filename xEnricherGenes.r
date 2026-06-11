EnricherGenes <- function(data, background=NULL, check.symbol.identity=FALSE, ontology=NA, ontology.customised=NULL, size.range=c(10,2000), min.overlap=5, which.distance=NULL, test=c("fisher","hypergeo","binomial"), background.annotatable.only=NULL, p.tail=c("one-tail","two-tails"), p.adjust.method=c("BH", "BY", "bonferroni", "holm", "hochberg", "hommel"), ontology.algorithm=c("none","pc","elim","lea"), elim.pvalue=1e-2, lea.depth=2, path.mode=c("all_paths","shortest_paths","all_shortest_paths"), true.path.rule=FALSE, verbose=TRUE, silent=FALSE, RData.location="http://galahad.well.ox.ac.uk/bigdata", guid=NULL)
{
    startT <- Sys.time()
    if(!silent){
    	message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=TRUE)
    	message("", appendLF=TRUE)
    }else{
    	verbose <- FALSE
    }
    ####################################################################################
    
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    #ontology <- match.arg(ontology)
    #ontology <- ontology[1]
    test <- match.arg(test)
    p.tail <- match.arg(p.tail)
    p.adjust.method <- match.arg(p.adjust.method)
    ontology.algorithm <- match.arg(ontology.algorithm)
    path.mode <- match.arg(path.mode)
    p.tail <- match.arg(p.tail)
    
    ############
    if(length(data)==0){
    	return(NULL)
    }
    ############
    
    if (is.vector(data)){
        data <- unique(data)
    }else{
        warnings("The input data must be a vector.\n")
        return(NULL)
    }
    
    data <- as.character(data)
    
    #################################
 	# aOnto <- xDefineOntology(ontology, ontology.customised=ontology.customised, verbose=verbose, RData.location=RData.location, guid=guid)
  aOnto <- ontology
  class(aOnto)
 	g <- aOnto$g
 	class(g)
 	anno <- aOnto$anno
 	class(anno)
 	if(is.null(g)){
		warnings("There is no input for the ontology.\n")
        return(NULL)
	}
    #################################

    if(is.null(ontology.customised)){
		## convert gene symbol to entrz gene for both input data of interest and the input background (if given)
		if(verbose){
			now <- Sys.time()
			message(sprintf("Do gene mapping from Symbols to EntrezIDs (%s) ...", as.character(now)), appendLF=TRUE)
		}
		data <- xSymbol2GeneID(data, check.symbol.identity=check.symbol.identity, verbose=verbose, RData.location=RData.location, guid=guid)
		data <- data[!is.na(data)]
		if(length(background)>0){
			background <- xSymbol2GeneID(background, check.symbol.identity=check.symbol.identity, verbose=verbose, RData.location=RData.location, guid=guid)
			background <- background[!is.na(background)]
		}
    }
    
    #############################################################################################
    
    if(verbose){
        now <- Sys.time()
        message(sprintf("\n#######################################################"))
        message(sprintf("'xEnricher' is being called (%s):", as.character(now)), appendLF=TRUE)
        message(sprintf("#######################################################"))
    }
    eTerm <- xEnricher(data=data, annotation=anno, g=g, background=background, size.range=size.range, min.overlap=min.overlap, which.distance=which.distance, test=test, background.annotatable.only=background.annotatable.only, p.tail=p.tail, p.adjust.method=p.adjust.method, ontology.algorithm=ontology.algorithm, elim.pvalue=elim.pvalue, lea.depth=lea.depth, path.mode=path.mode, true.path.rule=true.path.rule, verbose=verbose)
	
	# replace EntrezGenes with gene symbols	
	if(is.null(ontology.customised) & is(eTerm,"eTerm")){
		## load Enterz Gene information
		EG <- xRDataLoader(RData.customised=paste('org.Hs.eg', sep=''), RData.location=RData.location, guid=guid, verbose=verbose)
		allGeneID <- EG$gene_info$GeneID
		allSymbol <- as.vector(EG$gene_info$Symbol)
	
		## overlap
		overlap <- eTerm$overlap
		overlap_symbols <- lapply(overlap,function(x){
			ind <- match(x, allGeneID)
			allSymbol[ind]
		})
		eTerm$overlap <- overlap_symbols
		## data
		eTerm$data <- allSymbol[match(eTerm$data,allGeneID)]
		## background
		eTerm$background <- allSymbol[match(eTerm$background,allGeneID)]		
		
		## annotation
		annotation <- eTerm$annotation
		annotation_symbols <- lapply(annotation,function(x){
			ind <- match(x, allGeneID)
			allSymbol[ind]
		})
		eTerm$annotation <- annotation_symbols
	}
	
	if(verbose){
        now <- Sys.time()
        message(sprintf("#######################################################"))
        message(sprintf("'xEnricher' has been finished (%s)!", as.character(now)), appendLF=TRUE)
        message(sprintf("#######################################################\n"))
    }
    
    ####################################################################################
    endT <- Sys.time()
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    
    if(!silent){
    	message(paste(c("\nEnd at ",as.character(endT)), collapse=""), appendLF=TRUE)
    	message(paste(c("Runtime in total (xEnricherGenes): ",runTime," secs\n"), collapse=""), appendLF=TRUE)
    }
    
    invisible(eTerm)
}
