#' Function to prioritise pathways based on GSEA analysis of prioritised genes
#'
#' \code{xPierGSEA} is supposed to prioritise pathways given prioritised genes and the ontology in query. It is done via gene set enrichment analysis (GSEA). It returns an object of class "eGSEA". 
#
#' @param pNode an object of class "pNode" (or "sTarget" or "dTarget"). Alternatively, it can be a data frame with two columns ('priority' and 'rank')
#' @param priority.top the number of the top targets used for GSEA. By default, it is NULL meaning all targets are used
#' @param ontology the ontology supported currently. It can be "GOBP" for Gene Ontology Biological Process, "GOMF" for Gene Ontology Molecular Function, "GOCC" for Gene Ontology Cellular Component, "PS" for phylostratific age information, "PS2" for the collapsed PS version (inferred ancestors being collapsed into one with the known taxonomy information), "SF" for SCOP domain superfamilies, "Pfam" for Pfam domain families, "DO" for Disease Ontology, "HPPA" for Human Phenotype Phenotypic Abnormality, "HPMI" for Human Phenotype Mode of Inheritance, "HPCM" for Human Phenotype Clinical Modifier, "HPMA" for Human Phenotype Mortality Aging, "MP" for Mammalian Phenotype, "EF" for Experimental Factor Ontology (used to annotate GWAS Catalog genes), Drug-Gene Interaction database ("DGIdb") for drugable categories, tissue-specific eQTL-containing genes from GTEx ("GTExV4", "GTExV6p" and "GTExV7"), crowd extracted expression of differential signatures from CREEDS ("CreedsDisease", "CreedsDiseaseUP", "CreedsDiseaseDN", "CreedsDrug", "CreedsDrugUP", "CreedsDrugDN", "CreedsGene", "CreedsGeneUP" and "CreedsGeneDN"), KEGG pathways (including 'KEGG' for all, 'KEGGmetabolism' for 'Metabolism' pathways, 'KEGGgenetic' for 'Genetic Information Processing' pathways, 'KEGGenvironmental' for 'Environmental Information Processing' pathways, 'KEGGcellular' for 'Cellular Processes' pathways, 'KEGGorganismal' for 'Organismal Systems' pathways, and 'KEGGdisease' for 'Human Diseases' pathways), and the molecular signatures database (Msigdb, including "MsigdbH", "MsigdbC1", "MsigdbC2CGP", "MsigdbC2CPall", "MsigdbC2CP", "MsigdbC2KEGG", "MsigdbC2REACTOME", "MsigdbC2BIOCARTA", "MsigdbC3TFT", "MsigdbC3MIR", "MsigdbC4CGN", "MsigdbC4CM", "MsigdbC5BP", "MsigdbC5MF", "MsigdbC5CC", "MsigdbC6", "MsigdbC7")
#' @param customised.genesets a list each containing gene symbols. By default, it is NULL. If the list provided, it will overtake the previous parameter "ontology"
#' @param size.range the minimum and maximum size of members of each term in consideration. By default, it sets to a minimum of 10 but no more than 500
#' @param p.adjust.method the method used to adjust p-values. It can be one of "BH", "BY", "bonferroni", "holm", "hochberg" and "hommel". The first two methods "BH" (widely used) and "BY" control the false discovery rate (FDR: the expected proportion of false discoveries amongst the rejected hypotheses); the last four methods "bonferroni", "holm", "hochberg" and "hommel" are designed to give strong control of the family-wise error rate (FWER). Notes: FDR is a less stringent condition than FWER
#' @param path.mode the mode of paths induced by vertices/nodes with input annotation data. It can be "all_paths" for all possible paths to the root, "shortest_paths" for only one path to the root (for each node in query), "all_shortest_paths" for all shortest paths to the root (i.e. for each node, find all shortest paths with the equal lengths)
#' @param weight an integer specifying score weight. It can be "0" for unweighted (an equivalent to Kolmogorov-Smirnov, only considering the rank), "1" for weighted by input gene score (by default), and "2" for over-weighted, and so on
#' @param seed an integer specifying the seed
#' @param nperm the number of random permutations. For each permutation, gene-score associations will be permutated so that permutation of gene-term associations is realised
#' @param fast logical to indicate whether to fast calculate GSEA resulting. By default, it sets to true, but not necessarily does so. It will depend on whether the package "fgsea" has been installed
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true
#' @param silent logical to indicate whether the messages will be silent completely. By default, it sets to false. If true, verbose will be forced to be false
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @param guid a valid (5-character) Global Unique IDentifier for an OSF project. See \code{\link{xRDataLoader}} for details
#' @return
#' an object of class "eGSEA", a list with following components:
#' \itemize{
#'  \item{\code{df_summary}: a data frame of nTerm x 9 containing gene set enrichment analysis result, where nTerm is the number of terms/genesets, and the 9 columns are "setID" (i.e. "Term ID"), "name" (i.e. "Term Name"), "nAnno" (i.e. number in members annotated by a term), "nLead" (i.e. number in members as leading genes), "peak" (i.e. the rank at peak), "total" (i.e. the total number of genes analysed), "es" (i.e. enrichment score), "nes" (i.e. normalised enrichment score; enrichment score but after being normalised by gene set size), "pvalue" (i.e. nominal p value), "adjp" (i.e. adjusted p value; p value but after being adjusted for multiple comparisons), "distance" (i.e. term distance or metadata)}
#'  \item{\code{leading}: a list of gene sets, each storing leading gene info (i.e. the named vector with names for gene symbols and elements for priority rank). Always, gene sets are identified by "setID"}
#'  \item{\code{full}: a list of gene sets, each storing full info on gene set enrichment analysis result (i.e. a data frame of nGene x 6, where nGene is the number of genes, and the 6 columns are "GeneID", "Rank" for priority rank, "Score" for priority score, "RES" for running enrichment score,  "Hits" for gene set hits info with 1 for gene hit, 2 for leading gene hit, 3 for the point defining leading genes, 0 for no hit), and "Symbol" for gene symbols. Always, gene sets are identified by "setID"}
#'  \item{\code{cross}: a matrix of nTerm X nTerm, with an on-diagnal cell for the leading genes observed in an individaul term, and off-diagnal cell for the overlapped leading genes shared between two terms}
#' }
#' @note none
#' @export
#' @seealso \code{\link{xSymbol2GeneID}}, \code{\link{xRDataLoader}}, \code{\link{xDAGanno}}
#' @include xPierGSEA.r
#' @examples
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' \dontrun{
#' # a) provide the seed nodes/genes with the weight info
#' ## load ImmunoBase
#' ImmunoBase <- xRDataLoader(RData.customised='ImmunoBase', RData.location=RData.location)
#' ## get genes within 500kb away from AS GWAS lead SNPs
#' seeds.genes <- ImmunoBase$AS$genes_variants
#' ## seeds weighted according to distance away from lead SNPs
#' data <- 1- seeds.genes/500000
#'
#' # b) perform priority analysis
#' pNode <- xPierGenes(data=data, network="PCommonsDN_medium",restart=0.7, RData.location=RData.location)
#' 
#' # c) do pathway-level priority using GSEA
#' eGSEA <- xPierGSEA(pNode=pNode, ontology="DGIdb", nperm=2000, RData.location=RData.location)
#' bp <- xGSEAbarplot(eGSEA, top_num="auto", displayBy="nes")
#' gp <- xGSEAdotplot(eGSEA, top=1)
#' }

xPierGSEA <- function(pNode, priority.top=NULL, ontology=c("GOBP","GOMF","GOCC","PS","PS2","SF","Pfam","DO","HPPA","HPMI","HPCM","HPMA","MP", "EF", "MsigdbH", "MsigdbC1", "MsigdbC2CGP", "MsigdbC2CPall", "MsigdbC2CP", "MsigdbC2KEGG", "MsigdbC2REACTOME", "MsigdbC2BIOCARTA", "MsigdbC3TFT", "MsigdbC3MIR", "MsigdbC4CGN", "MsigdbC4CM", "MsigdbC5BP", "MsigdbC5MF", "MsigdbC5CC", "MsigdbC6", "MsigdbC7", "DGIdb", "GTExV4", "GTExV6p", "GTExV7", "CreedsDisease", "CreedsDiseaseUP", "CreedsDiseaseDN", "CreedsDrug", "CreedsDrugUP", "CreedsDrugDN", "CreedsGene", "CreedsGeneUP", "CreedsGeneDN", "KEGG","KEGGmetabolism","KEGGgenetic","KEGGenvironmental","KEGGcellular","KEGGorganismal","KEGGdisease"), customised.genesets=NULL, size.range=c(10,500), p.adjust.method=c("BH", "BY", "bonferroni", "holm", "hochberg", "hommel"), path.mode=c("all_paths","shortest_paths","all_shortest_paths"), weight=1, seed=825, nperm=2000, fast=TRUE, verbose=TRUE, silent=FALSE, RData.location="http://galahad.well.ox.ac.uk/bigdata", guid=NULL)
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
    ontology <- match.arg(ontology)
    path.mode <- match.arg(path.mode)
    
    weight <- as.integer(weight)
    
    if(is(pNode,"pNode")){
        df_priority <- pNode$priority[, c("seed","weight","priority","rank")]
    }else if(is(pNode,"sTarget") | is(pNode,"dTarget")){
    	df_priority <- pNode$priority[, c("name","rank","rating")]
    	df_priority$priority <- df_priority$rating
    }else if(is(pNode,"data.frame")){
    	df_priority <- pNode[,c(1:2)]
    	colnames(df_priority) <- c("priority","rank")
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
    
    ## convert gene symbols to entrez geneid
	name_GeneID <- xSymbol2GeneID(rownames(df_priority), check.symbol.identity=FALSE, verbose=verbose, RData.location=RData.location, guid=guid)
	
	ind <- which(!is.na(name_GeneID))
	GeneID <- name_GeneID[ind]
	priority <- df_priority$priority[ind]
	# also remove duplicated GeneID
	flag <- !duplicated(GeneID)
    data <- data.frame(priority=priority[flag], row.names=GeneID[flag], stringsAsFactors=FALSE)
    
    #############################################################################################

    if(is.null(customised.genesets)){
	
		if(verbose){
			now <- Sys.time()
			message(sprintf("Load the ontology %s and its gene annotations (%s) ...", ontology, as.character(now)), appendLF=TRUE)
		}

		#########
		## load GS information
		## flag the simplified version of PS
		flag_PS2 <- FALSE
		if(ontology=="PS2"){
			flag_PS2 <- TRUE
			ontology <- "PS"
		}
		GS <- xRDataLoader(RData.customised=paste('org.Hs.eg', ontology, sep=''), RData.location=RData.location, guid=guid, verbose=verbose)
		
		################
		if(flag_PS2){
			tmp <- as.character(unique(GS$set_info$name))
			inds <- sapply(tmp,function(x) which(GS$set_info$name==x))
		
			## new set_info
			set_info <- data.frame()
			for(i in 1:length(inds)){
				set_info<- rbind(set_info,as.matrix(GS$set_info[max(inds[[i]]),]))
			}
			## new gs
			gs <- list()
			for(i in 1:length(inds)){
				gs[[i]] <- unlist(GS$gs[inds[[i]]], use.names=FALSE)
			}
			names(gs) <- rownames(set_info)
		
			## new GS
			GS$set_info <- set_info
			GS$gs <- gs
		}
		################
		
		#########
		## get annotation information
		anno <- GS$gs

		#########
		## get ontology information
		## check the eligibility for the ontology
		all.ontologies <- c("GOBP","GOMF","GOCC","DO","HPPA","HPMI","HPCM","HPMA","MP","EF")
		flag_ontology <- ontology %in% all.ontologies
    	
    	if(flag_ontology){
			g <- xRDataLoader(paste0('ig.',ontology), RData.location=RData.location, guid=guid, verbose=verbose)
			true.path.rule <- TRUE
			
		}else{
			nodes <- data.frame(name=as.character(GS$set_info$setID), term_id=as.character(GS$set_info$setID), term_name=as.character(GS$set_info$name), term_distance=as.character(GS$set_info$distance), stringsAsFactors=FALSE)
			nodes <- rbind(nodes, c('root','root','root','root'))
			relations <- data.frame(from='root', to=nodes$name)
			g <- igraph::graph.data.frame(d=relations, directed=TRUE, vertices=nodes)
			true.path.rule <- FALSE
		}
	
		# obtain the induced subgraph according to the input annotation data
		subg <- xDAGanno(g=g, annotation=anno, path.mode=path.mode, true.path.rule=true.path.rule, verbose=verbose)
		anno <- V(subg)$anno
		names(anno) <- V(subg)$name
		g <- subg
	
	}else{
        if(is.list(customised.genesets)){
            if(is.null(names(customised.genesets))){
                names(customised.genesets) <- paste("C", 1:length(customised.genesets), sep="")
            }
            #######
            customised.genesets <- customised.genesets[lapply(customised.genesets,length)>0]
            if(length(customised.genesets)==0){
            	return(NULL)
            }
            #######
			anno <- lapply(customised.genesets, function(x){
				GeneID <- xSymbol2GeneID(x, check.symbol.identity=FALSE, verbose=verbose, RData.location=RData.location, guid=guid)
				GeneID <- GeneID[!is.na(GeneID)]
				return(GeneID)
			})
			nodes <- data.frame(name=names(customised.genesets), term_id=names(customised.genesets), term_name=names(customised.genesets), term_distance=rep(1, length(customised.genesets)), stringsAsFactors=FALSE)
			nodes <- rbind(nodes, c('root','root','root','root'))
			relations <- data.frame(from='root', to=nodes$name)
			g <- igraph::graph.data.frame(d=relations, directed=TRUE, vertices=nodes)
            true.path.rule <- FALSE
			subg <- xDAGanno(g=g, annotation=anno, path.mode=path.mode, true.path.rule=true.path.rule, verbose=verbose)
			anno <- V(subg)$anno
			names(anno) <- V(subg)$name
			g <- subg
			
        }else{
			stop("There is no input for the ontology.\n")
		}
	}
	
	#####
	ind <- which(sapply(anno,length) >= size.range[1] & sapply(anno,length) <= size.range[2])
	if(length(ind)>0){
		anno <- anno[ind]
	}else{
		return(NULL)
	}
	#####
	
    #############################################################################################
    
    flag_fgsea <- FALSE
    pkgs <- c("fgsea")
    if(all(pkgs %in% rownames(utils::installed.packages()))){
        tmp <- sapply(pkgs, function(pkg) {
            requireNamespace(pkg, quietly=TRUE)
        })
        if(all(tmp)){
        	flag_fgsea <- TRUE
        }
    }
    
    if(flag_fgsea & fast){
    
		if(verbose){
			now <- Sys.time()
			message(sprintf("\n#######################################################", appendLF=TRUE))
			message(sprintf("'fgsea' from the fgsea package is being called (%s):", as.character(now)), appendLF=TRUE)
			message(sprintf("#######################################################", appendLF=TRUE))
		}
		
		stats <- data$priority
		names(stats) <- rownames(data)
		fgseaRes <- fgsea::fgsea(pathway=anno, stats=stats, minSize=size.range[1], maxSize=size.range[2], gseaParam=weight, nperm=nperm)
		tab <- data.frame(setID         = fgseaRes$pathway,
						   ES           = fgseaRes$ES,
						   nES          = fgseaRes$NES,
						   pvalue       = fgseaRes$pval,
						   adjp         = fgseaRes$padj,
						   setSize      = fgseaRes$size,
						   stringsAsFactors=FALSE
						  )
		pvalue <- nES <- ES <- NULL
		res <- tab[with(tab,order(pvalue,-nES,-ES)),]
	
		if(verbose){
			now <- Sys.time()
			message(sprintf("#######################################################", appendLF=TRUE))
			message(sprintf("'fgsea' has been finished (%s)!", as.character(now)), appendLF=TRUE)
			message(sprintf("#######################################################\n", appendLF=TRUE))
		}
	
	}else{
	
		if(verbose){
			now <- Sys.time()
			message(sprintf("\n#######################################################", appendLF=TRUE))
			message(sprintf("'dGSEA' from the dnet package is being called (%s):", as.character(now)), appendLF=TRUE)
			message(sprintf("#######################################################", appendLF=TRUE))
		}
		
		if(!is.null(seed)) set.seed(seed)
		suppressMessages(eTerm <- dnet::dGSEA(data=data, identity="entrez", check.symbol.identity=FALSE, ontology="Customised", customised.genesets=anno, sizeRange=size.range, which_distance=NULL, weight=weight, nperm=nperm, fast=TRUE, sigTail="one-tail", p.adjust.method=p.adjust.method, verbose=verbose, RData.location=RData.location))
		
		if(!is.null(eTerm)){
			res <- dnet::dGSEAview(eTerm, which_sample=1, top_num=NULL, sortBy="pvalue", decreasing=TRUE, details=TRUE)
			if(nrow(res)>0){
				res <- res[,c("setID","ES","nES","pvalue","adjp","setSize")]
				rownames(res) <- NULL
			}else{
				return(NULL)
			}
		}else{
			return(NULL)
		}
		
		if(verbose){
			now <- Sys.time()
			message(sprintf("#######################################################", appendLF=TRUE))
			message(sprintf("'dGSEA' has been finished (%s)!", as.character(now)), appendLF=TRUE)
			message(sprintf("#######################################################\n", appendLF=TRUE))
		}
	}

    #############################################################################################
	### get Leading genes
    if(TRUE){

		# replace EntrezGenes with gene symbols	
		## load Enterz Gene information
		EG <- xRDataLoader('org.Hs.eg', RData.location=RData.location, guid=guid, verbose=verbose)
		allGeneID <- EG$gene_info$GeneID
		allSymbol <- as.vector(EG$gene_info$Symbol)

		########################
		geneid <- rownames(data)
		nGene <- nrow(data)
		## score rank
		rank.score <- data$priority
		ind <- order(rank.score, decreasing=TRUE)
		rank.score.sorted <- rank.score[ind]
		geneid.sorted <- geneid[ind]
    
		ls_df_leading <- lapply(res$setID, function(x){
			#df_leading <- visGSEA(eTerm, which_sample=1, which_term=x, plot=FALSE)
			
			## initialisation
			nHit <- length(anno[[x]])
			nMiss <- nGene - nHit
			## observed
			observed.point <- rep(-1/nMiss, nGene)
			flag <- match(anno[[x]], geneid.sorted)
			###### remove NA
			flag <- flag[!is.na(flag)]
			######
			if(weight==0) {
				observed.point[flag] <- 1/nHit
			}else if(weight==1){
				hit_tmp <- abs(rank.score.sorted[flag])
				observed.point[flag] <- hit_tmp/sum(hit_tmp)
			}else{
				hit_tmp <- abs(rank.score.sorted[flag] ** weight)
				observed.point[flag] <- hit_tmp/sum(hit_tmp)      
			}
			RES <- cumsum(observed.point)
			max.RES <- max(RES)
			min.RES <- min(RES)
			es.observed <- signif(ifelse(max.RES>abs(min.RES), max.RES, min.RES), digits=5)
			###########################################################
			if(0){
				## TODO: based on max.RES and min.RES
				es.position <- ifelse(max.RES>abs(min.RES), which.max(RES), which.min(RES))
			}else{
				## TODO: based on sign of ES
				if(res$ES[res$setID==x] > 0){
					es.position <- which.max(RES)
				}else{
					es.position <- which.min(RES)
				}
			}
			###########################################################
			## for leading genes
			if(RES[es.position]<0){
				ind <- which(flag >= es.position)
			}else{
				ind <- which(flag <= es.position)
			}
			hits <- rep(0, length(RES))
			hits[flag] <- 1
			hits[flag[ind]] <- 2
			hits[es.position] <- 3
			
			GeneID <- geneid.sorted
			ind <- match(GeneID, allGeneID)
			Symbol <- allSymbol[ind]
			df_leading <- data.frame(GeneID=GeneID, Rank=1:length(RES), Score=rank.score.sorted, RES=RES, Hits=hits, Symbol=Symbol, stringsAsFactors=FALSE)
			
		})
		names(ls_df_leading) <- res$setID
		
		ls_leadingGenes <- lapply(ls_df_leading, function(x){
			x$Symbol[x$Hits>=2]
		})
		
    }
	
	# leadingGenes
	if(1){
		
		leadingGenes <- lapply(ls_leadingGenes,function(x){
			ind <- match(x, rownames(df_priority))
			rank <- df_priority$rank[ind]
			names(rank) <- x
			return(rank)
		})
		
		## append leading genes
		res$nLead <- sapply(leadingGenes,length)
		## append peak
		res$peak <- sapply(leadingGenes,max)
	}
	
	## append "term_name" and "term_distance"
	ind <- match(res$setID, V(g)$name)
	summary <- data.frame(setID=res$setID, name=V(g)$term_name[ind], nAnno=res$setSize, nLead=res$nLead, peak=res$peak, total=rep(nrow(df_priority),length(res$peak)), es=res$ES, nes=res$nES, pvalue=res$pvalue, adjp=res$adjp, distance=V(g)$term_distance[ind], stringsAsFactors=FALSE)
	
	## scientific notation
	summary$es <- signif(summary$es, digits=3)
	summary$nes <- signif(summary$nes, digits=3)
	summary$pvalue <- signif(summary$pvalue, digits=2)
	summary$pvalue <- ifelse(summary$pvalue<0.01 & summary$pvalue!=0, as.numeric(format(summary$pvalue,scientific=TRUE)), summary$pvalue)
	summary$adjp <- signif(summary$adjp, digits=2)
	summary$adjp <- ifelse(summary$adjp<0.01 & summary$adjp!=0, as.numeric(format(summary$adjp,scientific=TRUE)), summary$adjp)
	
    ################################
    cross <- matrix(0, nrow=length(leadingGenes), ncol=length(leadingGenes))
    if(length(leadingGenes)>=2){
		for(i in seq(1,length(leadingGenes)-1)){
			x1 <- names(leadingGenes[[i]])
			for(j in seq(i+1,length(leadingGenes))){
				x2 <- names(leadingGenes[[j]])
				cross[i,j] <- length(intersect(x1, x2))
				cross[j,i] <- length(intersect(x1, x2))
			}
		}
		colnames(cross) <- rownames(cross) <- names(leadingGenes)
		diag(cross) <- sapply(leadingGenes, length)
    }
    ####################################################################################
	
    eGSEA <- list(df_summary = summary,
    			  leading = leadingGenes,
    			  full = ls_df_leading,
    			  cross = cross
                 )
    class(eGSEA) <- "eGSEA"
    
####################################################################################
    endT <- Sys.time()
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    
    if(!silent){
    	message(paste(c("\nEnd at ",as.character(endT)), collapse=""), appendLF=TRUE)
    	message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=TRUE)
    }
    
    invisible(eGSEA)
}
