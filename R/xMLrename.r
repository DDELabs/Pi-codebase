#' Function to rename predictors used in machine learning
#'
#' \code{xMLrename} is supposed to rename predictors used in machine learning. It returns an object of class "sTarget".
#'
#' @param sTarget an object of class "sTarget"
#' @param old_names a vector for the original names of predictors to be renamed
#' @param new_names a vector for the new names
#' @return an object of class "sTarget"
#' @note none
#' @export
#' @seealso \code{\link{xMLrename}}
#' @include xMLrename.r
#' @examples
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' \dontrun{
#' old_names <- colnames(sTarget$predictor)[-c(1,2)]
#' old_names_1 <- c('nGene_20000_constant', 'eGene_CD14', 'eGene_LPS2', 'eGene_LPS24', 'eGene_IFN', 'eGene_Bcell', 'eGene_CD4', 'eGene_CD8', 'eGene_Neutrophil', 'eGene_NK', 'eGene_Blood', 'cGene_Monocytes', 'cGene_Macrophages_M0', 'cGene_Macrophages_M1', 'cGene_Macrophages_M2', 'cGene_Neutrophils', 'cGene_Naive_CD4_T_cells', 'cGene_Total_CD4_T_cells', 'cGene_Naive_CD8_T_cells', 'cGene_Total_CD8_T_cells', 'cGene_Naive_B_cells', 'cGene_Total_B_cells', 'dGene', 'pGene', 'fGene')
#' old_names_2 <- c('nGene_20000_constant', 'eGene_Pi_eQTL_CD14', 'eGene_Pi_eQTL_LPS2', 'eGene_Pi_eQTL_LPS24', 'eGene_Pi_eQTL_IFN', 'eGene_Pi_eQTL_Bcell', 'eGene_Pi_eQTL_CD4', 'eGene_Pi_eQTL_CD8', 'eGene_Pi_eQTL_Neutrophil', 'eGene_Pi_eQTL_NK', 'eGene_Pi_eQTL_Blood', 'cGene_Monocytes', 'cGene_Macrophages_M0', 'cGene_Macrophages_M1', 'cGene_Macrophages_M2', 'cGene_Neutrophils', 'cGene_Naive_CD4_T_cells', 'cGene_Total_CD4_T_cells', 'cGene_Naive_CD8_T_cells', 'cGene_Total_CD8_T_cells', 'cGene_Naive_B_cells', 'cGene_Total_B_cells', 'dGene', 'pGene', 'fGene')
#' new_names <- c('nearbyGenes_nGene: nearby genes', 'eQTL_eGene: resting state (CD14+)', 'eQTL_eGene: activating state (CD14+ by LPS2h)', 'eQTL_eGene: activating state (CD14+ by LPS24h)', 'eQTL_eGene: activating state (CD14+ by IFN24h)', 'eQTL_eGene: B cells', 'eQTL_eGene: CD4+ T cells', 'eQTL_eGene: CD8+ T cells', 'eQTL_eGene: neutrophils', 'eQTL_eGene: NK cells','eQTL_eGene: peripheral blood', 'HiC_cGene: monocytes', 'HiC_cGene: macrophages (M0)', 'HiC_cGene: macrophages (M1)', 'HiC_cGene: macrophages (M2)', 'HiC_cGene: neutrophils', 'HiC_cGene: CD4+ T cells (naive)', 'HiC_cGene: CD4+ T cells (total)', 'HiC_cGene: CD8+ T cells (naive)', 'HiC_cGene: CD8+ T cells (total)', 'HiC_cGene: B cells (naive)', 'HiC_cGene: B cells (total)', 'Annotation_dGene: disease genes', 'Annotation_pGene: phenotype genes', 'Annotation_fGene: function genes')
#' sTarget_renamed <- xMLrename(sTarget, old_names_1, new_names)
#' sTarget_renamed <- xMLrename(sTarget_renamed, old_names_2, new_names)
#' }

xMLrename <- function(sTarget, old_names, new_names) 
{
    
    if(!is(sTarget,"sTarget")){
    	stop("The function must apply to a 'sTarget' object.\n")
    }
    sTarget_renamed <- sTarget
    
    ## predictor
	if(!is.null(sTarget_renamed$predictor)){
		ind <- match(colnames(sTarget_renamed$predictor), old_names)
		colnames(sTarget_renamed$predictor)[!is.na(ind)] <- new_names[ind[!is.na(ind)]]
	}
	
	## importance2fold
	if(!is.null(sTarget_renamed$importance2fold)){
		ind <- match(rownames(sTarget_renamed$importance2fold), old_names)
		rownames(sTarget_renamed$importance2fold)[!is.na(ind)] <- new_names[ind[!is.na(ind)]]
	}
		
	## roc2fold
	if(!is.null(sTarget_renamed$roc2fold)){
		ind <- match(rownames(sTarget_renamed$roc2fold), old_names)
		rownames(sTarget_renamed$roc2fold)[!is.na(ind)] <- new_names[ind[!is.na(ind)]]
	}
	
	## fmax2fold
	if(!is.null(sTarget_renamed$fmax2fold)){
		ind <- match(rownames(sTarget_renamed$fmax2fold), old_names)
		rownames(sTarget_renamed$fmax2fold)[!is.na(ind)] <- new_names[ind[!is.na(ind)]]
	}
	
	## importance
	if(!is.null(sTarget_renamed$importance)){
		ind <- match(rownames(sTarget_renamed$importance), old_names)
		rownames(sTarget_renamed$importance)[!is.na(ind)] <- new_names[ind[!is.na(ind)]]
	}
	
	## performance
	if(!is.null(sTarget_renamed$performance)){	
		ind <- match(rownames(sTarget_renamed$performance), old_names)
		rownames(sTarget_renamed$performance)[!is.na(ind)] <- new_names[ind[!is.na(ind)]]
	}
	
	invisible(sTarget_renamed)
}
