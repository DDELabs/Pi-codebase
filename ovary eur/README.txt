```{r}
data_lead <- read.table("lead_eur.txt", header = FALSE, sep = "\t")
str(data_lead)
# data_lead <- data_lead[data_lead$V2 <= 5e-8, ]
# str(data_lead)
data_ld <- read.table("ld_eur.txt", header = FALSE, sep = "\t")

load("n_GR.SNP_all.RData")
GR.SNP <- n_GR.SNP
load("TAD_ovary.RData")
include.TAD <- n_TAD
load("ovary_e097.RData")
n_chrom <- n_chrom

eqtl_c <- read.table("QTL_ov_eur_n.txt", header = F)

df_nodes_data <- read.table("df_nodes_ov.txt", sep="\t", header =T)
df_edges_data <- read.table("df_edges_ov.txt", sep="\t", header =T)
HiC <- graph_from_data_frame(d=df_edges_data, vertices=df_nodes_data, directed=TRUE)

# edges <- read.table("edges.txt", header = TRUE, sep = "\t")
# vertices <- read.delim("nodes.txt", header = TRUE, sep = "\t")
# g_n <- graph_from_data_frame(d = edges, vertices = vertices, directed = FALSE)
load("bg_f_wo_ubc.RData")
g_n <- bg
```

```{r}
source("xSNPscores.r")
source("xSNP2nGenes_w_chromhmm_w_tad.r")
source("xSNP2eGenes.r")
source("xDefineHIC.r")
source("xSNP2cGenes.r")
source("xPierSNPs.r")
source("xPierSNPsAdv.r")

ls_pNode_g <- PierSNPsAdv(data=data_lead, significance.threshold=5e-5, LD.customised=data_ld, LD.r2=0.5, distance.max=10000, decay.kernel="linear", decay.exponent=1, GR.SNP=GR.SNP, include.TAD=n_TAD, eQTL.customised=eqtl_c, cdf.function="empirical", scoring.scheme="max", include.HiC="Ovary", network.customised=g_n, normalise="laplacian", restart=0.6, normalise.affinity.matrix="none", RData.location=RData.location, verbose.details=T)
# str(ls_pNode_g)

save(ls_pNode_g, file = "ls_pNode_g_ov_eur_wo_ubc_n2_verify.RData")
```

```{r}
source("xPierAnno.r")
iA <- read.delim("i_anno.txt", header=TRUE, stringsAsFactors=FALSE)[,c("Symbol","Disease","Phenotype","Function")]
colnames(iA) <- c("Symbol","dGene","pGene","fGene")
ls_pNode_anno <- lapply(2:4, function(j){
    print(j)
    data_anno <- subset(data.frame(seed=iA$Symbol,weight=iA[,j],stringsAsFactors=F), weight>0)
    pNode <- xPierAnno(data_anno, list_pNode=ls_pNode_g, network.customised=g_n, seeds.inclusive=TRUE, normalise="laplacian", restart=0.6, normalise.affinity.matrix="none", RData.location=RData.location)
})
names(ls_pNode_anno) <- colnames(iA)[2:4]
# str(ls_pNode_anno)

save(ls_pNode_anno, file = "ls_pNode_anno_ov_eur_wo_ubc_n2_8_0.8.RData")
```

```{r}
ls_pNode <- c(ls_pNode_anno, ls_pNode_g)
save(ls_pNode, file = "ls_pNode_eur_wo_ubc_n2_8_0.8.RData")
# str(ls_pNode)
```

```{r}
source("xPierMatrix.r")

dTarget <- xPierMatrix(ls_pNode, displayBy="pvalue", aggregateBy="fishers", RData.location=RData.location)
class(dTarget)
# str(dTarget)

save(dTarget, file = "dTarget_ov_eur_wo_ubc_n2_8_0.5_verify.RData")

write.table(dTarget$priority, file="Pi_output_ov_eur_wo_ubc_n2_8_0.5_verify.txt", sep="\t", quote = F, row.names=FALSE)
```