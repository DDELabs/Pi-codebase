# Pi — Genetic Prioritisation of Drug Targets from GWAS SNPs

> Leverages GWAS SNPs to prioritise gene targets and pathways using network-based Random Walk with Restart (RWR).

---

## Overview

Pi takes a list of GWAS lead SNPs (with p-values) and, through a multi-step pipeline, produces a ranked list of gene targets, enriched pathways, and a network of top targets — all with performance evaluation.

---

## Pipeline at a Glance

```
[User Input: SNPs + p-values]
         │
         ▼
  xPierSNPsAdv          ← orchestrator; runs xPierSNPs once per predictor type
         │
         ├── xPierSNPs [distance predictor]    relative.importance = c(1,0,0)
         ├── xPierSNPs [eQTL predictor(s)]     relative.importance = c(0,1,0)
         └── xPierSNPs [HiC predictor(s)]      relative.importance = c(0,0,1)
                  │
                  ├── xSNPscores          (p-value → SNP score, LD expansion)
                  ├── xGR2xGeneScores     (SNPs → seed genes + scores)
                  │     ├── xGR2nGenes    (distance-based nearby genes)
                  │     ├── xDefineEQTL   (eQTL-linked genes)
                  │     ├── xDefineHIC    (Promoter Capture HiC-linked genes)
                  │     └── xAggregate    (aggregate scores per gene)
                  └── xPierGenes          (RWR on gene network)
                        ├── xDefineNet    (load STRING / PCommons network)
                        └── xPier → xRWR  (Random Walk with Restart)

         ▼
  returns: ls_pNode  (one pNode object per predictor)

         ▼
  xPierMatrix           ← aggregates ls_pNode into a priority matrix → dTarget

         ▼
  xPierPathways         ← enrichment analysis of top-ranked genes → eTerm
  xPierSubnet           ← maximum-scoring gene subnetwork → igraph

         ▼
  xPierROCR             ← ROC + Precision-Recall analysis, calculates AUC
  xPredictCompare       ← compares AUC against user-selected benchmarking method
```

---

## Input

A two-column data frame or named vector:

| Column | Content |
|--------|---------|
| SNP | dbSNP ID (e.g. `rs1234567`) or genomic position (`chr16:28525386`) |
| Pvalue | GWAS p-value or FDR for each lead SNP |

---

## Parameters

### SNP Scoring

| Parameter | What it controls | Example values |
|-----------|-----------------|----------------|
| `significance.threshold` | p-value cutoff; SNPs above this get zero score | `5e-5` |
| `LD.r2` | r² threshold for including LD SNPs | `0.8` (range: 0.8–1.0) |

### LD Population Panel

Choose one or more 1000 Genomes Project super-populations, **or** supply your own LD table.

| Code | Population |
|------|-----------|
| `AFR` | African |
| `AMR` | Admixed American |
| `EAS` | East Asian |
| `EUR` | European |
| `SAS` | South Asian |

For custom LD: provide `LD.customised` — a 3-column data frame (Lead SNP, LD SNP, r² value).

### Gene-SNP Linking

| Parameter | What it controls | Options |
|-----------|-----------------|---------|
| `distance.max` | Max kb distance for nearby genes | e.g. `2000` (2 Mb) |
| `decay.kernel` | How distance weight falls off | `slow`, `linear`, `rapid`, `constant` |
| `decay.exponent` | Exponent for the decay function | `2` (default) |
| `include.TAD` | TAD boundary filter for nearby genes | `none`, `GM12878`, `IMR90`, `MSC`, `TRO`, `H1`, `MES`, `NPC` |

### eQTL Dataset

| Parameter | Options |
|-----------|---------|
| `include.eQTL` | Pre-built context codes (see `xDefineEQTL`); e.g. `"JKng_mono"` |
| `eQTL.customised` | Your own 4-column table: SNP, Gene, Sig, Context — **or** a file path |

### HiC Dataset

| Parameter | Options |
|-----------|---------|
| `include.HiC` | Pre-built PCHiC cell type codes (see `xDefineHIC`); e.g. `"Monocytes"` |

### Network (for RWR)

| Parameter | Options |
|-----------|---------|
| `network` | `STRING_highest`, `STRING_high`, `STRING_medium`, `STRING_low`, `PCommonsUN_high`, `PCommonsUN_medium`, `PCommonsDN_high`, `PCommonsDN_medium`, `KEGG`, `REACTOME`, and more |
| `network.customised` | Your own `igraph` object (overrides built-in network) |

### RWR Parameters

| Parameter | What it controls | Default |
|-----------|-----------------|---------|
| `restart` | Restart probability (0–1); higher = closer to seeds | `0.7` |
| `normalise` | Adjacency matrix normalisation | `laplacian`, `row`, `column`, `none` |
| `normalise.affinity.matrix` | Output affinity matrix normalisation | `none`, `quantile` |

### Scoring Scheme

| Parameter | Options |
|-----------|---------|
| `cdf.function` | How SNP scores are transformed to gene weights | `empirical`, `exponential` |
| `scoring.scheme` | How multiple SNP scores are aggregated per gene | `max`, `sum`, `sequential` |

---

## Step-by-Step Usage

### Step 1 — Run xPierSNPsAdv

```r
RData.location <- "http://galahad.well.ox.ac.uk/bigdata"

ls_pNode <- xPierSNPsAdv(
  data                    = SNP_data,

  # SNP scoring
  significance.threshold  = 5e-5,

  # LD
  include.LD              = "EUR",       # or LD.customised = my_ld_df
  LD.r2                   = 0.8,

  # Gene-SNP distance
  distance.max            = 2000,
  decay.kernel            = "slow",
  decay.exponent          = 2,
  include.TAD             = "GM12878",

  # eQTL
  include.eQTL            = "JKng_mono", # or eQTL.customised = my_eqtl_df

  # HiC
  include.HiC             = "Monocytes",

  # Scoring
  cdf.function            = "empirical",
  scoring.scheme          = "max",

  # Network + RWR
  network                 = "PCommonsUN_medium",
  restart                 = 0.7,
  normalise               = "laplacian",
  normalise.affinity.matrix = "none",

  RData.location          = RData.location
)
# Returns: named list of pNode objects
# - "nGene_2000_slow"        (distance predictor)
# - "eGene_JKng_mono"        (eQTL predictor)
# - "cGene_Monocytes"        (HiC predictor)
```

### Step 2 — Aggregate into Priority Matrix

```r
dTarget <- xPierMatrix(
  list_pNode    = ls_pNode,
  displayBy     = "rank",          # or "score"
  RData.location = RData.location
)
# Returns: dTarget object with aggregated priority scores across all predictors
```

### Step 3 — Pathway Enrichment

```r
eTerm <- xPierPathways(
  pNode          = dTarget,
  priority.top   = 100,            # use top 100 prioritised genes
  ontology       = "MsigdbC2CPall", # or GOBP, KEGG, REACTOME, etc.
  RData.location = RData.location
)
xEnrichViewer(eTerm)
```

### Step 4 — Target Subnetwork

```r
subnet <- xPierSubnet(
  pNode          = dTarget,
  priority.quantile = 0.10,        # top 10% of genes as candidates
  subnet.size    = 30,             # desired number of nodes
  RData.location = RData.location
)
# Returns: igraph object of the maximum-scoring gene subnetwork
```

### Step 5 — Performance Evaluation (AUC)

```r
# Requires Gold Standard Positives (GSP) and Negatives (GSN)
pPerf <- xPierROCR(
  pNode  = dTarget,
  GSP    = my_GSP_genes,
  GSN    = my_GSN_genes
)
# Returns pPerf object with ROC curve and AUC
```

### Step 6 — Compare Against Benchmarking Method

```r
xPredictCompare(
  list_pPerf = list(Pi = pPerf, Baseline = baseline_pPerf),
  displayBy  = "ROC"
)
# Plots AUC comparison across methods
```

---

## What Each Function Does

| Function | Input | Output | Role |
|----------|-------|--------|------|
| `xPierSNPsAdv` | SNPs + p-values | `ls_pNode` (list) | Orchestrates all predictor runs |
| `xPierSNPs` | SNPs + params | `pNode` | One predictor: score SNPs → seed genes → RWR |
| `xSNPscores` | SNPs | Scored GRanges | p-value → score, LD expansion |
| `xGR2xGeneScores` | Scored SNPs | Seed gene scores | Maps SNPs to genes (distance/eQTL/HiC) |
| `xDefineEQTL` | SNPs | eGene-SNP pairs | Loads eQTL data |
| `xDefineHIC` | SNPs | cGene-SNP pairs | Loads PCHiC data |
| `xPierGenes` | Seed genes + scores | `pNode` | Runs RWR on network |
| `xDefineNet` | Network name | `igraph` | Loads gene interaction network |
| `xRWR` | Graph + seeds | Affinity scores | Core Random Walk with Restart |
| `xPierMatrix` | `ls_pNode` | `dTarget` | Aggregates all predictors into priority matrix |
| `xPierPathways` | `dTarget` | `eTerm` | Enrichment of top genes in pathways |
| `xPierSubnet` | `dTarget` | `igraph` | Maximum-scoring gene subnetwork |
| `xPierROCR` | `dTarget` + GSP/GSN | `pPerf` | ROC/PR curves and AUC |
| `xPredictCompare` | list of `pPerf` | Plot | AUC comparison across methods |
| `xRDataLoader` | RData name | R object | Fetches built-in data from remote server |

---

## Output Objects

**`pNode`** (one per predictor, inside `ls_pNode`):
- `priority` — gene × 6 matrix: name, node, seed, weight, priority score [0,1], rank
- `g` — the igraph network used
- `SNP` — input SNPs + LD SNPs with scores
- `Gene2SNP` — gene–SNP pairs with influence scores
- `nGenes` — nearby gene info (distance predictor)
- `eGenes` — eQTL gene info per context
- `cGenes` — HiC gene info per context

**`dTarget`** (from `xPierMatrix`):
- `priority` — aggregated gene priority with rank, p-value, FDR, star rating
- `predictor` — full predictor matrix (one column per `pNode`)

---

## Data Sources (fetched at runtime)

Large reference files are **not stored in this repository**. They are downloaded automatically by `xRDataLoader` from:

```
http://galahad.well.ox.ac.uk/bigdata
```

| Data | Used by |
|------|---------|
| `dbSNP_GWAS` — SNP genomic coordinates | `xSNPscores` |
| `UCSC_knownGene` — gene coordinates | `xGR2nGenes` |
| LD panels (1000 Genomes, per population) | `xSNPscores` |
| eQTL datasets (GTEx, Blueprint, etc.) | `xDefineEQTL` |
| Promoter Capture HiC datasets | `xDefineHIC` |
| TAD boundary datasets | `xGR2nGenes` |
| Gene networks (STRING, PCommons, KEGG) | `xDefineNet` |

---

## Installation

```r
# Stable release (Bioconductor)
BiocManager::install("Pi")

# Latest development version
BiocManager::install(c("Pi", "devtools"), dependencies = TRUE)
BiocManager::install("hfang-bristol/Pi")
```

**Dependencies:** R (≥ 3.5), Bioconductor packages: `GenomicRanges`, `BiocGenerics`, `igraph`, `dnet`

---

## Citation

Fang H et al. *A genetics-led approach defines the drug target landscape of 30 immune-related traits.* Nature Genetics, 2019.
