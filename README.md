# Pi — Genetic Prioritisation of Drug Targets from GWAS SNPs

> Leverages GWAS SNPs to prioritise gene targets and pathways using network-based Random Walk with Restart (RWR).
(https://github.com/r-forge/pi314/tree/master/pkg; https://github.com/hfang-bristol/PIC2; https://github.com/hfang-bristol/supraHex)
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
| `LD.r2` | r² threshold for including LD SNPs | `0.8` (range: 0–1.0) |

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

## Citation

Fang H et al. *A genetics-led approach defines the drug target landscape of 30 immune-related traits.* Nature Genetics, 2019.
