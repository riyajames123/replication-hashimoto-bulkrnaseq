# Replication Study: Bulk RNA-Seq Analysis in Hashimoto's Thyroiditis

This repository contains shell and R scripts used to replicate the upstream bulk RNA-seq analysis pipeline from the study:

> Yang, D., Wang, X., Sun, Y., Shao, Y., & Shi, X. (2024). Identification and experimental validation of genes associated with programmed cell death in dendritic cells of the thyroid tissue in Hashimoto’s thyroiditis. *International Immunopharmacology*, 142, 113083.

---

## 📦 Dataset

- RNA-seq samples from 10 human thyroid tissue samples (5 HT, 5 healthy controls)
- 8 female, 2 male
- Downloaded from NGDC under accession number **HRA001684**
- **Raw FASTQ files are NOT shared in this repository** due to size and data-sharing restrictions.

---

## 🔁 Workflow Summary

| Step | Description | Script |
|------|-------------|--------|
| **1. QC** | Adapter & quality trimming with Fastp | *not included* |
| **2. rRNA Removal** | Used Bowtie2 (unsuccessful), then Ribodetector (5–6% rRNA removed) | `scripts/ribo_depletion.sh` |
| **3. Alignment** | Aligned to GRCh38 with HISAT2 | `scripts/alignment_hisat2.sh` |
| **4. Quantification** | Gene-level quantification via FeatureCounts | `scripts/quantification_featurecounts.sh` |
| **5. DE Analysis** | `limma + voom` for normalization, DEG detection | `scripts/limma_de.R` |
| **6. GSEA & Marker Check** | GSEA + marker genes for dendritic cell check | `scripts/gse.R` |

---

## 🧬 DE Analysis Insights

- Used **TMM normalization** & `voom` transformation.
- PCA and MDS plots showed no sex-based separation — so sex was not included in the design matrix.
- DEGs were selected with:  
  `adjusted p-value < 0.05` and `|log2FC| > 1`
- DEGs were cross-referenced with:
  - PCD-related genes (from paper's supplementary)
  - DC markers (CD83, ITGAX), which were not significant.

---

## 📊 Interpretation

- Significant DEGs (e.g., **STAT1, TNFAIP3**) are known to affect dendritic cell function.
- No strong evidence of DC-specific markers differentially expressed — likely due to:
  - Small sample size (n=10)
  - Lack of cell-type resolution in bulk RNA-seq

---

## 🚫 Data Disclaimer

- Raw sequencing data was downloaded from **NGDC (HRA001684)**
- All scripts are shared, but raw data is **not uploaded** due to data size and policy.

---

## 📎 License

MIT License 

