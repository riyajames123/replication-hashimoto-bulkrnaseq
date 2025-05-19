# =============================
# Functional Enrichment Analysis of PCD-Related DEGs
# =============================

# 1. Load Required Libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(enrichplot)
library(ggplot2)

# Optional: If not installed, install them with:
# BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "DOSE"))

# 2. Load Your DEG Table with Gene Symbols (from earlier)
results <- read.csv("DE_results_all.csv", row.names = 1)
pcd_degs <- read.csv("PCD_DEGs.csv")  # This must have a 'GeneSymbol' column

# 3. Map Gene Symbols to Entrez IDs (needed for GO/KEGG)
entrez_map <- bitr(pcd_degs$GeneSymbol,
                   fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)

# Merge to keep logFC and adjusted p-value
pcd_degs_entrez <- merge(pcd_degs, entrez_map, by.x = "GeneSymbol", by.y = "SYMBOL")
gene_ids <- unique(pcd_degs_entrez$ENTREZID)

# =============================
# 4. GO Enrichment (Biological Process)
# =============================
go_enrich <- enrichGO(gene = gene_ids,
                      OrgDb = org.Hs.eg.db,
                      keyType = "ENTREZID",
                      ont = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.2,
                      readable = TRUE)

# View top enriched terms
head(go_enrich)

# =============================
# 5. KEGG Pathway Enrichment
# =============================
kegg_enrich <- enrichKEGG(gene = gene_ids,
                          organism = "hsa",
                          keyType = "kegg",
                          pvalueCutoff = 0.05)

# View top pathways
head(kegg_enrich)

# =============================
# 6. Plotting Enrichment Results
# =============================

## GO Barplot
barplot(go_enrich,
        showCategory = 10,
        title = "GO BP Enrichment: PCD DEGs",
        font.size = 12)

## GO Dotplot
dotplot(go_enrich,
        showCategory = 10,
        title = "GO BP Enrichment (Dotplot)",
        font.size = 10)

## KEGG Dotplot
dotplot(kegg_enrich,
        showCategory = 10,
        title = "KEGG Pathway Enrichment",
        font.size = 10)

# =============================
# 7. Save Enrichment Tables
# =============================

write.csv(as.data.frame(go_enrich), "GO_BP_enrichment_PCD_DEGs.csv", row.names = FALSE)
write.csv(as.data.frame(kegg_enrich), "KEGG_enrichment_PCD_DEGs.csv", row.names = FALSE)
