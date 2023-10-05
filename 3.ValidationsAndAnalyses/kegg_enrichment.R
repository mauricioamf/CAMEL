library(readr)
library(clusterProfiler)

data <- read_csv("enrichment_lowhigh_ratio.csv")

eco_kegg <- enrichKEGG(
  data$Yeast_kapp_low,
  organism = "sce",
  keyType = "ncbi-geneid",
  pvalueCutoff = 0.1,
  #qvalueCutoff = 0.25
)

results_DF = data.frame(eco_kegg)
ID_DF = data.frame(results_DF$ID,results_DF$Description,results_DF$p.adjust,results_DF$Count)

dotplot(eco_kegg, showCategory = 10, title = "Yeast_kapp_low", orderBy = "Count")
