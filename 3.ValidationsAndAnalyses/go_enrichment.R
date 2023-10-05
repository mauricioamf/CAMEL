library(readr)
library(clusterProfiler)
library(org.EcK12.eg.db)
library(org.Sc.sgd.db)

data <- read_csv("enrichment_lowhigh_ratio.csv")

eco_kegg <- enrichGO(
  data$Yeast_kcat_low,
  org.Sc.sgd.db,
  ont = "ALL",
  pvalueCutoff = 0.1,
  #qvalueCutoff = 0.25
)

results_DF = data.frame(eco_kegg)
ID_DF = data.frame(results_DF$ID,results_DF$Description,results_DF$p.adjust,results_DF$Count)

dotplot(eco_kegg, showCategory = 10, title = "Yeast_kcat_low", orderBy = "Count")
