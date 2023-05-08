
# Loading necessary Packages
library("clusterProfiler")
library("topGO")
library(org.EcK12.eg.db) 

#------Figure 3--------------------
# Read the csv file of the top20 DEGs in CIP-treatment experiment.
data<-read.csv('CIP-degs.csv')

# GO enrichment analysis of the 2h time group in CIP-treatment experiment.
ego <- enrichGO(
  gene = data$CIP_2h_n,
  OrgDb = org.EcK12.eg.db,
  keyType = 'SYMBOL',
  ont = "ALL", 
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05, 
  qvalueCutoff = 0.05)
dotplot(ego,font.size=12)+theme_classic()#Visualization of Go enrichment result

#------Figure 4--------------------
# Read the csv file of the top20 DEGs in AMP-treatment experiment.
data<-read.csv('AMP-degs.csv')

# GO enrichment analysis of the subcluster 6 in AMP-treatment experiment.
ego <- enrichGO(
  gene = data$X6_n,
  OrgDb = org.EcK12.eg.db,
  keyType = 'SYMBOL',
  ont = "ALL", 
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05, 
  qvalueCutoff = 0.05)
dotplot(ego,font.size=12)+theme_classic()#Visualization of Go enrichment result
# GO enrichment analysis of the subcluster 9 in AMP-treatment experiment.
ego <- enrichGO(
  gene = data$X9_n,
  OrgDb = org.EcK12.eg.db,
  keyType = 'SYMBOL',
  ont = "ALL", 
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05, 
  qvalueCutoff = 0.05)
dotplot(ego,font.size=12)+theme_classic()#Visualization of Go enrichment result
# GO enrichment analysis of the subcluster 10 in AMP-treatment experiment.
ego <- enrichGO(
  gene = data$X10_n,
  OrgDb = org.EcK12.eg.db,
  keyType = 'SYMBOL',
  ont = "ALL", 
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05, 
  qvalueCutoff = 0.05)
dotplot(ego,font.size=12)+theme_classic()#Visualization of Go enrichment result


