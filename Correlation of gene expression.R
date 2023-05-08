# Load necessary Packages
library(ggplot2)
library(ggpubr)

#---Fig. 2f------------------------------------------------
# Load data of the gene expression (Log10(Counts+1)) of E. coli samples with and without rRNA depletion.
data<-read.csv('Gene expression-rRNA depletion.csv')

# Create a scatter plot.
ggscatter(C, x="rRNA-depletion",y="No-rRNA-depletion",
          color = "gray25",size = 1, font.label = c(20, "plain"),cor.coef.size =7# Plot style.
          add = "reg.line",# Add regressin line.
          conf.int = TRUE,# Add confidence interval.
          cor.coef = TRUE, # Add correlation coefficient.
          cor.method = "pearson",# Method for computing correlation coefficient.
          )
  theme_bw()+theme(text = element_text(size=20))
# Save plot.
ggsave("rRNA-correlation.png",width = 8,height = 7) 

#---Supplementary Fig. 9b-----------------------------------
# Load data of the gene expression (Log10(Counts+1)) of repeated E. coli samples.
data<-read.csv()

# Create a scatter plot.
ggscatter(C, x="repeat1",y="repeat2",
          color = "gray25",size = 1, font.label = c(20, "plain"),cor.coef.size =7# Plot style.
          add = "reg.line",# Add regressin line.
          conf.int = TRUE,# Add confidence interval.
          cor.coef = TRUE, # Add correlation coefficient.
          cor.method = "pearson",# Method for computing correlation coefficient.
          )
  theme_bw()+theme(text = element_text(size=20))
# Save plot.
ggsave("technical repeatability-correlation.png",width = 8,height = 7) 