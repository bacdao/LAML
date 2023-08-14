#####Figure1

###R
##Loading R library
library(ggplot2)
library(data.table)
library(ggrepel)
library(readxl)
library(ggpubr)
library(viridis)
library(dplyr)
library(gghalves)
library(VennDiagram)

##-------------
##Loading data
##Clinical information obtaning from cBioportal 
clinical <- fread("~/Downloads/article/laml/data/laml_tcga/data_clinical_patient.txt", sep="\t", header=T)

##Labeling the 173 samples into 2 group 'FLT3 Mutant' and 'FLT3 wildtype'
designTSV <- fread("/home/bacdao/Downloads/article/laml/data/design.tsv", sep="\t", header=T)

# DEGs file from cbioportal
deg <- fread("/home/bacdao/Downloads/article/laml/data/deg.tsv", sep = "\t", header = T)
# dim(deg)
# [1] 19000   10

#Most anti-correlated probes file from cbioportal
acp <- fread("/home/bacdao/Downloads/article/laml/data/acp.tsv", header = T, sep = "\t")
# dim(acp)
# [1] 255  10

##Keeping the useful information
deg_filter <- data.frame(deg$Gene,deg$`Log Ratio`,deg$`q-Value`)
acp_filter <- data.frame(acp$Gene,acp$`Log Ratio`,acp$`q-Value`)

##Changing the colnames
colnames(deg_filter) <- c("gene","Log2FC","qvalue")
rownames(deg_filter) <- deg_filter$gene

colnames(acp_filter) <- c("gene","Delta","qvalue")


designTSV$Title <- substring(designTSV$Title, 1, 12)

##Only keeping 173 samples
clinical <- subset(clinical, clinical$`Patient Identifier` %in% designTSV$Title)

clinical_merged <- merge(clinical, designTSV, by.x = 'Patient Identifier', by.y = 'Title')
#---------------------
###Figure1.A
##Comparing age at diagnosis plot

colnames(clinical_merged) <- gsub('Diagnosis Age', 'Age', colnames(clinical_merged))

ggplot(clinical_merged, aes(x=Group, y=Age)) + theme_bw() + geom_half_point(aes(y=Age, fill=Group, color=Group), side="l", size=0.6, alpha=0.8) + geom_half_boxplot(aes(y=Age, color=Group, fill=Group), side="l", width=0.5, alpha=0.5, nudge=0.1, outlier.size=0) + geom_half_violin(aes(y=Age, fill=Group, color=Group), side="r") + stat_compare_means(label = "p.format", label.x = 1.4, size=10, vjust=0.7) + theme(text = element_text(size=20),axis.title.y= element_text(size=16, face="bold"), legend.position="top") + scale_x_discrete(labels=c("","")) + xlab("") + ylab("Age at diagnosis") + scale_fill_brewer(palette="Set1") + scale_color_brewer(palette="Set1")

#--------------------------------
##Comparing Gender percentage plot
ggplot(clinical_merged, aes(x = Group, fill = factor(Sex))) +
  geom_bar(position = "fill", alpha = 0.7) +
  scale_fill_brewer(palette="Dark2") +  # Set colors for females and males
  labs(x = "", y = "Gender", fill = "Group") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 16, face = "bold"),
    legend.position = "top",
    legend.text = element_text(size = 16)
  ) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = scales::percent(c(0, 0.25, 0.5, 0.75, 1)))

# Perform Chi-square test
chisq_result <- chisq.test(table(survival_merged$Group, survival_merged$Sex))
p_value <- chisq_result$p.value
# > p_value
# [1] 0.9475453

#------------------------------------------
###Comparing 'FAB' classification plot

ggplot(survival_merged, aes(x = Group, fill = factor(FAB))) +
  geom_bar(position = "fill", alpha = 0.7) +
  scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Accent"))(colourCount)) +  # Set colors for females and males
  labs(x = "", y = "FAB", fill = "Group") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 16, face = "bold"),
    legend.position = "top",
    legend.text = element_text(size = 16)
  ) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = scales::percent(c(0, 0.25, 0.5, 0.75, 1)))
dev.off()

fisher_result <- fisher.test(table(survival_merged$Group, survival_merged$FAB))                             
p_value <- fisher_result$p.value
# p_value
# [1] 0.6469185

#-----------------------
##Volcano plot

deg_filter$diffexpressed <- "NO"
deg_filter$diffexpressed[deg_filter$Log2FC >= 0.5 & deg_filter$qvalue < 0.05] <- "UP"
deg_filter$diffexpressed[deg_filter$Log2FC <= -0.5 & deg_filter$qvalue < 0.05] <- "DOWN"
deg_filter$diffexpressed[abs(deg_filter$Log2FC) < 0.5 & deg_filter$qvalue < 0.05] <- "FDR"
deg_filter$label <- ifelse(deg_filter$diffexpressed != c("NO","FDR"), deg_filter$gene, "")

ggplot(data = deg_filter, aes(x = Log2FC, y = -log10(qvalue), col = diffexpressed)) +
  geom_vline(xintercept = c(-0.5, 0.5), col = "black", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "black", linetype = 'dashed') + 
  geom_point(size = 4, alpha = 0.4) + 
  scale_color_manual(
    values = c("blue", "green", "grey", "red"),
    labels = c("Downregulated", "FDR", "NS", "Upregulated")
  ) + theme_bw() +
  theme(
    text = element_text(size = 12),
    legend.position = "top",
    legend.box = "vertical",
    legend.text = element_text(size = 14),
    axis.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 16, face = "bold")
  ) +
  labs(x = "Log2FC", y = "-Log10 FDR", color = "") +
  guides(alpha = FALSE, col = guide_legend(override.aes = list(label = ""))) +
  geom_text_repel(data = subset(deg_filter, abs(Log2FC) >= 0.5 & qvalue < 0.05), aes(label = gene), color = "black")

deg_filter$Exp <- ifelse(deg_filter$Log2FC>=0, "Up-regulated genes", "Down-regulated genes")

table(subset(deg_filter, deg_filter$qvalue<0.05)$Exp)
# Down-regulated genes   Up-regulated genes 
#                 1816                 1020 

table(subset(deg_filter, deg_filter$qvalue<0.05 & abs(deg_filter$Log2FC) >=0.5)$Exp)
# Down-regulated genes   Up-regulated genes 
#                  318                  193

#----------------------------------------
### Venndigram
#Build function for Venn diagram in R
display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}
#Venn diagram for overlapping between significant DEGs (qvalue <0.05) and genes contains differentially methylated most anti-correlated probes
venn1 <- list(A=subset(deg_filter,deg_filter$qvalue<0.05)$gene, B=acp_filter$gene) 

display_venn(
  venn1,
  category.names = c("DEGs" , "ACPs"), lwd=1, lty ="blank",
  fill = c("purple", "red"), cex= 1.5, fontface= "italic", dist=c(0.05,0.05), cat.cex =1.5, cat.fontface="bold", cat.pos=c(0,0), cat.dist=c(0.05, 0.05) )

#-----------------------------------
###Bar chart
#Merge DEG with qvalue < 0.05 and acp
deg.acp <- merge(subset(deg_filter, deg_filter$qvalue<0.05), acp_filter, by.x="gene", by.y="gene")
table(deg.acp$Exp)
# Down-regulated genes   Up-regulated genes 
#                   34                   33 

deg.acp$Meth <- ifelse(deg.acp$Delta >=0, "Hyper-methylated probes", "Hypo-methylated probes")
table(deg.acp$Meth)
# Hyper-methylated probes  Hypo-methylated probes 
#                       7                      60 

ggplot(deg.acp, aes(x = Exp, fill = Meth)) +
  geom_bar(alpha = 0.7) +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  theme(
    text = element_text(size = 12),
    legend.position = "top",
    legend.box = "vertical",
    legend.text = element_text(size = 14),
    axis.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 16, face = "bold")
  ) +
  labs(x = "", y = "Number of genes", fill = "")

#-----------------------
###Scatter plot between Log2FC expression gene and Delta Beta value

deg.acp$rlogqvalue.x <- -log10(deg.acp$qvalue.x)
deg.acp$rlogqvalue.y <- -log10(deg.acp$qvalue.y)

ggplot(deg.acp, aes(x = Delta, y = Log2FC)) +
  geom_point(aes(color = rlogqvalue.x, size = rlogqvalue.y, alpha = 0.5)) +
  labs(
    y = "Log2FC",
    x = "Delta Beta value",
    color = "-log10 FDR DEGs",
    size = "-log10 FDR ACPs"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 12),
    legend.position = "right",
    legend.box = "vertical",
    legend.text = element_text(size = 14),
    legend.title = element_text(face = "bold"),
    axis.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 16, face = "bold")
  ) +
  geom_smooth(method = "lm") +
  stat_cor(
    aes(label = paste(tolower(..r.label..), ..p.label.., sep = "~`,`~")),
    size = 6,
    label.x = -0.09,
    label.y = 2
  ) +
  scale_color_gradientn(colours = rainbow(8)) +
  geom_text_repel(aes(label = gene)) + guides(alpha = FALSE)
