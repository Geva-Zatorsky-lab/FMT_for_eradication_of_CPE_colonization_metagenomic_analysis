

res <- as.data.frame(res)

#volcano plot to visualize the differences between the factor we check
ggplot() +
  geom_point(data = res, mapping = aes(x = log2FoldChange, y = -log10(pvalue)), 
             col = "grey80", pch = 21, size = 4) + 
  ggtitle("taxa post-fmt non responders vs responders, \n FC>2, alpha<0.05") +
  geom_point(data = subset(res, log2FoldChange < -2 & padj < 0.05),
             aes(log2FoldChange, -log10(pvalue)), color = "red", size =4)  +
  geom_point(data = subset(res, log2FoldChange > 2 & padj < 0.05),
             aes(log2FoldChange, -log10(pvalue)), color = "green3", size =4)  +
  geom_text(data = subset(res, abs(log2FoldChange) > 2 & padj < 0.05),
          aes(log2FoldChange, -log10(pvalue)), label = rownames(subset(res, abs(log2FoldChange) > 2 & padj < 0.05)), vjust =1)

#normalizing counts
sizeFactors <- estimateSizeFactors(DES_psonly_after)
normCounts.deseq <- counts(sizeFactors ,normalized=TRUE)

#subsetting the data to include only the genes with highest DE
highestDE <- subset(res, abs(log2FoldChange) > 1 & padj < 0.05)
nrow(highestDE)
highestDEnames <- rownames(highestDE)

allColData_filtered_samples = meta(psonly_after)

#creating a vector to filter the donors out
subsetDE2ColData <- allColData_filtered_samples[ grep("Post-FMT|Donor", allColData_filtered_samples$Timepoint), ]
subsetDE2Samples <- as.vector(rownames(subsetDE2ColData))

#preparing the data for heatmap analysis, subsetting it by the most DE modules
highestDEexp <- as.data.frame(subset(normCounts.deseq, rownames(normCounts.deseq) %in% highestDEnames,
                                     colnames(normCounts.deseq) %in% subsetDE2Samples))


highestDEexp_Names <- read.csv(file = "significant taxa highestDE post-fmt NAMES.csv", row.names = 1)

Family_annotation <- read.csv(file = "sig_taxa_annotation_family.csv", row.names = 1)

#num of clusters can be changed according to heatmap results
num_clusters <- 2

annotation_for_cols = data.frame(subsetDE2ColData[0], Group = subsetDE2ColData$group)
var1 = list(Group = c("Donor" = "royalblue1","Responders, post-FMT" = "palegreen","Non-responders, post-FMT" = "red1"))


annotation_for_rows = data.frame(Family_annotation[0], Taxa = Family_annotation$Family_annotation_lite)




heatmap.numClusters <- pheatmap::pheatmap(highestDEexp_Names, 
                                          scale = "row",
                                          cluster_rows = T,
                                          annotation_col = annotation_for_cols,
                                          #annotation_row = annotation_for_rows,
                                          annotation_colors = var1,
                                          annotation_legend = F,
                                          treeheight_row = 0,
                                          show_rownames = T, 
                                          #main = "Pre-FMT with donors functional (modules) \n sig modules fc>4, alpha<0.05",
                                          labels_col = as.character(subsetDE2ColData$Case.study),
                                          fontsize_row = 10,
                                          border_color = "NA",
                                          #cutree_cols = num_clusters,
                                          #cutree_rows = num_clusters,
                                          clustering_method = "complete",
                                          color = colorRampPalette(rev(RColorBrewer::brewer.pal(10, "RdYlBu")))(256))




#ps_entero_sub = subset_taxa(psonly_after_for_DE, Family == "Enterobacteriaceae")

#plot_heatmap(ps_entero_sub, method= "PCoA", ditance ="bray", sample.order = "group", taxa.label="Species")


