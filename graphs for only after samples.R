#########################
# After + Donor samples #
#########################

####################################
# filtering and preparing the data #
####################################

#filtering taxonomy files
psonly_after = subset_samples(psfiltered, (Timepoint == "Post-FMT" | Timepoint == "Donor"))
filter <- phyloseq::genefilter_sample(psonly_after, filterfun_sample(function(x) x >= 100), A = 0.2*nsamples(psonly_after))
psonly_after = prune_taxa(filter, psonly_after)
psonly_after_norm = phyloseq_standardize_otu_abundance(psonly_after, method = "total")

#extracting all after samples from functional
subsetColData <- allColData_filtered_samples[grep("Post-FMT|Donor", allColData_filtered_samples$Timepoint), ]
subsetSamples <- as.vector(rownames(subsetColData))
subsetCounts <- allCountData[, subsetSamples]

#creating DEseq objects - don't use normalized data

#taxonomy

psonly_after_for_DE = subset_samples(psfiltered, (Timepoint == "Post-FMT"))
filter <- phyloseq::genefilter_sample(psonly_after_for_DE, filterfun_sample(function(x) x >= 100), A = 0.2*nsamples(psonly_after_for_DE))
psonly_after_for_DE = prune_taxa(filter, psonly_after_for_DE)

DES_psonly_after = phyloseq_to_deseq2(psonly_after_for_DE, ~ Outcome)
DES_psonly_after = DESeq(DES_psonly_after)


#functional
countdata.deseq <- DESeq2::DESeqDataSetFromMatrix(countData = subsetCounts, 
                                                  colData = subsetColData, 
                                                  design = ~ Outcome)
#filtering the data
filter <- rowSums(counts(countdata.deseq) >= 4) >= (0.2*nsamples(psonly_after))
countdata.deseq <- countdata.deseq[filter, ]




 
#######################
#     Oridinations    #
#######################

#ordinate normalized filtered data by method and distance
PS.ord <- ordinate(psonly_after_norm, method = "PCoA",distance =  "bray")

percentVar <- round(100 * PS.ord$values$Relative_eig, digits = 1)

#antibiotic treatment taxonomy


 
plot_ordination(psonly_after_norm, PS.ord, type="samples", color="Antibiotic_treatment_1_month") +
  coord_fixed() +
  #ggtitle("antibiotics 1 month - only after") +
  geom_point(size = 2) +
  stat_ellipse(geom = "polygon", type="norm", alpha=0.1, aes(fill=Antibiotic_treatment_1_month)) +

  #geom_text(aes(label = Case_study), size = 4, vjust = 1) + 
  scale_color_manual(values = c("Donor" = "blue", "No" = "green3", "Yes" = "red"),labels = c("Donors","No", "Yes"))+
  scale_fill_manual(values = c("Donor" = "blue", "No" = "green3", "Yes" = "red"),labels = c("Donors","No", "Yes"))+
  geom_text_repel(aes(label = (ifelse(Case.study == "Case 8 ", "Case 8", NA))), 
                  color = "black", size = 3,  show.legend = FALSE,
                  nudge_y = -0.6, nudge_x = -0.4, segment.size = 0.2, hjust = 0.5, direction = "x")+
  xlab(paste0("PC1: ", percentVar[1], "%")) + 
  ylab(paste0("PC2: ", percentVar[2], "%")) +
  labs(color = "Antibiotic treatment", fill = "Antibiotic treatment") +
  theme_bw()


#outcome taxonomy

plot_ordination(psonly_after_norm, PS.ord, type="samples", color="Outcome") +
  coord_fixed() +
  #ggtitle("Outcome - only after") +
  geom_point(size = 2) +
  stat_ellipse(geom = "polygon", type="norm", alpha=0.1, aes(fill=Outcome)) +
  scale_color_manual(values = c("Donor" = "blue", "Responders" = "green3", "Non-responders" = "red"),
                     labels = c("Donors","Failures", "Responders"))+
  scale_fill_manual(values = c("Donor" = "blue", "Responders" = "green3", "Non-responders" = "red"),
                    labels = c("Donors","Failures", "Responders"))+
  geom_text_repel(aes(label = (ifelse(Case.study == "Case 8 ", "Case 8", NA))), 
                  color = "black", size = 3,  show.legend = FALSE,
                  nudge_y = -0.5, nudge_x = -0.3, segment.size = 0.2, hjust = 0.5, direction = "x")+
  xlab(paste0("PC1: ", percentVar[1], "%")) + 
  ylab(paste0("PC2: ", percentVar[2], "%")) +
  theme_bw()
 



#outcome functional

vsd.deseq <- varianceStabilizingTransformation(countdata.deseq, blind = TRUE)
pcaData <- plotPCA(vsd.deseq, intgroup = c( "Outcome", "Timepoint"), returnData = TRUE)
pcaData 
percentVar <- round(100 * attr(pcaData, "percentVar"),digits = 1)

#sampleDists <- dist(t(assay(vsd.deseq)), method = "euclidean")

#adonis(sampleDists ~ Outcome, data = metadata)




ggplot(pcaData, aes(x = PC1, y = PC2, color = Outcome)) +
  coord_fixed()+
  #ggtitle("Post-FMT with donors functional (Modules) \n Outcome") +
  #geom_text(label = subsetColData$Case_study, vjust = 1) +
  geom_point(size = 2) + #the size of the points on the graph 
  xlab(paste0("PC1: ", percentVar[1], "%")) + 
  ylab(paste0("PC2: ", percentVar[2], "%")) +
  stat_ellipse(geom = "polygon", type="norm", alpha=0.1, aes(fill=Outcome)) +
  scale_colour_manual(values = c("Donor" = "blue", "Responders" = "green3", "Non-responders" = "red"),
                      labels = c("Donors","Failures", "Responders")) + 
  scale_fill_manual(values = c("Donor" = "blue", "Responders" = "green3", "Non-responders" = "red"),
                    labels = c("Donors","Failures", "Responders")) + 
  geom_text_repel(aes(label = (ifelse(rownames(pcaData) == "case8_after", "Case 8", NA))), 
                  color = "black", size = 2,  show.legend = FALSE,
                  nudge_y = 10, nudge_x = 8, segment.size = 0.2, hjust = 0.5, direction = "x")+
  theme_bw()





###########################################
#     differential analysis - functional  #
###########################################

#creating a DE analysis object
DE.countdata.deseq <- DESeq(countdata.deseq)
DEres <- results(DE.countdata.deseq, contrast=c("Outcome", "Responders", "Non-responders"))
DEres <- as.data.frame(DEres)

#volcano plot to visualize the differences between the factor we check
ggplot() +
  geom_point(data = DEres, mapping = aes(x = log2FoldChange, y = -log10(pvalue)), 
             col = "grey80", pch = 21, size = 4) + 
  ggtitle("Modules, FC>2, alpha<0.05") +
  geom_point(data = subset(DEres, log2FoldChange < -2 & padj < 0.05),
             aes(log2FoldChange, -log10(pvalue)), color = "red", size =4) +
  geom_point(data = subset(DEres, log2FoldChange > 2 & padj < 0.05),
             aes(log2FoldChange, -log10(pvalue)), color = "green3", size =4)+
  theme_bw()
  

#+
  #geom_text(data = subset(DEres, abs(log2FoldChange) > 2 & padj < 0.05),
        #  aes(log2FoldChange, -log10(pvalue)), label = rownames(subset(DEres, abs(log2FoldChange) > 2 & padj < 0.05)), vjust =1)

#normalizing counts
sizeFactors <- estimateSizeFactors(countdata.deseq)
normCounts.deseq <- counts(sizeFactors ,normalized=TRUE)

#subsetting the data to include only the genes with highest DE
highestDE <- subset(DEres, abs(log2FoldChange) > 2 & padj < 0.05)
nrow(highestDE)
highestDEnames <- rownames(highestDE)

#creating a vector to filter the donors out if needed
subsetDE2ColData <- allColData_filtered_samples[grep("Post-FMT|Donor", allColData_filtered_samples$Timepoint), ]
subsetDE2Samples <- as.vector(rownames(subsetDE2ColData))

#preparing the data for heatmap analysis, subsetting it by the most DE modules
highestDEexp <- as.data.frame(subset(normCounts.deseq, rownames(normCounts.deseq) %in% highestDEnames,
                                     colnames(normCounts.deseq) %in% subsetDE2Samples))
#num of clusters can be changed according to heatmap results
num_clusters <- 2

annotation_for_cols = data.frame(subsetDE2ColData[0], Group = subsetDE2ColData$group)
var1 = list(Group = c("Donor" = "royalblue1","Responders, post-FMT" = "palegreen","Non-responders, post-FMT" = "red1"))



heatmap.numClusters <- pheatmap::pheatmap(highestDEexp, 
                                          annotation_col = annotation_for_cols,
                                          annotation_colors = var1,
                                          scale = "row",
                                          show_rownames = T, 
                                          annotation_legend = F,
                                          treeheight_row = 0,
                                          labels_col = as.character(subsetDE2ColData$Case.study),
                                          #main = "Post-FMT with donors functional (modules) \n sig modules fc>4, alpha<0.05",
                                          border_color = "NA",
                                          cutree_cols = num_clusters,
                                          #cutree_rows = num_clusters,
                                          clustering_method = "complete",
                                          color = colorRampPalette(rev(RColorBrewer::brewer.pal(10, "RdYlBu")))(256))





#extracting to file the most diff expressed modules
alpha = 0.05
sigtab = DEres[which(DEres$padj < alpha), ]
#order results by pvalues
sigtab = sigtab[order(sigtab$pvalue),]
head(sigtab)
write.csv(sigtab,"only after/significant modules post-fmt.csv", row.names = TRUE)

#########################################
#     differential analysis - taxonomy  #
#########################################

res = results(DES_psonly_after,contrast=c("Outcome", "Responders", "Non-responders"))
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(psonly_after_for_DE)[rownames(sigtab), ], "matrix"))

#order results by pvalues
sigtab = sigtab[order(sigtab$pvalue),]
head(sigtab)

write.csv(sigtab,"corrected graphs/significant OTUs only after responders vs non.csv", row.names = TRUE)
