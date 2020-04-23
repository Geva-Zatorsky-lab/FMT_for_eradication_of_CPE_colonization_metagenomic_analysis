##################
# Before samples #
##################

####################################
# filtering and preparing the data #
####################################

#filtering taxonomy files
psonly_before = subset_samples(psfiltered, (Timepoint == "Pre-FMT"))
filter <- phyloseq::genefilter_sample(psonly_before, filterfun_sample(function(x) x >= 100), A = 0.2*nsamples(psonly_before))
psonly_before = prune_taxa(filter, psonly_before)
psonly_before_norm = phyloseq_standardize_otu_abundance(psonly_before, method = "total")

#extracting all after samples from functional
subsetColData <- allColData_filtered_samples[grep("Pre-FMT", allColData_filtered_samples$Timepoint), ]
subsetSamples <- as.vector(rownames(subsetColData))
subsetCounts <- allCountData[, subsetSamples]

#creating DEseq objects - don't use normalized data

#taxonomy
psonly_before_for_DE = subset_samples(psfiltered, (Timepoint == "Pre-FMT"))
filter <- phyloseq::genefilter_sample(psonly_before_for_DE, filterfun_sample(function(x) x >= 100), A = 0.2*nsamples(psonly_before_for_DE))
psonly_before_for_DE = prune_taxa(filter, psonly_before_for_DE)

DES_psonly_before = phyloseq_to_deseq2(psonly_before_for_DE, ~ Outcome)
DES_psonly_before = DESeq(DES_psonly_before)



#functional
countdata.deseq <- DESeq2::DESeqDataSetFromMatrix(countData = subsetCounts, 
                                                  colData = subsetColData, 
                                                  design = ~ Outcome)

#filtering the data
filter <- rowSums(counts(countdata.deseq) >= 4) >= (0.2*nsamples(psonly_before))
countdata.deseq <- countdata.deseq[filter, ]
 
#######################
#     Oridinations    #
#######################

#ordinate normalized filtered data by method and distance
PS.ord <- ordinate(psonly_before_norm, method = "PCoA",distance =  "bray")

percentVar <- round(100 * PS.ord$values$Relative_eig, digits = 1)

#antibiotic treatment taxonomy
 
plot_ordination(psonly_before_norm, PS.ord, type="samples", color="Antibiotic_treatment_1_month") +
  coord_fixed() +
  #ggtitle("antibiotics 1 month - only before") +
  geom_point(size = 2) +
  stat_ellipse(geom = "polygon", type="norm", alpha=0.1, aes(fill=Antibiotic_treatment_1_month)) +

  #geom_text(aes(label = Case_study), size = 4, vjust = 1) + 
  scale_color_manual(values = c("Donor" = "blue", "No" = "green3", "Yes" = "red"))+
  scale_fill_manual(values = c("Donor" = "blue", "No" = "green3", "Yes" = "red"))+
  xlab(paste0("PC1: ", percentVar[1], "% variance")) + 
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  labs(color = "Antibiotic treatment", fill = "Antibiotic treatment") +
  theme_bw()


metadata <- as(sample_data(psonly_before_for_DE), "data.frame")

adonis(phyloseq::distance(psonly_before_for_DE, method="bray") ~ Antibiotic_treatment_1_month, data = metadata)

#outcome taxonomy

plot_ordination(psonly_before_norm, PS.ord, type="samples", color="Outcome") +
  coord_fixed() +
  #ggtitle("Outcome - only before") +
  geom_point(size = 2) +
  stat_ellipse(geom = "polygon", type="norm", alpha=0.1, aes(fill=Outcome)) +
  scale_color_manual(values = c("Responders" = "green3", "Non-responders" = "red"),
                     labels = c("Failures", "Responders"))+
  scale_fill_manual(values = c("Responders" = "green3", "Non-responders" = "red"),
                    labels = c("Failures", "Responders"))+
  #geom_text_repel(aes(label = (ifelse(Case.study == "case8 ", "Case 8", NA))), 
   #               color = "black", size = 3,  show.legend = FALSE,
    #              nudge_y = -0.5, nudge_x = -0.3, segment.size = 0.2, hjust = 0.5, direction = "x")+
  xlab(paste0("PC1: ", percentVar[1], "%")) + 
  ylab(paste0("PC2: ", percentVar[2], "%")) +
  theme_bw()
 


metadata <- as(sample_data(psonly_before_for_DE), "data.frame")

adonis(phyloseq::distance(psonly_before_for_DE, method="jaccard") ~ Outcome, data = metadata)


#outcome functional

vsd.deseq <- varianceStabilizingTransformation(countdata.deseq, blind = TRUE)
pcaData <- plotPCA(vsd.deseq, intgroup = c( "Outcome", "Timepoint"), returnData = TRUE)
pcaData 

percentVar <- round(100 * attr(pcaData, "percentVar"),digits = 1)

ggplot(pcaData, aes(x = PC1, y = PC2, color = Outcome)) +
  coord_fixed()+
  #ggtitle("Pre-FMT with donors functional (Modules) \n Outcome") +
  #geom_text(label = subsetColData$Case_study, vjust = 1) +
  geom_point(size = 2) + #the size of the points on the graph 
  xlab(paste0("PC1: ", percentVar[1], "% variance")) + 
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  stat_ellipse(geom = "polygon", type="norm", alpha=0.1, aes(fill=Outcome)) +
  scale_color_manual(values = c("Responders" = "green3", "Non-responders" = "red"),
                     labels = c("Failures", "Responders"))+
  scale_fill_manual(values = c("Responders" = "green3", "Non-responders" = "red"),
                    labels = c("Failures", "Responders"))+
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
             col = "grey80", pch = 21) + 
  ggtitle("Modules, FC>2, alpha<0.05") +
  geom_point(data = subset(DEres, abs(log2FoldChange) > 2 & padj < 0.05),
             aes(log2FoldChange, -log10(pvalue)), color = "red", size =4)  
geom_text(data = subset(DEres, abs(log2FoldChange) > 2 & padj < 0.05),
          aes(log2FoldChange, -log10(pvalue)), label = rownames(subset(DEres, abs(log2FoldChange) > 2 & padj < 0.05)), vjust =1)

#normalizing counts
sizeFactors <- estimateSizeFactors(countdata.deseq)
normCounts.deseq <- counts(sizeFactors ,normalized=TRUE)

#subsetting the data to include only the genes with highest DE
highestDE <- subset(DEres, abs(log2FoldChange) > 4 & padj < 0.05)
nrow(highestDE)
highestDEnames <- rownames(highestDE)

#creating a vector to filter the donors out
subsetDE2ColData <- allColData_filtered_samples[ grep("Responders"|"Non-responders", allColData_filtered_samples$Outcome), ]
subsetDE2Samples <- as.vector(rownames(subsetDE2ColData))

#preparing the data for heatmap analysis, subsetting it by the most DE modules
highestDEexp <- as.data.frame(subset(normCounts.deseq, rownames(normCounts.deseq) %in% highestDEnames,
                                     colnames(normCounts.deseq) %in% subsetDE2Samples))
#num of clusters can be changed according to heatmap results
num_clusters <- 2

heatmap.numClusters <- pheatmap::pheatmap(highestDEexp, 
                                          scale = "row",
                                          show_rownames = T, 
                                          #main = "Pre-FMT with donors functional (modules) \n sig modules fc>4, alpha<0.05",
                                          border_color = "NA",
                                          #cutree_cols = num_clusters,
                                          #cutree_rows = num_clusters,
                                          clustering_method = "complete",
                                          color = colorRampPalette(rev(RColorBrewer::brewer.pal(10, "RdBu")))(256))

#extracting to file the most diff expressed modules
alpha = 0.05
sigtab = DEres[which(DEres$padj < alpha), ]
#order results by pvalues
sigtab = sigtab[order(sigtab$pvalue),]
head(sigtab)
write.csv(sigtab,"only before/significant modules pre-fmt.csv", row.names = TRUE)

#########################################
#     differential analysis - taxonomy  #
#########################################

res = results(DES_psonly_before, contrast=c("Outcome", "Responders", "Non-responders"))
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(psonly_before)[rownames(sigtab), ], "matrix"))

#order results by pvalues
sigtab = sigtab[order(sigtab$pvalue),]
head(sigtab)

write.csv(sigtab,"only before/significant OTUs only before responders vs non.csv", row.names = TRUE)
