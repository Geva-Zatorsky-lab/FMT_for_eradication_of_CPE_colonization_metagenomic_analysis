#load needed libraries, assuming all of them are installed
library(phyloseq)
library(knitr)
library(readr)
library(microbiome)
library(pheatmap)
library(RColorBrewer)
library("metagMisc")
library(dplyr)
library(forcats)
library(ggplot2)
library(tidytext)
library(DESeq2)
library(vsn)
library(sjPlot)
library(stringr)
library(ggrepel)
library(vegan)
library(reshape2)
library(pairwiseAdonis)


#create a phyloseq object for taxonomy
otu_table <- read.csv(file = "Counts just bact archea_new cases.csv", row.names = 1)
otu_matrix <- read.csv(file = "taxa bact and archea cleaned.csv", row.names = 1)
otu_matrix <- as.matrix(otu_matrix)

metadata <- read.csv(file = "sample_metadata_short.csv", row.names = 1)


OTU = otu_table(otu_table, taxa_are_rows = TRUE)
TAX = tax_table(otu_matrix)
meta = sample_data(metadata)

ps = phyloseq(OTU, TAX, meta)

#filtering samples and normalizing (TSS)

psfiltered = subset_samples(ps, 
                            sample_names(ps)  != "case9_after" & sample_names(ps)  != "case1_before" 
                            & sample_names(ps)  != "donor_2" & sample_names(ps)  != "donor_3")

sample_names(psfiltered)



#the filter used is minimum 100 reads and 20% prevalance in samples
filter <- phyloseq::genefilter_sample(psfiltered, filterfun_sample(function(x) x >= 100), 
                                      A = 0.2*nsamples(psfiltered))
psfiltered_taxa = prune_taxa(filter, psfiltered)



psfiltered_norm = phyloseq_standardize_otu_abundance(psfiltered_taxa, method = "total")


#open files for functional analysis

allCountData <- read.csv(file ="funl_modules.csv", row.names = 1)

allColData <- read.csv(file ="sample_metadata_short.csv", row.names = 1)

allColData_filtered_samples <- subset(allColData, row.names(allColData) != ("case1_before"))
allColData_filtered_samples <- subset(allColData_filtered_samples, row.names(allColData_filtered_samples) != ("donor_2"))
allColData_filtered_samples <- subset(allColData_filtered_samples, row.names(allColData_filtered_samples) != ("case9_after"))
allColData_filtered_samples <- subset(allColData_filtered_samples, row.names(allColData_filtered_samples) != ("donor_3"))

