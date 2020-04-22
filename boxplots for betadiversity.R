dist <- phyloseq::distance(psonly_after_nor_for_DE, method="jaccard")

metadata <- as(sample_data(psonly_after_nor_for_DE), "data.frame")

vect <- metadata$group

bd<-betadisper(dist, vect )
boxplot(bd)


df <- data.frame(Distance_to_centroid=bd$distances,Group=bd$group)
groups <- bd$group

p<- ggplot(data=df,aes(x=Group,y=Distance_to_centroid,colour=groups))
p <-p + geom_point(alpha=0.5) + geom_text_repel(aes(label = row.names(df)))
p <-p + theme(axis.text.x=element_text(angle = -90, hjust = 0,vjust=0.5))
p

psonly_after_nor_for_DE = subset_samples(psonly_after_norm, (Timepoint == "Post-FMT"))

adonis(phyloseq::distance(psonly_after_nor_for_DE, method="jaccard") ~ Outcome, data = metadata)
