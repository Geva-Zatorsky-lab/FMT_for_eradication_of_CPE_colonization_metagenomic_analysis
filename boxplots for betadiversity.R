dist <- phyloseq::distance(psfiltered_norm, method="bray")

metadata <- as(sample_data(psfiltered_norm), "data.frame")

psonly_after_nor_for_DE = subset_samples(psonly_after_norm, (Timepoint == "Post-FMT"))

vect <- metadata$group

bd<-betadisper(dist, vect )
boxplot(bd)

df <- data.frame(Distance_to_centroid=bd$distances,Group=bd$group)
groups <- bd$group

grouping_colors <- c("Donor" = "blue", 
                     "Responders, pre-FMT" = "darkorange", 
                     "Non-responders, pre-FMT" = "purple",
                     "Responders, post-FMT" = "green3", 
                     "Non-responders, post-FMT" = "red")

group_labels = c("Donors","Failures, post-FMT", "Failures, pre-FMT", "Responders, post-FMT", "Responders, pre-FMT")

positions <- c( "Donor","Responders, post-FMT","Responders, pre-FMT","Non-responders, post-FMT", "Non-responders, pre-FMT")
positions_names <- c( "Donors","Responders\npost-FMT","Responders\npre-FMT","Failures\npost-FMT", "Failures\npre-FMT")

ggplot(data=df,aes(x=Group,y=Distance_to_centroid,colour=groups)) +
        #geom_point(alpha=0.5) + 
        geom_boxplot(alpha=0.5, aes(fill = groups)) +
        scale_color_manual(values = grouping_colors,labels = group_labels) +
        scale_fill_manual(values = grouping_colors,labels = group_labels) +   
        #geom_text_repel(aes(label = row.names(df))) +
        theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) + 
        scale_x_discrete(limits = positions, labels = positions_names) +
        ylab("Distance to centroid") +
        xlab(NULL) +
        theme(strip.background = element_blank(), strip.text.x = element_blank(), legend.position = "none")



permutest(bd, permutations = 999)
anova(bd)
plot(bd)
boxplot(bd)
plot(TukeyHSD(bd))
