############
# barplots #
############

#tax_glom merges species that have the same taxonomy at a certain taxaonomic rank, you can change the rank to fit your needs
#later, normalizing the data

ps_bar <- tax_glom(psfiltered_taxa, "Family")
ps_samples_timepoint <- merge_samples(ps_bar, "group")


ps_bar0 <- transform_sample_counts(ps_bar, function(x) x / sum(x))
#plot_bar(ps_bar0, fill="Species")

#melts the table to data frame
df <- psmelt(ps_bar0)

#finds the most abundant taxa by mean abundance and arranging them by order
top_Family <- df %>%
  group_by(Sample, Family, Timepoint) %>%
  summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean)
top_Family


positions <- c( "Donor","Responders, post-FMT","Responders, pre-FMT","Non-responders, post-FMT", "Non-responders, pre-FMT")
positions_names <- c( "Donor","Responders\npost-FMT","Responders\npre-FMT","Failuers\npost-FMT", "Failures\npre-FMT")


#takes top 10 taxa from the original filtered data frame, and plots
top10 <- top_Family$Family[1:62]
df0 <- df %>%
  mutate(Family = fct_other(Family, top10))

correct_group_names <- list(
  'A'="Donors",
  'B'="Responders\nPost-FMT",
  'C'="Responders\nPre-FMT",
  'D'="Failures\nPost-FMT",
  'E'="Failures\nPre-FMT"
)

name_labeller <- function(variable,value){
  return(correct_group_names[value])
}



ggplot(df0, aes(x = reorder(Case.number, sort(as.numeric(Case.number))), Abundance, fill = Family)) +  
  scale_fill_brewer(palette = "Paired") +
  geom_col() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 7)) + 
  #ggtitle("barplots") +
  
  facet_grid(.~ordering, scales = "free",labeller = name_labeller, switch="both") +
  #aes(x = reorder(Case.number, sort(metadata$Case.number))) +
  xlab(NULL) +
  ylab("Relative abundance")
  #scale_x_discrete(limits = sample_names(ps_bar0), labels = meta_barplot$Case.number) +
  theme(strip.background = element_blank(), strip.text.x = element_blank())
  
#saves last figure

save_plot(
  filename = "Graphs/Barplots responders vs non vs donors order grouped.tif",
  fig = last_plot(),
  width = 20,
  height = 9,
  dpi = 600,
  theme = theme_get(),
  label.color = "black",
  label.size = 0.8,
  axis.textsize = 0.8,
  axis.titlesize = 0.8,
  legend.textsize = 0.8,
  legend.titlesize = 1,
  legend.itemsize = 0.4
)

