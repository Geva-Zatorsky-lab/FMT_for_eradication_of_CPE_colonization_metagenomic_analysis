############
# barplots #
############

#tax_glom merges species that have the same taxonomy at a certain taxaonomic rank, you can change the rank to fit your needs
#later, normalizing the data

ps_bar <- tax_glom(psfiltered, "Family")
ps_samples_timepoint <- merge_samples(ps_bar, "group")

ps_bar0 <- transform_sample_counts(ps_samples_timepoint, function(x) x / sum(x))
#plot_bar(ps_bar0, fill="Family")

#melts the table to data frame
df <- psmelt(ps_bar0)

#finds the most abundant taxa by mean abundance and arranging them by order
top_Family <- df %>%
  group_by(Sample, Family, Timepoint) %>%
  summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean)
top_Family

positions <- c( "Donor","Responders, post-FMT","Responders, pre-FMT","Non-responders, post-FMT", "Non-responders, pre-FMT")
positions_names <- c( "Donors","Responders\npost-FMT","Responders\npre-FMT","Failures\npost-FMT", "Failures\npre-FMT")

#takes top 10 taxa from the original filtered data frame, and plots
top10 <- top_Family$Family[1:32]
df0 <- df %>%
  mutate(Family = fct_other(Family, top10))
ggplot(df0, aes(Sample, Abundance, fill = Family)) +  scale_fill_brewer(palette = "Paired") +
  geom_col() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) + 
  #ggtitle("barplots") +
  #facet_grid(col = vars(group), scales = "free") +
  ylab("Relative Abundance") +
  xlab(NULL)+
  scale_x_discrete(limits = positions, labels = positions_names) +
  theme(strip.background = element_blank(), strip.text.x = element_blank())
  
#saves last figure

save_plot(
  filename = "Graphs/Barplots responders vs non vs donors family.tif",
  fig = last_plot(),
  width = 12,
  height = 9,
  dpi = 600,
  theme = theme_get(),
  label.color = "black",
  label.size = 2.4,
  axis.textsize = 0.7,
  axis.titlesize = 0.75,
  legend.textsize = 0.6,
  legend.titlesize = 0.65,
  legend.itemsize = 0.5
)

