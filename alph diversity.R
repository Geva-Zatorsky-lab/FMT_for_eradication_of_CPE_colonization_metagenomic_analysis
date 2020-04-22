###############################
# Alpha diversity all samples #
###############################

#use filtered but not normalized data, to account for singeltons.

correct_group_names <- list(
        'Donor'="Donors",
        'Responders'="Responders",
        'Non-responders'="Failures"
)

name_labeller <- function(variable,value){
        return(correct_group_names[value])
}

positions <- c("Donor", "Responders", "Failures")


p = plot_richness(psfiltered_taxa, measures=c("Shannon"), color = "Timepoint", title = "Shannon Diversity Index")

p$data$Outcome_n = factor(p$data$Outcome, levels = c("Donor", "Responders", "Non-responders"))

p + geom_point(size=2, alpha=0.7) +
        aes(x = reorder(Case.number, sort(as.numeric(sample.number)))) +
        facet_grid(~ Outcome_n, scales = "free", labeller = name_labeller, switch="both") +
        xlab("Case") +
        scale_color_manual(values = c("Donor" = "blue", "Post-FMT" = "green3", "Pre-FMT" = "red"))+
        theme_bw() +
        theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 6))

save_plot(
        filename = "Graphs/Shannon alpha diversity.tif",
        fig = last_plot(),
        width = 12,
        height = 9,
        dpi = 600,
        theme = theme_get(),
        label.color = "black",
        label.size = 2.4,
        axis.textsize = 0.8,
        axis.titlesize = 0.75,
        legend.textsize = 0.6,
        legend.titlesize = 0.65,
        legend.itemsize = 0.5
)
