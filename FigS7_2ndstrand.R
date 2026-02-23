# This script combines the reference-based and contig-based results of 260116_nextseq and plots using facet by conditions
library(data.table)
library(readxl)
library(ggpubr)
library(openxlsx)
library(dplyr)
library(patchwork)
library(plotly)
library(ggplot2)
library(ggh4x)
library(paletteer)

# Import the data
setwd("/Users/jduan/bushman/virome_methods/Second-strand synthesis/260117_nextseq/R analysis")
ref_df <-read_excel("260116_nextseq_R_analysis.xlsx", sheet="ref_viral_analysis")
ref_reads_df <-read_excel("260116_nextseq_R_analysis.xlsx", sheet="ref_viral_percentage")
cenote_df <-read_excel("260116_nextseq_cenote_R_analysis.xlsx", sheet="cenote_result")
cenote_reads_df <-read_excel("260116_nextseq_cenote_R_analysis.xlsx", sheet="viral_percentage")

# Combine the two viral percentages tables into one
cenote_reads_df <- cenote_reads_df %>% mutate(contig_type="viral")
ref_reads_df <- ref_reads_df %>% mutate(contig_type="ref viral") %>% rename(percent_viral_reads=percent_ref_reads, sum_viral_reads=sum_viralref_reads)
viral_reads_abundance <- rbind(cenote_reads_df, ref_reads_df)


# Create a bar plot for the ref_rpkm of each replicate (+ ref rpkm annotated on graph)
i_ref_abundance_plot <- ref_df %>% select(Experiment, Sample_ID, Sample_type, Mock_Community, Treatment, RT, thermocycling_condition, Klenow, purification, reference_genome, ref_rpkm, mapped_ref_reads, total_reads)
i_ref_abundance_plot <- i_ref_abundance_plot %>% left_join(ref_reads_df %>% ungroup() %>% select (Experiment, Sample_ID, sum_replicate_rpkm, percent_viral_reads),
                                                           by = c("Sample_ID" = "Sample_ID", "Experiment"="Experiment")) #match the sum rpkm to their respective Sample_ID so that each input reference in each group has a matching total rpkm
i_ref_abundance_plot <- i_ref_abundance_plot %>% mutate(ref_relative_abundance=ref_rpkm/sum_replicate_rpkm) #calculate the relative abundance of each input reference

############### Stacked bar graph for reference-based result
Exp1_ref_abundance <- i_ref_abundance_plot %>% filter(Mock_Community=="mock community spiked-in")

p_Exp1_ref <- ggplot(Exp1_ref_abundance, aes(x=Sample_ID, y=ref_relative_abundance, fill=reference_genome)) +
  geom_bar(stat="identity") + 
  scale_fill_paletteer_d("MoMAColors::Lupi") +
  scale_x_discrete(drop = T) +
  facet_nested(. ~ RT + thermocycling_condition + Klenow + purification, 
               scales = 'free_x', 
               nest_line = element_line(linetype=1)) +
  labs(x="Sample", y="relative abundance based on RPKM", 
       #title = "relative abundance based on reference genomes (2nd-strand synthesis experiment)"
  ) +
  theme_minimal(base_family="Helvetica") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background=element_blank(),
        ggh4x.facet.nestline=element_line(color="black"))

############### Experiment 1 percent viral reads result
# Make the total reads numeric
viral_reads_abundance$total_reads <- as.numeric(viral_reads_abundance$total_reads)

Exp1_reads_percentage <- viral_reads_abundance %>% filter(Mock_Community=="mock community spiked-in")
Exp1_reads_percentage <- Exp1_reads_percentage %>% mutate(condition=paste(RT, "|", thermocycling_condition, "|", Klenow, "|", purification))
q_Exp1 <- ggplot(Exp1_reads_percentage, aes(x = condition, y = percent_viral_reads)) +
  geom_jitter(aes(fill = contig_type, color = contig_type, size=total_reads),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
              #size=3, 
              alpha = 0.7, shape = 21) +
  stat_summary(fun = mean, 
               geom = "crossbar",
               aes(group = contig_type),
               position = position_dodge(width = 0.75), width = 0.25, fatten = 1.5, color = "black", show.legend = FALSE) +
  scale_x_discrete(drop = TRUE) +
  facet_nested(. ~ RT + thermocycling_condition + Klenow + purification,
               scales = "free_x", nest_line = element_line(linetype = 1)) +
  scale_y_continuous(limits = c(0, 15)) +
  scale_size_continuous(name = "Total reads", labels = scales::scientific, range = c(1, 6)) +
  # Optional: Better color palettes
  scale_fill_brewer(palette = "Set2") + scale_color_brewer(palette = "Set2") +
  guides(size = guide_legend(override.aes = list(shape = 21, fill = "black",color = "black"))) +
  labs(
    x = "Sample", y = "Percent in total reads", 
    #title = "Viral reads percentage (2nd-strand synthesis experiment)",
    fill = "Condition", color = "Condition") +
  theme_bw(base_family="Helvetica") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_blank(),
    ggh4x.facet.nestline = element_line(color = "black")
  )


# get the plots that are going to be saved
plots <- list(p_Exp1_ref, q_Exp1)

# specify the destination folder
dest_folder <- "/Users/jduan/bushman/virome_methods/paper/manuscript/msystems_submission2/figure"
if(!dir.exists(dest_folder)) dir.create(dest_folder, recursive = TRUE)

# loop over plots and save as PDF
for (i in seq_along(plots)) {
  print(plots[[i]])   # ensures the plot is drawn
  ggsave(
    filename = file.path(dest_folder, paste0("Resub_SuppFig7_", i, ".pdf")),
    plot = plots[[i]],
    width = 6.5,
    height = 7.5,
    units = "in"
    # no need to specify device; showtext handles font embedding
  )
}

