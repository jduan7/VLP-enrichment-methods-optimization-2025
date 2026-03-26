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

# Get the metadata
setwd("/Users/jduan/bushman/virome_methods/250429_nextseq")
metadata_250505 <- read_excel("metadata.xlsx")
setwd("/Users/jduan/bushman/virome_methods/250616_nextseq")
metadata_250617 <- read_excel("metadata.xlsx")
metadata <- bind_rows(metadata_250505,metadata_250617)

# Get the output excel file for 250505_nextseq
setwd("/Users/jduan/bushman/virome_methods/250429_nextseq/R analysis")
ref_df_250505 <-read_excel("250505_nextseq_R_analysis.xlsx", sheet="ref_viral_analysis")
ref_reads_df_250505 <-read_excel("250505_nextseq_R_analysis.xlsx", sheet="ref_viral_percentage")
cenote_df_250505 <-read_excel("250505_nextseq_cenote_R_analysis.xlsx", sheet="cenote_result")
cenote_reads_df_250505 <-read_excel("250505_nextseq_cenote_R_analysis.xlsx", sheet="viral_percentage")

# Add columns for 250505_nextseq
ref_df_250505 <- ref_df_250505 %>% mutate(Amplification="none", VC_aliquot_used="250415")
ref_reads_df_250505 <- ref_reads_df_250505 %>% mutate(Amplification="none", VC_aliquot_used="250415")
cenote_df_250505 <- cenote_df_250505 %>% mutate(Amplification="none", VC_aliquot_used="250415")
cenote_reads_df_250505 <- cenote_reads_df_250505 %>% mutate(Amplification="none", VC_aliquot_used="250415")

# Get the output excel file for 250617_nextseq
setwd("/Users/jduan/bushman/virome_methods/250616_nextseq/R analysis")
ref_df_250617 <-read_excel("250617_nextseq_R_analysis.xlsx", sheet="ref_viral_analysis")
ref_reads_df_250617 <-read_excel("250617_nextseq_R_analysis.xlsx", sheet="ref_viral_percentage")
cenote_df_250617 <-read_excel("250617_nextseq_cenote_R_analysis.xlsx", sheet="cenote_result")
cenote_reads_df_250617 <-read_excel("250617_nextseq_cenote_R_analysis.xlsx", sheet="viral_percentage")

# Combine the two tables
ref_df <- rbind(ref_df_250505, ref_df_250617)
ref_reads_df <- rbind(ref_reads_df_250505, ref_reads_df_250617)
cenote_df <- rbind(cenote_df_250505, cenote_df_250617)
cenote_reads_df <- rbind(cenote_reads_df_250505, cenote_reads_df_250617)

# Combine the two viral percentages tables into one
cenote_reads_df <- cenote_reads_df %>% mutate(contig_type="viral")
ref_reads_df <- ref_reads_df %>% mutate(contig_type="ref viral") %>% rename(percent_viral_reads=percent_ref_reads, sum_viral_reads=sum_viralref_reads)
viral_reads_abundance <- rbind(cenote_reads_df, ref_reads_df)


# Plot the data
# Create a bar plot for the ref_rpkm of each replicate (+ ref rpkm annotated on graph)
i_ref_abundance_plot <- ref_df %>% 
  select(Experiment, Sample_ID, Sample_type, Mock_Community, Treatment, Nuclease, Amplification, VC_aliquot_used, reference_genome, ref_rpkm, mapped_ref_reads, total_reads) %>% 
  left_join(ref_reads_df %>% 
              ungroup() %>% 
              select (Experiment, Sample_ID, sum_replicate_rpkm, percent_viral_reads),
            by = c("Sample_ID" = "Sample_ID", "Experiment"="Experiment")) #match the sum rpkm to their respective Sample_ID so that each input reference in each group has a matching total rpkm
i_ref_abundance_plot <- i_ref_abundance_plot %>% mutate(ref_relative_abundance=ref_rpkm/sum_replicate_rpkm) #calculate the relative abundance of each input reference

i_ref_abundance_plot <- i_ref_abundance_plot %>%
  mutate(ref_relative_abundance = ifelse(is.na(ref_relative_abundance), 0, ref_relative_abundance)) # Replace NAs with 0 to show empty bars


############### Experiment 1 percent viral reads result
# Make the total reads numeric
viral_reads_abundance$total_reads <- as.numeric(viral_reads_abundance$total_reads)

# Add the condition for plotting purpose
Exp1_reads_percentage <- viral_reads_abundance %>% filter(Experiment %in% c("saliva/OP/BAL spike-in")) %>% 
  filter(Mock_Community=="mock community spiked-in") %>% 
  mutate(condition=paste(Sample_type, Treatment, sep=" | "))

q_Exp1 <- ggplot(Exp1_reads_percentage, aes(x = condition, y = percent_viral_reads)) +
  geom_jitter(aes(fill = contig_type, color = contig_type),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
              size=3, alpha = 0.7, shape = 21) +
  stat_summary(fun = mean, 
               geom = "crossbar",
               aes(group = contig_type),
               position = position_dodge(width = 0.75), width = 0.25, fatten = 1.5, color = "black", show.legend = FALSE) +
  scale_x_discrete(drop = TRUE) +
  facet_nested(. ~ Mock_Community + Sample_type + Treatment,
               scales = "free_x", nest_line = element_line(linetype = 1)) +
  scale_y_continuous(limits = c(0, 50)) +
  # Optional: Better color palettes
  scale_fill_brewer(palette = "Set2") + scale_color_brewer(palette = "Set2") +
  labs(
    x = "Sample", y = "Percent in total reads", 
    fill = "Contig type", color = "Contig type") +
  theme_bw(base_family="Helvetica") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_blank(),
    ggh4x.facet.nestline = element_line(color = "black")
  )

# Calculate for average percent viral reads and stdev for saliva/OP wash/BAL samples (w/spike-in)
Exp1_reads_percentage %>% filter(Mock_Community=="mock community spiked-in") %>% group_by(Treatment, contig_type) %>%
  summarize(mean_perc_viral=mean(percent_viral_reads),
            sd_perc_viral=sd(percent_viral_reads))

####################################### Save plots
# enable showtext
#library(showtext)
#font_add("Arial", regular = "/System/Library/Fonts/Supplemental/Arial.ttf")  # macOS path
#showtext_auto()  # enables showtext for all plots

# get the plots that are going to be saved
plots <- list(q_Exp1)

# specify the destination folder
dest_folder <- "/Users/jduan/bushman/virome_methods/paper/manuscript/msystems_submission3/figures"
if(!dir.exists(dest_folder)) dir.create(dest_folder, recursive = TRUE)

# loop over plots and save as PDF
for (i in seq_along(plots)) {
  print(plots[[i]])   # ensures the plot is drawn
  ggsave(
    filename = file.path(dest_folder, paste0("FigS2_", i, ".pdf")),
    plot = plots[[i]],
    width = 6.5,
    height = 7.5,
    units = "in"
    # no need to specify device; showtext handles font embedding
  )
}
