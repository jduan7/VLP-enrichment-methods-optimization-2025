# Figure S2. Comparison of viral sequence recovery after virus particle enrichment followed by sequencing, 
# versus direct metagenomic sequencing of total DNA or RNA.(geNomad-based viral annotation) 
library(data.table)
library(readxl)
library(ggpubr)
library(openxlsx)
library(dplyr)
library(patchwork)
library(plotly)
library(stringr)
library(ggh4x)
library(paletteer)
library(purrr)

# Get the metadata
setwd("/Users/jduan/bushman/virome_methods/metagenomic sequencing/neat_stoolBALsaliva_genomad")
metadata <- read_excel("stoolBALsaliva_comparison_metadata.xlsx")

# Import the data
stool_VP1VP2_cov <- read.table("stool_VP1VP2_allcontig_genomad_cov.tsv", header=TRUE, sep="\t") %>%
  mutate(seq_type="virome")
saliva_BAL_stoolVP3_cov <- read.table("saliva_BAL_stoolVP3_allcontig_genomad_cov.tsv", header=TRUE, sep="\t") %>%
  mutate(seq_type="virome")
metaDNA_stoolBALsaliva_cov <- read.table("stoolBALsaliva_neat_DNA_genomad_cov.tsv", header=TRUE, sep="\t") %>%
  mutate(seq_type="metagenomic DNAseq")
metaRNA_stoolBALsaliva_cov <- read.table("stoolBALsaliva_neat_RNA_genomad_cov.tsv", header=TRUE, sep="\t") %>%
  mutate(seq_type="metagenomic RNAseq")

# Combine the tables
cov_df <- rbind(stool_VP1VP2_cov, saliva_BAL_stoolVP3_cov, metaDNA_stoolBALsaliva_cov, metaRNA_stoolBALsaliva_cov) %>% rename(sample=raw_read)
cov_df <- cov_df %>%
  left_join(metadata %>% select(Experiment, Sample_ID, Mock_Community, Sample_type, Treatment, Sequencing), by = c("sample"="Sample_ID", "seq_type"="Sequencing")) %>% 
  rename(num_total_reads=total, num_mapped_reads=mapped)

# calculate the viral percentage
viral_percentage <- cov_df %>% group_by(Experiment, sample, Mock_Community, Sample_type, Treatment, seq_type) %>% 
  summarize(sum_viral_reads = sum(num_mapped_reads), total_reads=first(num_total_reads),
            percent_viral_reads = sum(num_mapped_reads)/first(num_total_reads)*100)

viral_percentage_plot <- viral_percentage %>% 
  mutate(condition=paste(Sample_type, Treatment, sep=" | "))

# Plot for stool viral percentage comparison
stool_viral_reads_abundance <- viral_percentage_plot %>% filter(Sample_type=="stool")

q_stool <- ggplot(stool_viral_reads_abundance, aes(x=condition, y=percent_viral_reads, fill=seq_type)) +
  geom_jitter(aes(fill = seq_type, color = seq_type),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
              alpha = 0.7, shape = 21) +
  stat_summary(fun = mean, 
               geom = "crossbar",
               aes(group = seq_type),
               position = position_dodge(width = 0.75), width = 0.25, fatten = 1.5, color = "black", show.legend = FALSE) +
  scale_x_discrete(drop = TRUE) +
  facet_nested(. ~ Mock_Community + Sample_type + Treatment,
               scales = "free_x", nest_line = element_line(linetype = 1)) +
  scale_y_continuous(limits = c(0, 100)) +
  #scale_size_continuous(name = "Total reads", labels = scales::scientific, range = c(1, 6)) +
  # Optional: Better color palettes
  scale_fill_brewer(palette = "Set1") + scale_color_brewer(palette = "Set1") +
  guides(size = guide_legend(override.aes = list(shape = 21, fill = "black",color = "black"))) +
  labs(
    x = "Sample", y = "Percent in total reads", 
    #title = "Viral reads percentage (viromeprep vs metagenomic sequencing)",
    fill = "Analytical Method", color = "Analytical Method") +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_blank(),
    ggh4x.facet.nestline = element_line(color = "black")
  )

SalivaBAL_viral_reads_abundance <- viral_percentage_plot %>% filter(Sample_type %in% c("saliva","BAL"))
q_SalivaBAL <- ggplot(SalivaBAL_viral_reads_abundance, aes(x=condition, y=percent_viral_reads, fill=seq_type)) +
  geom_jitter(aes(fill = seq_type, color = seq_type),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
              alpha = 0.7, shape = 21) +
  stat_summary(fun = mean, 
               geom = "crossbar",
               aes(group = seq_type),
               position = position_dodge(width = 0.75), width = 0.25, fatten = 1.5, color = "black", show.legend = FALSE) +
  scale_x_discrete(drop = TRUE) +
  facet_nested(. ~ Mock_Community + Sample_type + Treatment,
               scales = "free_x", nest_line = element_line(linetype = 1)) +
  scale_y_continuous(limits = c(0, 100)) +
  #scale_size_continuous(name = "Total reads", labels = scales::scientific, range = c(1, 6)) +
  # Optional: Better color palettes
  scale_fill_brewer(palette = "Set1") + scale_color_brewer(palette = "Set1") +
  guides(size = guide_legend(override.aes = list(shape = 21, fill = "black",color = "black"))) +
  labs(
    x = "Sample", y = "Percent in total reads", 
    #title = "Viral reads percentage (viromeprep vs metagenomic sequencing)",
    fill = "Analytical Method", color = "Analytical Method") +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_blank(),
    ggh4x.facet.nestline = element_line(color = "black")
  )


####################################### Save plots
# enable showtext
#library(showtext)
#font_add("Arial", regular = "/System/Library/Fonts/Supplemental/Arial.ttf")  # macOS path
#showtext_auto()  # enables showtext for all plots

# get the plots that are going to be saved
plots <- list(q_stool, q_SalivaBAL)

# specify the destination folder
dest_folder <- "/Users/jduan/bushman/virome_methods/paper/manuscript/msystems_submission2/figure"
if(!dir.exists(dest_folder)) dir.create(dest_folder, recursive = TRUE)

# loop over plots and save as PDF
for (i in seq_along(plots)) {
  print(plots[[i]])   # ensures the plot is drawn
  ggsave(
    filename = file.path(dest_folder, paste0("Resub_FigS2_", i, ".pdf")),
    plot = plots[[i]],
    width = 6.5,
    height = 7.5,
    units = "in"
    # no need to specify device; showtext handles font embedding
  )
}
