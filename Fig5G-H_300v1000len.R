#Figure 5G-H and Table S6. Comparison of Illumina sequencing with the Illumina 1000-cycle kit and 300-cycle kit platforms. 
#Fig5G: Wilcoxon rank-sum test was performed to compare the density of the length of the viral contigs output by the 1000-cycle kit and 300-cycle kit. 
#Fig5H: Wilcoxon rank-sum test was performed to compare the density of the length of all contigs output by the 1000-cycle kit and 300-cycle kit.


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
library(tidyr)
library(stringr)
library(patchwork)
library(tidyverse)

######## Get the data
# Get the metadata
setwd("/Users/jduan/bushman/virome_methods/250429_nextseq")
metadata <-read_excel("metadata.xlsx")

########################### Plot for viral contig lengths
# Import the cenote result
setwd("/Users/jduan/bushman/virome_methods/250429_nextseq/300bp_vs_1000bp")
cenote_df <- read_excel("300bpvs1000bp_cenote_R_analysis.xlsx", sheet="cenote_result")
cenote_df_300bp <- cenote_df %>% filter(cycle_number=="300")
cenote_df_1000bp <- cenote_df %>% filter(cycle_number=="1000")

# Create density plot for viral contig length
library(viridis)
library(hrbrthemes)

# get tsv file containing total number of bases per sample
setwd("/Users/jduan/bushman/virome_methods/250429_nextseq/300bp_vs_1000bp")
contig_length_df <- read.table("per_sample_bases.tsv", header=TRUE, sep="\t")
contig_length_df <- contig_length_df %>%
  separate(sample, into = c("Sample_ID", "cycle_number"), sep = "_", extra = "drop", remove = FALSE) # separate the sample name to make fields Sample_ID and cycle_number

# normalize the viral contig length
vcontig_len_cenote_df <- cenote_df %>% left_join(contig_length_df %>% select (Sample_ID, cycle_number, total_bases),
                                                   by = c("Sample_ID" = "Sample_ID", "cycle_number"="cycle_number")) #match the sum rpkm to their respective group&type so that each input reference in each group has a matching total rpkm
vcontig_len_cenote_df <- vcontig_len_cenote_df %>%
  mutate(virus_seq_length = ifelse(is.na(virus_seq_length), 0, virus_seq_length)) # replace the NA with 0 for plotting

vcontig_len_cenote_df <- vcontig_len_cenote_df %>%
  mutate(norm_vcontig_len=virus_seq_length/(total_bases/10^6)) #normalize based on total number of bp sequenced per sample

vcontig_len_cenote_df <- vcontig_len_cenote_df %>% filter(grepl("47|48|49|50", Sample_ID))

# make the density plot (normalized vcontig len)
vcontig_density_df <- vcontig_len_cenote_df %>% filter(Sample_type=="stool")

#################################### Statistical Tests for viral contigs
# Perform Wilcoxon rank sum test
library(rstatix)
cycle_comparison <- list(c("300", "1000")) #define the comparison first
vcontig_wilcox_result <- vcontig_density_df %>% group_by(Mock_Community, Sample_type) %>% 
  wilcox_test(norm_vcontig_len ~ cycle_number,
              comparisons = cycle_comparison) %>%
  adjust_pvalue(method="BH") %>% #optional multiple testing correction
  add_significance("p.adj")

# Plot with annotation
vcontig_density_compare_plot <- ggplot(vcontig_density_df, aes(x = cycle_number, y = norm_vcontig_len, fill = cycle_number)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.3) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.format", size=3, # shows formatted p-value
    comparisons = cycle_comparison,
    label.y = log10(max(vcontig_density_df$norm_vcontig_len) * 1.5) #push label a bit above the highest point
  ) +
  facet_nested(. ~ Mock_Community + Sample_type, 
               scales = 'free_x', 
               nest_line = element_line(linetype=1)) +
  scale_y_log10(labels=scales::label_number(accuracy=1)) +
  #ggtitle("Viral contig length distribution comparison (1000-cycles vs 300-cycles)") +
  labs(x="cycle number", y="normalized viral contig length\n(bp per million sample bases sequenced)") +
  theme_minimal() +
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE)

# Calculate the median
vcontig_median_df <- vcontig_len_cenote_df %>%
  group_by(Mock_Community, Sample_type, cycle_number) %>%
  summarise(median_len = median(virus_seq_length),
            median_normlen = median(norm_vcontig_len), .groups = "drop",
            vcontig_num=n())

######################### Do the same for all contigs assembled by megahit
setwd("/Users/jduan/bushman/virome_methods/250429_nextseq/300bp_vs_1000bp")
allcontig_300bp <-read.table("300bp_megahit_contig_table.tsv", header=TRUE, sep="\t") %>% 
  mutate(sample = sub("(_[^_]*){2}$", "", contig.id)) %>% #removes everything after second-to-last underline
  mutate(cycle_number="300")
allcontig_300bp <- allcontig_300bp %>%
  left_join(metadata %>% select(Index, Sample_ID, Mock_Community, Sample_type, Experiment, Treatment, Nuclease), by = c("sample"="Index"))
allcontig_300bp <- allcontig_300bp %>% left_join(contig_length_df %>% select (Sample_ID, cycle_number, total_bases),
                                                 by = c("Sample_ID" = "Sample_ID", "cycle_number"="cycle_number"))

allcontig_1000bp <-read.table("1000bp_megahit_contig_table.tsv", header=TRUE, sep="\t") %>% 
  mutate(sample = sub("_.*", "", contig.id)) %>% #removes everything after the first underline
  mutate(cycle_number="1000")
allcontig_1000bp <- allcontig_1000bp %>%
  left_join(metadata %>% select(Index, Sample_ID, Mock_Community, Sample_type, Experiment, Treatment, Nuclease), by = c("sample"="Sample_ID")) %>% 
  rename("Sample_ID"="Index")
allcontig_1000bp <- allcontig_1000bp %>% left_join(contig_length_df %>% select (Sample_ID, cycle_number, total_bases),
                                                   by = c("sample" = "Sample_ID", "cycle_number"="cycle_number"))

allcontig_df <- rbind(allcontig_300bp, allcontig_1000bp) #combine the two tables

# normalize the viral contig length
allcontig_normlen_df <- allcontig_df %>%
  mutate(norm_contig_len=virus_seq_length/(total_bases/10^6)) #normalize based on total number of bp sequenced per sample

# make the density plot (normalized contig len)
allcontig_density_df <- allcontig_normlen_df %>% filter(Sample_type=="stool")

#################################### Statistical Tests for all contigs
# Perform Wilcoxon rank sum test
allcontig_wilcox_result <- allcontig_density_df %>% group_by(Mock_Community, Sample_type) %>% 
  wilcox_test(norm_contig_len ~ cycle_number,
              comparisons = cycle_comparison) %>%
  adjust_pvalue(method="BH") %>% #optional multiple testing correction
  add_significance("p.adj")

# Plot with annotation
allcontig_density_compare_plot <- ggplot(allcontig_density_df, aes(x = cycle_number, y = norm_contig_len, fill = cycle_number)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.3, outlier.shape = NA, alpha = 0.3) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.format", size=3, # shows formatted p-value
    comparisons = cycle_comparison,
    label.y = log10(max(allcontig_density_df$norm_contig_len) * 1.5) #push label a bit above the highest point
  ) +
  facet_nested(. ~ Mock_Community + Sample_type, 
               scales = 'free_x', 
               nest_line = element_line(linetype=1)) +
  scale_y_log10(labels=scales::label_number(accuracy=1)) +
  #ggtitle("Contig length distribution comparison (1000-cycles vs 300-cycles)") +
  labs(x="cycle number", y="normalized contig length\n(bp per million sample bases sequenced)") +
  theme_minimal() +
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE)


# Calculate the median
allcontig_median_df <- allcontig_normlen_df %>%
  group_by(Mock_Community, Sample_type, cycle_number) %>%
  summarise(median_len = median(virus_seq_length),
            median_normlen = median(norm_contig_len), .groups = "drop",
            total_contig_num=n())



####################################### Save plots
# enable showtext
#library(showtext)
#font_add("Arial", regular = "/System/Library/Fonts/Supplemental/Arial.ttf")  # macOS path
#showtext_auto()  # enables showtext for all plots

# get the plots that are going to be saved
plots <- list(vcontig_density_compare_plot, allcontig_density_compare_plot)

# specify the destination folder
dest_folder <- "/Users/jduan/bushman/virome_methods/paper/figures_pdf"
if(!dir.exists(dest_folder)) dir.create(dest_folder, recursive = TRUE)

# loop over plots and save as PDF
for (i in seq_along(plots)) {
  print(plots[[i]])   # ensures the plot is drawn
  ggsave(
    filename = file.path(dest_folder, paste0("Fig5_", i, ".pdf")),
    plot = plots[[i]],
    width = 6.5,
    height = 7.5,
    units = "in"
    # no need to specify device; showtext handles font embedding
  )
}


################## Export table results as excel
# Create a blank workbook
OUT <- createWorkbook()

# Add some sheets to the workbook
addWorksheet(OUT, "vcontig_len_median")
addWorksheet(OUT, "allcontig_len_median")

# Write the data to the sheets
writeData(OUT, sheet = "vcontig_len_median", x = vcontig_median_df)
writeData(OUT, sheet = "allcontig_len_median", x = allcontig_median_df)

# Export the file
setwd("/Users/jduan/bushman/virome_methods/250429_nextseq/300bp_vs_1000bp")
saveWorkbook(OUT, "300bp_vs_1000bp_contiglen.xlsx")


