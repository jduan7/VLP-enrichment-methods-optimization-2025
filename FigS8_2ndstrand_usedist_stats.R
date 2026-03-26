# This script calculates Bray-Curtis distance among samples, then calculates tests whether Klenow vs no Klenow, 
# then groups Pairwise Bray–Curtis distances by Klenow treatment, extracts within-treatment distances, 
# and uses linear regression to test whether mean within-group dissimilarity differed between Klenow and no-Klenow samples.

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
library(tidyr)
library(vegan)

# Import the data
setwd("/Users/jduan/bushman/virome_methods/Second-strand synthesis/260117_nextseq")
metadata <-read_excel("metadata.xlsx")
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

# Construct a sample-by-species table (row=sample, column=species) with relative abundance
ref_abun_mat <- i_ref_abundance_plot %>% filter(Mock_Community=="mock community spiked-in") %>% 
  select(Sample_ID, reference_genome, ref_relative_abundance) %>% 
  pivot_wider(names_from=reference_genome,values_from=ref_relative_abundance)
ref_abun_mat <- data.frame(ref_abun_mat) #convert to dataframe
rownames(ref_abun_mat) <- ref_abun_mat$Sample_ID #make sample names as row names
ref_abun_mat <- ref_abun_mat %>% select(-Sample_ID) #get rid of the sample column
ref_abun_mat[is.na(ref_abun_mat)] <- 0 #set the relative abundance of the votu undetected in the sample as 0
ref_abun_mat_clean <- ref_abun_mat[rowSums(ref_abun_mat) > 0, ] #remove rows with all zeros for bray-curtis

# Compute Bray-Curtis distance to represent differences in abundance patterns
ref_bray <- vegdist(ref_abun_mat_clean, method="bray")

############# Linear fit model comparing within-group variance using Bray-Curtis distances
library(usedist)

# Make a sample table regarding Klenow treatment
Klenow_metadata <- i_ref_abundance_plot %>% filter(Mock_Community=="mock community spiked-in") %>% 
  group_by(Sample_ID) %>% summarize(Klenow=first(Klenow))

# Sort the distance matrix in the same order as the sample
ref_bray_sorted <- dist_subset(ref_bray, Klenow_metadata$Sample_ID)

# Create a data frame of distances between groups of items (group=Klenow or no Klenow; items=samples)
ref_bray_distance_groups <- dist_groups(ref_bray_sorted, Klenow_metadata$Klenow)

# Use linear fit to test: Does Klenow treatment change the average within-group Bray–Curtis dissimilarity?
lm_output <- ref_bray_distance_groups %>% filter(Group1==Group2) |> #filter for within-group distances
  lm(Distance~Group1, data=_) %>% summary() #modeling the distance as a function of condition (Group1 and 2 are the same so just do Group1)

# Visualize the result
lm_df <- ref_bray_distance_groups %>%
  filter(Group1 == Group2) %>%  #within-group only
  mutate(Klenow = Group1)

lm_plot <- ggplot(lm_df, aes(x = Klenow, y = Distance)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(aes(color = Klenow),
              width = 0.15, alpha = 0.3, size = 1.5) +
  stat_compare_means(
    method = "t.test",
    label = "p.format",
    comparisons = list(c("no Klenow", "with Klenow")),
    label.y = max(lm_df$Distance) * 1.05
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
  labs(
    x = "Klenow treatment",
    y = "Within-group Bray–Curtis distance"
  ) +
  theme_bw()

# get the plots that are going to be saved
plots <- list(lm_plot)

# specify the destination folder
dest_folder <- "/Users/jduan/bushman/virome_methods/paper/manuscript/msystems_submission2/figure"
if(!dir.exists(dest_folder)) dir.create(dest_folder, recursive = TRUE)

# loop over plots and save as PDF
for (i in seq_along(plots)) {
  print(plots[[i]])   # ensures the plot is drawn
  ggsave(
    filename = file.path(dest_folder, paste0("Resub_SuppFigS7_O", i, ".pdf")),
    plot = plots[[i]],
    width = 6.5,
    height = 7.5,
    units = "in"
    # no need to specify device; showtext handles font embedding
  )
}
