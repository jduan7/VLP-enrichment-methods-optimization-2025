#Fig5A-B. Comparison of Illumina sequencing with the Illumina 1000-cycle kit and 300-cycle kit platforms.
#A: Relative abundance of each reference virus in VirMock1 spiked-in stool that underwent VP1 and was then sequenced using the 1000-cycle kit or the 300-cycle kit 
#B: Relative abundance of classes of viruses as annotated by Cenote-Taker2 in stool with or without VirMock1 spiked-in that underwent VP1 and was then sequenced using the 1000-cycle kit or the 300-cycle kit.

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

####################

# Get the metadata
setwd("/Users/jduan/bushman/virome_methods/paper/codes/codes for publication/")
metadata <-read_excel("metadata_250505.xlsx")

#Add the reference results
setwd("/Users/jduan/bushman/virome_methods/paper/codes/codes for publication/")
ref_df <-read_excel("300bp_vs_1000bp_R_analysis.xlsx", sheet="ref_viral_analysis") %>% rename("cycle_number"="read_len")
ref_reads_df <-read_excel("300bp_vs_1000bp_R_analysis.xlsx", sheet="ref_viral_percentage") %>% rename("cycle_number"="read_len")

cenote_df <- read_excel("300bpvs1000bp_cenote_R_analysis.xlsx", sheet="cenote_result")
cenote_df_300bp <- cenote_df %>% filter(cycle_number=="300")
cenote_df_1000bp <- cenote_df %>% filter(cycle_number=="1000")

cenote_reads_df <-read_excel("300bpvs1000bp_cenote_R_analysis.xlsx", sheet="viral_percentage")
cenote_reads_df_300bp <- cenote_df %>% filter(cycle_number=="300")
cenote_reads_df_1000bp <- cenote_df %>% filter(cycle_number=="1000")

# Combine the ref viral reads and total viral reads
ref_reads_df <- ref_reads_df %>% mutate(contig_type="ref_viral") %>% rename(percent_viral_reads=percent_ref_reads, sum_viral_reads=sum_viralref_reads)
cenote_reads_df <- cenote_reads_df %>% mutate(contig_type="viral")
viral_reads <- rbind(ref_reads_df, cenote_reads_df)


####### Make relative abundance plot for all viral contigs
# Manually add back the samples that are missing due to the lack of assembled viral contigs from cenote result
i_cenote_df_300bp <- metadata %>% left_join(cenote_df_300bp, by = c("Sample_ID"="Sample_ID", "Experiment"="Experiment", "Mock_Community"="Mock_Community", "Sample_type"="Sample_type", "Treatment"="Treatment", "Nuclease"="Nuclease")) %>% 
  mutate(cycle_number = ifelse(is.na(cycle_number), "300", cycle_number))
i_reads_cenote_df_1000bp <- metadata %>% left_join(cenote_reads_df_1000bp, by = c("Sample_ID"="Sample_ID", "Experiment"="Experiment", "Mock_Community"="Mock_Community", "Sample_type"="Sample_type", "Treatment"="Treatment", "Nuclease"="Nuclease")) %>% 
  mutate(cycle_number = ifelse(is.na(cycle_number), "1000", cycle_number))
i_cenote_df <- rbind(i_cenote_df_300bp, i_reads_cenote_df_1000bp)
cenote_reads_df <- cenote_reads_df %>%
  mutate(percent_viral_reads = ifelse(is.na(percent_viral_reads), 0, percent_viral_reads)) # replace the NA with 0 for plotting

# filter for only the samples we want to compare
reads_cenote_Class_abundance <- i_cenote_df %>% filter(grepl("47|48|49|50", Sample_ID)) %>%
  select(Sample_ID, cycle_number, Mock_Community, Sample_type, Experiment, Treatment, organism, vcontig_rpkm, num_mapped_reads, num_total_reads, Class) #make a table with all the reads data and the Class output by cenote
reads_cenote_Class_abundance <- reads_cenote_Class_abundance %>% left_join(cenote_reads_df %>% ungroup() %>% select (Experiment, Sample_ID, cycle_number, sum_replicate_rpkm, percent_viral_reads),
                                                                           by = c("Sample_ID" = "Sample_ID", "Experiment" = "Experiment", "cycle_number"="cycle_number")) #match the sum rpkm to their respective group&type so that each input reference in each group has a matching total rpkm
reads_cenote_Class_abundance <- reads_cenote_Class_abundance %>% mutate(vcontig_relative_abundance=vcontig_rpkm/sum_replicate_rpkm) %>% #calculate the relative abundance of each input reference
  mutate(vcontig_relative_abundance = ifelse(is.na(vcontig_relative_abundance), 0, vcontig_relative_abundance)) # fill in the all-NA rows with 0 for later plotting

# based on grouping all the contigs of the same class together
reads_cenote_Class_Sum_abundance <- reads_cenote_Class_abundance %>% 
  group_by(Sample_ID, Class, Mock_Community, Treatment, Sample_type, cycle_number) %>% summarize(class_sum_vrb=sum(vcontig_relative_abundance))
reads_cenote_Class_Sum_abundance <- reads_cenote_Class_Sum_abundance %>% filter(Sample_type=="stool")
ClassSum_relative_abundance_plot <- ggplot(reads_cenote_Class_Sum_abundance, aes(x=Sample_ID, y=class_sum_vrb, fill=Class)) + 
  geom_bar(stat="identity") + 
  scale_fill_paletteer_d("palettetown::wurmple", na.value = "grey80") +
  facet_nested(cycle_number ~ Mock_Community + Sample_type, 
               scales = 'free_x', 
               nest_line = element_line(linetype=1)) +
  labs(x="Sample", y="relative abundance based on RPKM", 
       title = "relative abundance based on Class of viral contigs (1000-cycles vs 300-cycles)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background=element_blank(),
        ggh4x.facet.nestline=element_line(color="black"),
        strip.text.x=element_text(size=8))


# Add the reference relative abundance plot
i_ref_abundance_plot <- ref_df %>% select(Experiment, Sample_ID, Sample_type, Mock_Community, Treatment, Nuclease, reference_genome, ref_rpkm, mapped_ref_reads, total_reads, cycle_number)
i_ref_abundance_plot <- i_ref_abundance_plot %>% left_join(ref_reads_df %>% ungroup() %>% select (Experiment, Sample_ID, sum_replicate_rpkm, percent_viral_reads, cycle_number),
                                                           by = c("Sample_ID" = "Sample_ID", "Experiment"="Experiment", "cycle_number"="cycle_number")) #match the sum rpkm to their respective Sample_ID so that each input reference in each group has a matching total rpkm
i_ref_abundance_plot <- i_ref_abundance_plot %>% mutate(ref_relative_abundance=ref_rpkm/sum_replicate_rpkm) #calculate the relative abundance of each input reference
i_ref_abundance_plot <- i_ref_abundance_plot %>% filter(Mock_Community=="mock community spiked-in") %>% filter(Sample_type=="stool")
p1 <- ggplot(i_ref_abundance_plot, aes(x=Sample_ID, y=ref_relative_abundance, fill=reference_genome)) +
  geom_bar(stat="identity") + 
  scale_fill_paletteer_d("MoMAColors::Lupi") +
  scale_x_discrete(drop = T) +
  facet_nested(cycle_number ~ Mock_Community + Sample_type, 
               scales = 'free_x', 
               nest_line = element_line(linetype=1)) +
  labs(x="Sample", y="relative abundance based on RPKM", 
       title = "relative abundance based on reference genomes (1000-cycles vs 300-cycles)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust =1)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background=element_blank(),
        ggh4x.facet.nestline=element_line(color="black"))

