# Figure 2. Comparison of methods for purifying viral particles from stool.
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
setwd("/Users/jduan/bushman/virome_methods/paper/codes/codes for publication/")
metadata_250505 <- read_excel("metadata_250505.xlsx")
metadata_250617 <- read_excel("metadata_250617.xlsx")
metadata <- bind_rows(metadata_250505,metadata_250617)

# Get the output excel file for 250505_nextseq
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
cenote_reads_df <- cenote_reads_df %>% mutate(contig_type="All viruses (de novo assembled)")
ref_reads_df <- ref_reads_df %>% mutate(contig_type="VirMock1 only") %>% rename(percent_viral_reads=percent_ref_reads, sum_viral_reads=sum_viralref_reads)
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

############### Experiment 1 reference-based result
Exp1_ref_abundance <- i_ref_abundance_plot %>% filter(Experiment=="stool spike-in") %>% filter(Mock_Community=="mock community spiked-in")

p_Exp1_ref <- ggplot(Exp1_ref_abundance, aes(x=Sample_ID, y=ref_relative_abundance, fill=reference_genome)) +
  geom_bar(stat="identity") + 
  scale_fill_paletteer_d("MoMAColors::Lupi") +
  scale_x_discrete(drop = T) +
  facet_nested(. ~ Mock_Community + Sample_type + Treatment, 
               scales = 'free_x', 
               nest_line = element_line(linetype=1)) +
  labs(x="Sample", y="relative abundance based on RPKM", 
       #title = "relative abundance based on reference genomes (stool spike-in experiment)"
  ) +
  theme_minimal(base_family="Helvetica") +
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust =1)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background=element_blank(),
        ggh4x.facet.nestline=element_line(color="black"))

############### Experiment 1 percent viral reads result
# Make the total reads numeric
viral_reads_abundance$total_reads <- as.numeric(viral_reads_abundance$total_reads)

# Add the condition for plotting purpose
Exp1_reads_percentage <- viral_reads_abundance %>% filter(Experiment=="stool spike-in") %>% 
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
  scale_y_continuous(limits = c(0, 100)) +
  #scale_size_continuous(name = "Total reads", labels = scales::scientific, range = c(1, 6)) +
  scale_fill_brewer(palette = "Set2") + scale_color_brewer(palette = "Set2") +
  #guides(size = guide_legend(override.aes = list(shape = 21, fill = "black",color = "black"))) +
  labs(
    x = "Sample", y = "Percent in total reads", 
    title = "Viral reads percentage (stool spike-in experiment)",
    fill = "Alignment target", color = "Alignment target") +
  theme_bw(base_family="Helvetica") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_blank(),
    ggh4x.facet.nestline = element_line(color = "black")
  )

# Calculate for average percent viral reads and stdev for stool samples (w/spike-in)
Exp1_reads_percentage %>% filter(Sample_type=="stool") %>% filter(Mock_Community=="mock community spiked-in") %>% group_by(Treatment, contig_type) %>%
  summarize(mean_perc_viral=mean(percent_viral_reads),
            sd_perc_viral=sd(percent_viral_reads))

############ Experiment 1 cenote-based result

# Manually add back the samples that are missing due to the lack of assembled viral contigs from cenote result
i_cenote_df <- metadata %>% left_join(cenote_df, by = c("Sample_ID"="Sample_ID", "Experiment"="Experiment", "Mock_Community"="Mock_Community", "Sample_type"="Sample_type", "Treatment"="Treatment", "Nuclease"="Nuclease"))
cenote_reads_df <- cenote_reads_df %>%
  mutate(percent_viral_reads = ifelse(is.na(percent_viral_reads), 0, percent_viral_reads)) # replace the NA with 0 for plotting

Exp1_reads_cenote_Class_abundance <- i_cenote_df %>% filter(Experiment=="stool spike-in") %>% filter(Mock_Community=="mock community spiked-in") %>% 
  filter(Sample_type=="stool") %>%
  select(Sample_ID, Mock_Community, Sample_type, Experiment, Treatment, organism, vcontig_rpkm, num_mapped_reads, num_total_reads, Class) #make a table with all the reads data and the Class output by cenote
Exp1_reads_cenote_Class_abundance <- Exp1_reads_cenote_Class_abundance %>% left_join(cenote_reads_df %>% ungroup() %>% select (Experiment, Sample_ID, sum_replicate_rpkm, percent_viral_reads),
                                                                                     by = c("Sample_ID" = "Sample_ID", "Experiment" = "Experiment")) #match the sum rpkm to their respective group&type so that each input reference in each group has a matching total rpkm
Exp1_reads_cenote_Class_abundance <- Exp1_reads_cenote_Class_abundance %>% mutate(vcontig_relative_abundance=vcontig_rpkm/sum_replicate_rpkm) %>% #calculate the relative abundance of each input reference
  mutate(vcontig_relative_abundance = ifelse(is.na(vcontig_relative_abundance), 0, vcontig_relative_abundance)) # fill in the all-NA rows with 0 for later plotting

####### based on grouping all the contigs of the same class together
Exp1_reads_cenote_Class_Sum_abundance <- Exp1_reads_cenote_Class_abundance %>% 
  group_by(Sample_ID, Class, Mock_Community, Treatment, Sample_type) %>% summarize(class_sum_vrb=sum(vcontig_relative_abundance))

Exp1_Class_relative_abundance_plot <- ggplot(Exp1_reads_cenote_Class_Sum_abundance, aes(x=Sample_ID, y=class_sum_vrb, fill=Class)) + 
  geom_bar(stat="identity") + 
  scale_fill_paletteer_d("palettetown::wurmple", na.value = "grey80") +
  facet_nested(. ~ Mock_Community + Sample_type + Treatment, 
               scales = 'free_x', 
               nest_line = element_line(linetype=1)) +
  labs(x="Sample", y="relative abundance based on RPKM", 
       title = "relative abundance based on Class of viral contigs (stool spike-in experiment)"
  ) +
  theme_minimal(base_family="Helvetica") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background=element_blank(),
        ggh4x.facet.nestline=element_line(color="black"),
        strip.text.x=element_text(size=8))


# Put all panels together
(p_Exp1_ref/Exp1_Class_relative_abundance_plot) | q_Exp1
