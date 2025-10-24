#Figure1. Comparison of viral sequence recovery after virus particle enrichment followed by sequencing, 
#versus direct metagenomic sequencing of total DNA or RNA. 

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

# Get the output excel file for 250505_nextseq
setwd("/Users/jduan/bushman/virome_methods/paper/codes/codes for publication/")
stool_VP1VP2_cenote_df <-read_excel("250505_nextseq_cenote_R_analysis.xlsx", sheet="cenote_result") %>%
  filter(Experiment=="stool spike-in"|Experiment=="saliva/OP/BAL spike-in", 
         Sample_type=="BAL"|Sample_type=="saliva"|Sample_type=="stool",
         Mock_Community=="no mock community")
stool_VP1VP2_cenote_reads_df <-read_excel("250505_nextseq_cenote_R_analysis.xlsx", sheet="viral_percentage") %>%
  filter(Experiment=="stool spike-in"|Experiment=="saliva/OP/BAL spike-in", 
         Sample_type=="BAL"|Sample_type=="saliva"|Sample_type=="stool",
         Mock_Community=="no mock community")

# Get the output excel file for 250617_nextseq
setwd("/Users/jduan/bushman/virome_methods/paper/codes/codes for publication/")
saliva_BAL_stoolVP3_cenote_df <-read_excel("250617_nextseq_cenote_R_analysis.xlsx", sheet="cenote_result") %>%
  filter(Experiment=="stool spike-in"|Experiment=="saliva/OP/BAL spike-in", 
         Sample_type=="BAL"|Sample_type=="saliva"|Sample_type=="stool",
         Mock_Community=="no mock community")
saliva_BAL_stoolVP3_cenote_reads_df <-read_excel("250617_nextseq_cenote_R_analysis.xlsx", sheet="viral_percentage") %>%
  filter(Experiment=="stool spike-in"|Experiment=="saliva/OP/BAL spike-in", 
         Sample_type=="BAL"|Sample_type=="saliva"|Sample_type=="stool",
         Mock_Community=="no mock community")

# Get the output excel file for metagenomic
setwd("/Users/jduan/bushman/virome_methods/paper/codes/codes for publication/")
mtg_cenote_df <- read_excel("CHOPMC484-485_cenote_R_analysis.xlsx", sheet="cenote_result") %>%
  filter(Storage_Buffer=="neat", 
         Sample_type=="BAL"|Sample_type=="saliva"|Sample_type=="stool",
         Mock_Community=="no mock community")
mtg_cenote_reads_df <-read_excel("CHOPMC484-485_cenote_R_analysis.xlsx", sheet="viral_percentage") %>%
  filter(Storage_Buffer=="neat", 
         Sample_type=="BAL"|Sample_type=="saliva"|Sample_type=="stool",
         Mock_Community=="no mock community") %>%
  mutate(Sample_ID=gsub("\\.", "_", Sample_ID))

# Format the tables
stool_VP1VP2_cenote_df <- stool_VP1VP2_cenote_df %>% select(-Experiment, -Nuclease) %>% mutate(Sequencing_type="virome")
stool_VP1VP2_cenote_reads_df <- stool_VP1VP2_cenote_reads_df %>% select(-Experiment, -Nuclease) %>% mutate(Sequencing_type="virome")
saliva_BAL_stoolVP3_cenote_df <- saliva_BAL_stoolVP3_cenote_df %>% select(-Experiment, -Nuclease, -Amplification, -VC_aliquot_used) %>% mutate(Sequencing_type="virome")
saliva_BAL_stoolVP3_cenote_reads_df <- saliva_BAL_stoolVP3_cenote_reads_df %>% select(-Experiment, -Nuclease, -Amplification, -VC_aliquot_used) %>% mutate(Sequencing_type="virome")
mtg_cenote_df <- mtg_cenote_df %>% select(-Storage_Buffer) %>% mutate(Treatment="direct extraction")
mtg_cenote_reads_df <- mtg_cenote_reads_df %>% select(-Storage_Buffer) %>% mutate(Treatment="direct extraction")

# Combine the two viral percentages tables into one
cenote_df <- rbind(stool_VP1VP2_cenote_df, saliva_BAL_stoolVP3_cenote_df, mtg_cenote_df)
cenote_reads_df <- rbind(stool_VP1VP2_cenote_reads_df, saliva_BAL_stoolVP3_cenote_reads_df, mtg_cenote_reads_df)


############### Percent viral reads result
# Make the total reads numeric
cenote_reads_df$total_reads <- as.numeric(cenote_reads_df$total_reads)

viral_reads_abundance <- cenote_reads_df %>% mutate(condition=paste(Sample_type, Treatment, sep=" | "))

# Plot for stool viral percentage comparison
stool_viral_reads_abundance <- viral_reads_abundance %>% filter(Sample_type=="stool")
q_stool <- ggplot(stool_viral_reads_abundance, aes(x = condition, y = percent_viral_reads)) +
  geom_jitter(aes(fill = Sequencing_type, color = Sequencing_type),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
              size=3, alpha = 0.7, shape = 21) +
  stat_summary(fun = mean, 
               geom = "crossbar",
               aes(group = Sequencing_type),
               position = position_dodge(width = 0.75), width = 0.25, fatten = 1.5, color = "black", show.legend = FALSE) +
  scale_x_discrete(drop = TRUE) +
  facet_nested(. ~ Mock_Community + Treatment + Sample_type, 
               scales = "free_x", nest_line = element_line(linetype = 1)) +
  scale_y_continuous(limits = c(0, 100)) +
  #scale_size_continuous(name = "Total reads", labels = scales::scientific, range = c(1, 6)) +
  scale_fill_brewer(palette = "Set1") + scale_color_brewer(palette = "Set1") +
  guides(size = guide_legend(override.aes = list(shape = 21, fill = "black",color = "black"))) +
  labs(
    x = "Sample", y = "Percent in total reads", 
    #title = "Viral reads percentage (metagenomic vs virome)",
    fill = "Sequencing type", color = "Sequencing type") +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_blank(),
    ggh4x.facet.nestline = element_line(color = "black")
  )

# Plot for saliva & BAL viral percentage comparison
SalivaBAL_viral_reads_abundance <- viral_reads_abundance %>% filter(Sample_type %in% c("saliva","BAL"))
q_SalivaBAL <- ggplot(SalivaBAL_viral_reads_abundance, aes(x = condition, y = percent_viral_reads)) +
  geom_jitter(aes(fill = Sequencing_type, color = Sequencing_type),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
              size=3, alpha = 0.7, shape = 21) +
  stat_summary(fun = mean, 
               geom = "crossbar",
               aes(group = Sequencing_type),
               position = position_dodge(width = 0.75), width = 0.25, fatten = 1.5, color = "black", show.legend = FALSE) +
  scale_x_discrete(drop = TRUE) +
  facet_nested(. ~ Mock_Community + Treatment + Sample_type, 
               scales = "free_x", nest_line = element_line(linetype = 1)) +
  scale_y_continuous(limits = c(0, 100)) +
  #scale_size_continuous(name = "Total reads", labels = scales::scientific, range = c(1, 6)) +
  scale_fill_brewer(palette = "Set1") + scale_color_brewer(palette = "Set1") +
  guides(size = guide_legend(override.aes = list(shape = 21, fill = "black",color = "black"))) +
  labs(
    x = "Sample", y = "Percent in total reads", 
    #title = "Viral reads percentage (metagenomic vs virome)",
    fill = "Sequencing type", color = "Sequencing type") +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_blank(),
    ggh4x.facet.nestline = element_line(color = "black")
  )

