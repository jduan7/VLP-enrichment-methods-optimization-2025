#Figure S1. Relative abundance of class of viruses as annotated by Cenote-Taker2 in saliva, stool or BAL that was analyzed by direct extraction of DNA (red) or RNA (blue), 
#or after viral particle enrichment and analysis using VP1, VP2, VP3, or VP4
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
mtg_cenote_df <- mtg_cenote_df %>% select(-Storage_Buffer) %>% mutate(Treatment="extraction")
mtg_cenote_reads_df <- mtg_cenote_reads_df %>% select(-Storage_Buffer) %>% mutate(Treatment="extraction")

# Combine the two viral percentages tables into one
cenote_df <- rbind(stool_VP1VP2_cenote_df, saliva_BAL_stoolVP3_cenote_df, mtg_cenote_df)
cenote_reads_df <- rbind(stool_VP1VP2_cenote_reads_df, saliva_BAL_stoolVP3_cenote_reads_df, mtg_cenote_reads_df)

# Calculate relative abundance of each viral contig in the respective sample
cenote_relative_abundance <- cenote_df %>% select(Sample_ID, Treatment, Mock_Community, Sample_type, Sequencing_type, organism, 
                                                  vcontig_rpkm, num_mapped_reads, num_total_reads, Class) %>%
  left_join(cenote_reads_df %>% ungroup() %>% select (Sample_ID, Treatment, Mock_Community, Sample_type, Sequencing_type, sum_replicate_rpkm, percent_viral_reads),
            by = c("Sample_ID" = "Sample_ID", "Sequencing_type" = "Sequencing_type", "Mock_Community"="Mock_Community", 
                   "Sample_type"="Sample_type", "Treatment"="Treatment")) #match the sum rpkm to their respective group&type so that each input reference in each group has a matching total rpkm
cenote_relative_abundance <- cenote_relative_abundance %>% 
  mutate(vcontig_relative_abundance=vcontig_rpkm/sum_replicate_rpkm) %>% #calculate the relative abundance of each input reference
  mutate(vcontig_relative_abundance = ifelse(is.na(vcontig_relative_abundance), 0, vcontig_relative_abundance)) # fill in the all-NA rows with 0 for later plotting

# based on grouping all the contigs of the same class together
cenote_Class_Sum_abundance <- cenote_relative_abundance %>% 
  group_by(Sample_ID, Treatment, Mock_Community, Sample_type, Sequencing_type, Class) %>% summarize(class_sum_vrb=sum(vcontig_relative_abundance))

# Make a plot
ClassSum_relative_abundance_plot <- ggplot(cenote_Class_Sum_abundance, aes(x=Sample_ID, y=class_sum_vrb, fill=Class)) + 
  geom_bar(stat="identity") + 
  scale_fill_paletteer_d("palettetown::wurmple", na.value = "grey80") +
  facet_nested(. ~ Mock_Community  + Sequencing_type + Treatment + Sample_type, 
               scales = 'free_x', 
               nest_line = element_line(linetype=1)) +
  labs(x="Sample", y="relative abundance based on RPKM", 
       title = "relative abundance based on Class of viral contigs (metagenomic vs virome)"
  ) +
  theme_minimal() +
  theme(strip.text.x = element_text(size=8), #change facet title font size
        text = element_text(family="Helvetica"), #change text font
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background=element_blank(),
        ggh4x.facet.nestline=element_line(color="black"))

