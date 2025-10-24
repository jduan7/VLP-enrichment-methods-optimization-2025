#Figure 4. Assessing the stability of viral genomes within viral particles in the presence of nucleases.
library(data.table)
library(readxl)
library(ggpubr)
library(openxlsx)
library(dplyr)
library(patchwork)
library(plotly)
library(paletteer)


############Sequencing plot
# Get sequencing data
setwd("/Users/jduan/bushman/virome_methods/paper/codes/codes for publication/")
ref_viral_analysis<-read_excel("250505_nextseq_R_analysis.xlsx", sheet="ref_viral_analysis")
ref_viral_percentage<-read_excel("250505_nextseq_R_analysis.xlsx", sheet="ref_viral_percentage")

ref_viral_abundance <- ref_viral_analysis %>% left_join(ref_viral_percentage %>% ungroup() %>% select (Experiment, Sample_ID, sum_replicate_rpkm),
                                                        by = c("Sample_ID" = "Sample_ID", "Experiment"="Experiment")) #match the sum rpkm to their respective Sample_ID so that each input reference in each group has a matching total rpkm
ref_viral_abundance <- ref_viral_abundance %>% mutate(ref_relative_abundance=ref_rpkm/sum_replicate_rpkm) #calculate the relative abundance of each input reference

nuc_ref_viral_abundance <- ref_viral_abundance %>% filter(Experiment=="nuclease titration")

nuc_ref_viral_abundance$Nuclease <- factor(nuc_ref_viral_abundance$Nuclease,
                                           levels = c("none", "0.1X", "1X", "10X"))  # or your actual order

summary_nuc_ref_viral_abundance <- nuc_ref_viral_abundance %>%
  filter(is.finite(ref_relative_abundance)) %>%
  group_by(Nuclease, reference_genome) %>%
  summarise(mean_abundance = mean(ref_relative_abundance), .groups = "drop")

# Calculate and plot for relative abundance normalized to the no treatment
norm_abundance <- summary_nuc_ref_viral_abundance %>% filter(Nuclease=="none") %>% rename(no_nuclease_abundance=mean_abundance) #get the no_nuclease relative abundance
nuc_norm_abundance <- nuc_ref_viral_abundance %>% left_join(norm_abundance %>% ungroup() %>% select (reference_genome,no_nuclease_abundance),
                                                            by = c("reference_genome" = "reference_genome")) #join the table by matching the no_nuclease relative abundance value for each nuclease group according to reference genome
nuc_norm_abundance <- nuc_norm_abundance %>% mutate(normalized_abundance = ref_relative_abundance/no_nuclease_abundance) #normalize the relative abundance of each nuclease treatment group by the no_treatment relative abundance value

summary_norm_abundance <- nuc_norm_abundance %>%
  filter(is.finite(ref_relative_abundance)) %>%
  group_by(Nuclease, reference_genome) %>%
  summarise(mean_norm_abundance = mean(normalized_abundance), .groups = "drop") #make a summary table of the mean of normalized abundance for triplicates

# Make plot
nuc_VC_seq_plot <- ggplot(nuc_norm_abundance %>% filter(is.finite(normalized_abundance)),
                          aes(x = Nuclease, y = normalized_abundance, fill = reference_genome)) +
  geom_jitter(aes(color = reference_genome),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.25),
              size = 1.5, alpha = 0.7) +
  geom_line(data = summary_norm_abundance %>% filter(mean_norm_abundance>0),
            aes(x = Nuclease, y = mean_norm_abundance, group = reference_genome, color = reference_genome),
            position = position_dodge(width = 0.25),
            size = 0.75, alpha=0.5) +
  geom_errorbar(data = summary_norm_abundance %>% filter(mean_norm_abundance>0),
                inherit.aes=FALSE,
                aes(x = Nuclease,
                    ymin = mean_norm_abundance, 
                    ymax = mean_norm_abundance,
                    color=reference_genome),
                position = position_dodge(width = 0.25), 
                linewidth=1,
                width=0.6) +
  #facet_wrap(~reference_genome, ncol=1) +
  scale_fill_paletteer_d("MoMAColors::Lupi") +
  scale_color_paletteer_d("MoMAColors::Lupi") +
  scale_y_log10() +
  labs(x = "Nuclease Treatment", y = "Relative Abundance normalized to no treatment (log10 scale)", 
       title = "nuclease titration on mock community only \n(sequencing result)",
       fill = "Reference Genome", color = "Reference Genome") + 
  theme_minimal(base_family="Helvetica") +
  theme(strip.text.x = element_text(size = 10),
        strip.background = element_blank(),
        axis.text.y = element_text(size=10))





############ qPCR plot
# Get qpcr data
setwd("/Users/jduan/bushman/virome_methods/paper/codes/codes for publication/")
nuc_qpcr<-read_excel("250516_nuc_titration_qpcr.xlsx")
nuc_qpcr$Nuclease <- factor(nuc_qpcr$Nuclease,
                            levels = c("none", "0.1X", "1X", "10X"))  # or your actual order
nuc_qpcr_df <- nuc_qpcr %>% group_by(Group, Nuclease, Virus) %>% summarize(avg_VCN=mean(VCN)) # get the average across the qpcr values
summary_nuc_qpcr <- nuc_qpcr_df %>% group_by(Nuclease, Virus) %>% summarize(avg_VCN_across_tri=mean(avg_VCN), sd_VCN=sd(avg_VCN)) #get the avg VCN across triplicates and the stdev

# Calculate and plot for VCN normalized to the no treatment
norm_nuc_qpcr <- summary_nuc_qpcr %>% filter(Nuclease=="none") %>% rename(no_nuclease_VCN=avg_VCN_across_tri) #get the no_nuclease avg VCN across triplicates
norm_nuc_qpcr_df <- nuc_qpcr_df %>% left_join(norm_nuc_qpcr %>% ungroup() %>% select (Virus, no_nuclease_VCN),
                                              by = c("Virus" = "Virus")) #join the table by matching the no_nuclease avg VCN for each nuclease group according to reference genome
norm_nuc_qpcr_df <- norm_nuc_qpcr_df %>% mutate(normalized_VCN = avg_VCN/no_nuclease_VCN) #normalize the avg VCN of each nuclease treatment group by the no_treatment avg VCN value
summary_norm_VCN <- norm_nuc_qpcr_df %>%
  group_by(Nuclease, Virus) %>%
  summarise(mean_norm_VCN = mean(normalized_VCN), .groups = "drop") # make a summary table of the mean of normalized VCN for triplicates
summary_norm_VCN$percent_VCN <- paste0(round(summary_norm_VCN$mean_norm_VCN*100,1), "%") # convert to percentage

# Make plot
nuc_VC_qpcr_plot <- ggplot(norm_nuc_qpcr_df,
                           aes(x = Nuclease, y = normalized_VCN, fill = Virus)) +
  geom_jitter(aes(color = Virus),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.25),
              size = 1.5, alpha = 0.7) +
  geom_line(data = summary_norm_VCN,
            aes(x = Nuclease, y = mean_norm_VCN, group = Virus, color = Virus),
            position = position_dodge(width = 0.25),
            size = 0.75, alpha=0.5) +
  geom_errorbar(data = summary_norm_VCN %>% filter(mean_norm_VCN>0),
                inherit.aes=FALSE,
                aes(x = Nuclease, 
                    ymin = mean_norm_VCN, 
                    ymax = mean_norm_VCN,
                    color=Virus),
                position = position_dodge(width = 0.25), 
                linewidth=1,
                width=0.6) +
  #facet_wrap(~Virus, ncol=1) +
  scale_y_log10() +
  scale_fill_paletteer_d("MoMAColors::Lupi") +
  scale_color_paletteer_d("MoMAColors::Lupi") +
  labs(x = "Nuclease Treatment", y = "Virus Copy Number normalized to no treatment (log10 scale)", 
       title = "nuclease titration on mock community only \n(qPCR result)",
       fill = "Virus", color = "Virus") + 
  theme_minimal(base_family="Helvetica") +
  theme(strip.text.x = element_text(size = 10),
        strip.background = element_blank(),
        axis.text.y = element_text(size=10))




############Sequencing plot (w/stool spike-in)
# Get sequencing data
setwd("/Users/jduan/bushman/virome_methods/250616_nextseq/R analysis")
stool_ref_viral_analysis<-read_excel("250617_nextseq_R_analysis.xlsx", sheet="ref_viral_analysis")
stool_ref_viral_percentage<-read_excel("250617_nextseq_R_analysis.xlsx", sheet="ref_viral_percentage")

stool_ref_viral_abundance <- stool_ref_viral_analysis %>% left_join(stool_ref_viral_percentage %>% ungroup() %>% select (Experiment, Sample_ID, sum_replicate_rpkm),
                                                                    by = c("Sample_ID" = "Sample_ID", "Experiment"="Experiment")) #match the sum rpkm to their respective Sample_ID so that each input reference in each group has a matching total rpkm
stool_ref_viral_abundance <- stool_ref_viral_abundance %>% mutate(ref_relative_abundance=ref_rpkm/sum_replicate_rpkm) #calculate the relative abundance of each input reference

nuc_stool_ref_viral_abundance <- stool_ref_viral_abundance %>% filter(Experiment=="nuclease titration w/ stool spike-in")

nuc_stool_ref_viral_abundance$Nuclease <- factor(nuc_stool_ref_viral_abundance$Nuclease,
                                                 levels = c("none", "0.1X", "1X", "10X"))  # or your actual order

summary_nuc_stool_ref_viral_abundance <- nuc_stool_ref_viral_abundance %>%
  filter(is.finite(ref_relative_abundance)) %>%
  group_by(Nuclease, reference_genome) %>%
  summarise(mean_abundance = mean(ref_relative_abundance), .groups = "drop")


#### Calculate and plot for relative abundance normalized to the no treatment
stool_norm_abundance <- summary_nuc_stool_ref_viral_abundance %>% filter(Nuclease=="none") %>% rename(no_nuclease_abundance=mean_abundance) #get the no_nuclease relative abundance
nuc_stool_norm_abundance <- nuc_stool_ref_viral_abundance %>% left_join(stool_norm_abundance %>% ungroup() %>% select (reference_genome,no_nuclease_abundance),
                                                                        by = c("reference_genome" = "reference_genome")) #join the table by matching the no_nuclease relative abundance value for each nuclease group according to reference genome
nuc_stool_norm_abundance <- nuc_stool_norm_abundance %>% mutate(normalized_abundance = ref_relative_abundance/no_nuclease_abundance) #normalize the relative abundance of each nuclease treatment group by the no_treatment relative abundance value

summary_stool_norm_abundance <- nuc_stool_norm_abundance %>%
  filter(is.finite(ref_relative_abundance)) %>%
  group_by(Nuclease, reference_genome) %>%
  summarise(mean_norm_abundance = mean(normalized_abundance), .groups = "drop") #make a summary table of the mean of normalized abundance for triplicates

nuc_VC_stool_seq_plot <- ggplot(nuc_stool_norm_abundance %>% filter(is.finite(normalized_abundance)),
                                aes(x = Nuclease, y = normalized_abundance, fill = reference_genome)) +
  geom_jitter(aes(color = reference_genome),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.25),
              size = 1.5, alpha = 0.7) +
  geom_line(data = summary_stool_norm_abundance %>% filter(mean_norm_abundance>0),
            aes(x = Nuclease, y = mean_norm_abundance, group = reference_genome, color = reference_genome),
            position = position_dodge(width = 0.25),
            size = 0.75, alpha=0.5) +
  geom_errorbar(data = summary_stool_norm_abundance %>% filter(mean_norm_abundance>0),
                inherit.aes=FALSE,
                aes(x = Nuclease, 
                    ymin = mean_norm_abundance, 
                    ymax = mean_norm_abundance,
                    color=reference_genome),
                position = position_dodge(width = 0.25), 
                linewidth=1,
                width=0.6) +
  #facet_wrap(~reference_genome, ncol=1) +
  scale_y_log10() +
  scale_fill_paletteer_d("MoMAColors::Lupi") +
  scale_color_paletteer_d("MoMAColors::Lupi") +
  labs(x = "Nuclease Treatment", y = "Relative Abundance normalized to no treatment (log10 scale)", 
       title = "nuclease titration on mock community spiked into stool \n(sequencing result)",
       fill = "Reference Genome", color = "Reference Genome") + 
  theme_minimal(base_family="Helvetica") +
  theme(strip.text.x = element_text(size = 10),
        strip.background = element_blank(),
        axis.text.y = element_text(size=10))

