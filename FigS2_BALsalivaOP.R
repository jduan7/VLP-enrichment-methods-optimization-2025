#Figure S2. Relative abundance of each reference virus after VirMock1 was spiked into saliva, OP wash, or BAL 
#and followed by VP4 and amplification by PTA or remained unamplified.
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

# Get the output excel file for 250617_nextseq
setwd("/Users/jduan/bushman/virome_methods/paper/codes/codes for publication/")
ref_df <-read_excel("250617_nextseq_R_analysis.xlsx", sheet="ref_viral_analysis")
ref_reads_df <-read_excel("250617_nextseq_R_analysis.xlsx", sheet="ref_viral_percentage")
ref_reads_df <- ref_reads_df %>% mutate(contig_type="ref viral") %>% rename(percent_viral_reads=percent_ref_reads, sum_viral_reads=sum_viralref_reads)

# Calculate relative abundance and make bar plots
i_ref_abundance_plot <- ref_df %>% select(Sample_ID, Experiment, Amplification, Treatment, Mock_Community, Sample_type, VC_aliquot_used, reference_genome, ref_rpkm, mapped_ref_reads, total_reads)
i_ref_abundance_plot <- i_ref_abundance_plot %>% 
  left_join(ref_reads_df %>% ungroup() %>% 
              select (Sample_ID, sum_replicate_rpkm, percent_viral_reads),
            by = c("Sample_ID" = "Sample_ID")) #match the sum rpkm to their respective group&type so that each input reference in each group has a matching total rpkm
i_ref_abundance_plot <- i_ref_abundance_plot %>% mutate(ref_relative_abundance=ref_rpkm/sum_replicate_rpkm) #calculate the relative abundance of each input reference


# Making relative abundance plots
Exp1_ref_abundance <- i_ref_abundance_plot %>% filter(grepl("saliva/OP/BAL spike-in", Experiment)) %>%
  filter(Mock_Community=="mock community spiked-in")

p_Exp1_ref <- ggplot(Exp1_ref_abundance, aes(x=Sample_ID, y=ref_relative_abundance, fill=reference_genome)) + 
  geom_bar(stat="identity") + 
  scale_fill_paletteer_d("MoMAColors::Lupi") +
  scale_x_discrete(drop = T) +
  facet_nested(. ~ Mock_Community + Treatment + Sample_type + Amplification, 
               scales = 'free_x', 
               nest_line = element_line(linetype=1)) +
  labs(x="Sample", y="relative abundance based on RPKM", 
       #title = "relative abundance based on reference genomes (saliva/OP wash/BAL spike-in experiment)",
       fill = "reference_genome", color = "reference_genome") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background=element_blank(),
        ggh4x.facet.nestline=element_line(color="black"))
