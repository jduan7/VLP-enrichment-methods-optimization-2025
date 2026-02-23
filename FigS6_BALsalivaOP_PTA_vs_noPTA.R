#Figure S6. Relative abundance of each reference virus after VirMock1 was spiked into saliva, OP wash, or BAL 
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

# Get the metadata
setwd("/Users/jduan/bushman/virome_methods/250616_nextseq")
metadata <- read_excel("metadata.xlsx")

# Get the output excel file for 250617_nextseq
setwd("/Users/jduan/bushman/virome_methods/250616_nextseq/R analysis")
ref_df <-read_excel("250617_nextseq_R_analysis.xlsx", sheet="ref_viral_analysis")
ref_reads_df <-read_excel("250617_nextseq_R_analysis.xlsx", sheet="ref_viral_percentage")
cenote_df <-read_excel("250617_nextseq_cenote_R_analysis.xlsx", sheet="cenote_result")
cenote_reads_df <-read_excel("250617_nextseq_cenote_R_analysis.xlsx", sheet="viral_percentage")

# Combine the two viral percentages tables into one
cenote_reads_df <- cenote_reads_df %>% mutate(contig_type="viral")
ref_reads_df <- ref_reads_df %>% mutate(contig_type="ref viral") %>% rename(percent_viral_reads=percent_ref_reads, sum_viral_reads=sum_viralref_reads)
viral_reads_abundance <- rbind(cenote_reads_df, ref_reads_df)

# Calculate relative abundance and make bar plots
i_ref_abundance_plot <- ref_df %>% select(Sample_ID, Experiment, Amplification, Treatment, Mock_Community, Sample_type, VC_aliquot_used, reference_genome, ref_rpkm, mapped_ref_reads, total_reads)
i_ref_abundance_plot <- i_ref_abundance_plot %>% 
  left_join(ref_reads_df %>% ungroup() %>% 
              select (Sample_ID, sum_replicate_rpkm, percent_viral_reads),
            by = c("Sample_ID" = "Sample_ID")) #match the sum rpkm to their respective group&type so that each input reference in each group has a matching total rpkm
i_ref_abundance_plot <- i_ref_abundance_plot %>% mutate(ref_relative_abundance=ref_rpkm/sum_replicate_rpkm) #calculate the relative abundance of each input reference

total_reads_data <- i_ref_abundance_plot %>% ungroup %>% group_by(Sample_ID) %>% summarize(TotalReads=first(total_reads), .groups="drop")
total_reads_data$TotalReads <- formatC(total_reads_data$TotalReads, format="e", digits=1) #change to scientific annotation

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



####################################### Save plots
# enable showtext
#library(showtext)
#font_add("Arial", regular = "/System/Library/Fonts/Supplemental/Arial.ttf")  # macOS path
#showtext_auto()  # enables showtext for all plots

# get the plots that are going to be saved
plots <- list(p_Exp1_ref)

# specify the destination folder
dest_folder <- "/Users/jduan/bushman/virome_methods/paper/supplemental"

# loop over plots and save as PDF
for (i in seq_along(plots)) {
  print(plots[[i]])   # ensures the plot is drawn
  ggsave(
    filename = file.path(dest_folder, paste0("FigS6_", i, ".pdf")),
    plot = plots[[i]],
    width = 6.5,
    height = 7.5,
    units = "in"
    # no need to specify device; showtext handles font embedding
  )
}

