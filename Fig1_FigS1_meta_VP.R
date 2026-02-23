#Figure1 and Fig S1. Comparison of viral sequence recovery after virus particle enrichment followed by sequencing, 
#versus direct metagenomic sequencing of total DNA or RNA. (Cenote-Taker2 based viral annotation)

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
setwd("/Users/jduan/bushman/virome_methods/250429_nextseq/R analysis")
stool_VP1VP2_cenote_df <-read_excel("250505_nextseq_cenote_R_analysis.xlsx", sheet="cenote_result") %>%
  filter(Experiment=="stool spike-in"|Experiment=="saliva/OP/BAL spike-in", 
         Sample_type=="BAL"|Sample_type=="saliva"|Sample_type=="stool",
         Mock_Community=="no mock community")
stool_VP1VP2_cenote_reads_df <-read_excel("250505_nextseq_cenote_R_analysis.xlsx", sheet="viral_percentage") %>%
  filter(Experiment=="stool spike-in"|Experiment=="saliva/OP/BAL spike-in", 
         Sample_type=="BAL"|Sample_type=="saliva"|Sample_type=="stool",
         Mock_Community=="no mock community")

# Get the output excel file for 250617_nextseq
setwd("/Users/jduan/bushman/virome_methods/250616_nextseq/R analysis")
saliva_BAL_stoolVP3_cenote_df <-read_excel("250617_nextseq_cenote_R_analysis.xlsx", sheet="cenote_result") %>%
  filter(Experiment=="stool spike-in"|Experiment=="saliva/OP/BAL spike-in", 
         Sample_type=="BAL"|Sample_type=="saliva"|Sample_type=="stool",
         Mock_Community=="no mock community")
saliva_BAL_stoolVP3_cenote_reads_df <-read_excel("250617_nextseq_cenote_R_analysis.xlsx", sheet="viral_percentage") %>%
  filter(Experiment=="stool spike-in"|Experiment=="saliva/OP/BAL spike-in", 
         Sample_type=="BAL"|Sample_type=="saliva"|Sample_type=="stool",
         Mock_Community=="no mock community")

# Get the output excel file for metagenomic
setwd("/Users/jduan/bushman/virome_methods/metagenomic sequencing/R analysis")
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
  facet_nested(. ~ Mock_Community + Sample_type + Sequencing_type + Treatment, 
               scales = 'free_x', 
               nest_line = element_line(linetype=1)) +
  labs(x="Sample", y="relative abundance based on RPKM", 
       #title = "relative abundance based on Class of viral contigs (metagenomic vs virome)"
  ) +
  theme_minimal() +
  theme(strip.text.x = element_text(size=8), #change facet title font size
        text = element_text(family="Helvetica"), #change text font
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background=element_blank(),
        ggh4x.facet.nestline=element_line(color="black"))


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
  # Optional: Better color palettes
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
  # Optional: Better color palettes
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


####### summarize for the average %viral reads of each group
stool_viral_reads_abundance %>% group_by(Treatment, Sequencing_type, Sample_type) %>% 
  summarize(avg_percent_viral_reads=mean(percent_viral_reads),
            sd_perc_viral=sd(percent_viral_reads))
SalivaBAL_viral_reads_abundance %>% group_by(Treatment, Sequencing_type, Sample_type) %>% 
  summarize(avg_percent_viral_reads=mean(percent_viral_reads),
            sd_perc_viral=sd(percent_viral_reads))

####################################### Save plots
# enable showtext
#library(showtext)
#font_add("Arial", regular = "/System/Library/Fonts/Supplemental/Arial.ttf")  # macOS path
#showtext_auto()  # enables showtext for all plots

# get the plots that are going to be saved
plots <- list(q_stool, q_SalivaBAL)
sup_plots <- list(ClassSum_relative_abundance_plot)

# specify the destination folder
dest_folder <- "/Users/jduan/bushman/virome_methods/paper/figures_pdf"
if(!dir.exists(dest_folder)) dir.create(dest_folder, recursive = TRUE)

# loop over plots and save as PDF
for (i in seq_along(plots)) {
  print(plots[[i]])   # ensures the plot is drawn
  ggsave(
    filename = file.path(dest_folder, paste0("Fig1_", i, ".pdf")),
    plot = plots[[i]],
    width = 6.5,
    height = 7.5,
    units = "in"
    # no need to specify device; showtext handles font embedding
  )
}

dest_folder <- "/Users/jduan/bushman/virome_methods/paper/supplemental"
for (i in seq_along(sup_plots)) {
  print(sup_plots[[i]])   # ensures the plot is drawn
  ggsave(
    filename = file.path(dest_folder, paste0("FigS1_", i, ".pdf")),
    plot = sup_plots[[i]],
    width = 6.5,
    height = 7.5,
    units = "in"
    # no need to specify device; showtext handles font embedding
  )
}

