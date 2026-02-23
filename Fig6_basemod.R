#Figure 6. Assessing the potential inhibitory effect of gluycosyl-hydroxymethylcytosine and hydroxymethylcytosine on recovery of modified T4 viral genomes. 
# Showing relative abundance of T4 and the lambda control by RPKM from sequencing normalized to that of the input genome copy number from qPCR
# Using Kruskal-Wallis to test whether any condition differs overall (running for each species separately).

library(data.table)
library(readxl)
library(ggpubr)
library(openxlsx)
library(dplyr)
library(patchwork)
library(plotly)
library(paletteer)


# Get the sequencing data
setwd("/Users/jduan/bushman/virome_methods/250910_nextseq/R analysis")
ref_rpkm_new <- read_excel("250918_nextseq_R_analysis_v3.xlsx", sheet="ref_viral_analysis")
read_percentage <- read_excel("250918_nextseq_R_analysis_v3.xlsx", sheet="ref_viral_percentage")
base_mod_ref_df <- ref_rpkm_new %>% filter(Experiment=="base modification", reference_genome %in% c("Enterobacteria phage T4","Enterobacteria phage lambda")) %>% 
  left_join(read_percentage %>% ungroup() %>% select (Sample_ID, sum_replicate_rpkm, percent_ref_reads),
            by = c("Sample_ID" = "Sample_ID")) %>%
  mutate(ref_relative_abundance=ref_rpkm/sum_replicate_rpkm) #calculate the relative abundance of each input reference

# Get the qPCR data
setwd("/Users/jduan/bushman/virome_methods/250910_nextseq/post-sequencing qPCR")
base_mod_qpcr <-read_excel("250922_T4+lambda_qPCR_master.xlsx") #get the qpcr data
base_mod_qpcr <- base_mod_qpcr %>% rename(reference_genome=Virus) %>% #change the column name to match up with the ref_viral table
  mutate(reference_genome=if_else(reference_genome=="lambda", "Enterobacteria phage lambda", reference_genome)) %>% 
  mutate(reference_genome=if_else(reference_genome=="T4", "Enterobacteria phage T4", reference_genome)) #change the virus name to match up with the ref_viral table


base_mod_qpcr_relative_abundance <- base_mod_qpcr %>% 
  group_by(Sample_ID) %>% #calculate relative abundance of lambda and T4 in each replicate
  mutate(qPCR_relative_abundance = VCN / sum(VCN))

base_mod_seq_qPCR_abundance <- base_mod_ref_df %>% 
  left_join(
    base_mod_qpcr_relative_abundance %>% 
      ungroup() %>% 
      select (Sample_ID, Group, reference_genome,
              qPCR_relative_abundance), #select the ones that will be joined to the table
    by = c("Sample_ID" = "Sample_ID", "reference_genome"="reference_genome") #left side of the equal sign is from the original table, right side of the equal sign is the columns that will be added
  ) 

base_mod_seq_qPCR_abundance <- base_mod_seq_qPCR_abundance %>% select(
  Experiment, Sample_ID, Mock_Community, reference_genome, ref_relative_abundance, 
  qPCR_relative_abundance) #make the table more concise

base_mod_seq_qPCR_abundance <- base_mod_seq_qPCR_abundance %>% 
  mutate(
    qpcr_ra_norm=qPCR_relative_abundance/qPCR_relative_abundance,
    seq_qpcr_ra_norm=ref_relative_abundance/qPCR_relative_abundance,
  ) #normalization based on qpcr relative abundance

base_mod_seq_qPCR_abundance_plot <- base_mod_seq_qPCR_abundance %>%
  pivot_longer(cols=c(ref_relative_abundance,
                      qPCR_relative_abundance, 
                      qpcr_ra_norm, 
                      seq_qpcr_ra_norm, 
  ), 
  names_to = "data_type", values_to = "relative_abundance") %>%
  filter(is.finite(relative_abundance))#reshape the data so that there is a column that is either qpcr_ra_norm or seq_ra_norm

# Calculate the normalized relative abundance of T4 and lambda from qPCR and sequencing data
norm_base_mod_seq_qPCR_abundance_plot <- base_mod_seq_qPCR_abundance_plot %>% filter(data_type %in% c("qpcr_ra_norm", "seq_qpcr_ra_norm")) 

  
norm_base_mod_seq_qPCR_abundance_plot <- norm_base_mod_seq_qPCR_abundance_plot %>% 
  filter(Mock_Community != "lambda", Mock_Community != "SM buffer", data_type=="seq_qpcr_ra_norm")

# Perform statistical test - H0: The median normalized relative abundance is the same across T4C, T4hmC, and T4ghmC.
kruskal_T4 <- norm_base_mod_seq_qPCR_abundance_plot %>% 
  filter(reference_genome=="Enterobacteria phage T4") #filter for just the T4 data only

# Perform Kruskal-Wallis and get the p-value
kw_p_T4 <- kruskal.test(relative_abundance ~ Mock_Community, data = kruskal_T4)$p.value

# Make barplot (noramlized to qpcr)
summary_seq_qpcr_bar_data <- norm_base_mod_seq_qPCR_abundance_plot %>% 
  group_by(Mock_Community, reference_genome) %>%
  summarise(mean_ra = mean(relative_abundance),
            se_ra = sd(relative_abundance) / sqrt(n()),
            .groups = "drop")

# Combine the bar plot and the scatter plot
plot_VCNnorm<-ggplot(summary_seq_qpcr_bar_data,
                     aes(x = factor(Mock_Community, levels=c("lambda + T4C", "lambda + T4hmC", "lambda + T4ghmC")), 
                         y = mean_ra, fill = reference_genome)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), alpha = 0.7) +
  geom_errorbar(aes(ymin = mean_ra - se_ra, ymax = mean_ra + se_ra),
                position = position_dodge(width = 0.75), width = 0.2) +
  geom_point(data=norm_base_mod_seq_qPCR_abundance_plot,
             aes(x = Mock_Community, y = relative_abundance, color= reference_genome),
             size = 2, 
             position = position_jitterdodge(jitter.width = 0.1,dodge.width=0.9)) +
  scale_y_continuous(limits=c(0,2)) +
  labs(x = "Mock Community", y = "Normalized relative abundance (VCN based)", fill = "Reference genome", color ="Reference genome") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust =1))

####################################### Save plots

# get the plots that are going to be saved
plots <- list(plot_VCNnorm)

# specify the destination folder
dest_folder <- "/Users/jduan/bushman/virome_methods/paper/manuscript/msystems_submission2/figure"
if(!dir.exists(dest_folder)) dir.create(dest_folder, recursive = TRUE)

# loop over plots and save as PDF
for (i in seq_along(plots)) {
  print(plots[[i]])   # ensures the plot is drawn
  ggsave(
    filename = file.path(dest_folder, paste0("Resub_Fig6_", i, ".pdf")),
    plot = plots[[i]],
    width = 6.5,
    height = 7.5,
    units = "in"
    # no need to specify device; showtext handles font embedding
  )
}

