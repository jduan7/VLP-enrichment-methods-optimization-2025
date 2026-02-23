# This script performs Spearmam's Correlations on the nuclease treatment experiment
# This test detects any monotonic trend (even if nonlinear) and reports one p-value per virus

library(data.table)
library(readxl)
library(ggpubr)
library(openxlsx)
library(dplyr)
library(patchwork)
library(plotly)
library(paletteer)
library(ggrepel)

############Sequencing plot
# Get sequencing data
setwd("/Users/jduan/bushman/virome_methods/250429_nextseq/R analysis")
ref_viral_analysis<-read_excel("250505_nextseq_R_analysis.xlsx", sheet="ref_viral_analysis")
ref_viral_percentage<-read_excel("250505_nextseq_R_analysis.xlsx", sheet="ref_viral_percentage")

ref_viral_abundance <- ref_viral_analysis %>% left_join(ref_viral_percentage %>% ungroup() %>% select (Experiment, Sample_ID, sum_replicate_rpkm),
                                                        by = c("Sample_ID" = "Sample_ID", "Experiment"="Experiment")) #match the sum rpkm to their respective Sample_ID so that each input reference in each group has a matching total rpkm
ref_viral_abundance <- ref_viral_abundance %>% mutate(ref_relative_abundance=ref_rpkm/sum_replicate_rpkm) #calculate the relative abundance of each input reference

nuc_ref_viral_abundance <- ref_viral_abundance %>% filter(Experiment=="nuclease titration")

nuc_ref_viral_abundance$Nuclease <- factor(nuc_ref_viral_abundance$Nuclease,
                                           levels = c("none", "0.1X", "1X", "10X"))  # or your actual order

# Compute log2 fold change relative to untreated
nuc_norm_abundance <- nuc_ref_viral_abundance %>% filter(Mock_Community=="mock community spiked-in") %>%
  group_by(reference_genome) %>%
  mutate(untreated_mean = mean(ref_relative_abundance[Nuclease == "none"])) %>% #calculate the mean of untreated
  ungroup() %>%
  mutate(log2FC = log2((ref_relative_abundance+1e-6)/untreated_mean)) #add pseudo value to avoid zero, and normalize other values to the mean of untreated

# Encode nuclease numerically for regression
nuc_norm_abundance <- nuc_norm_abundance %>%
  mutate(Dose = case_when(
    Nuclease == "none" ~ 0,
    Nuclease == "0.1X"     ~ 0.1,
    Nuclease == "1X"       ~ 1,
    Nuclease == "10X"      ~ 10
  ),
  logDose = log10(Dose + 0.01))  # small offset to avoid log10(0)


# Compute mean and standard error per virus and condition
summary_nuc_norm_abundance <- nuc_norm_abundance %>%
  group_by(reference_genome, Dose) %>%
  summarise(
    mean_log2FC = mean(log2FC, na.rm = TRUE),
    se_log2FC = sd(log2FC, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(Dose_log = ifelse(Dose == 0, 0.01, Dose))  # tiny value for log scale


# Calculate the overall effect of nuclease for each virus
nuc_VC_seq_spearman <- nuc_norm_abundance %>%
  group_by(reference_genome) %>%
  summarise(
    rho = cor(log2FC, logDose, method = "spearman"),
    p_value = cor.test(log2FC, logDose, method = "spearman")$p.value
  ) %>%
  mutate(star=case_when(p_value < 0.001 ~ "***",
                        p_value < 0.01  ~ "**",
                        p_value < 0.05  ~ "*",
                        TRUE ~ "ns"))

# Make log2 fold-change line plot
nuc_VC_seq_plot <- ggplot(summary_nuc_norm_abundance, aes(x = Dose_log, y = mean_log2FC, color = reference_genome)) +
  geom_line(aes(group = reference_genome), size = 1) + # lines connecting means
  geom_point(size = 2) +  # mean points
  geom_errorbar(aes(ymin = mean_log2FC - se_log2FC,
                    ymax = mean_log2FC + se_log2FC),
                width = 0.05, 
                size = 0.5) + 
  scale_x_log10(
    limits=c(0.01, 20),
    breaks=c(0.01, 0.1, 1, 10),
    labels=c("untreated", "0.1X", "1X", "10X") # log10 scale for dose
  ) +                                       
  labs(x = "Nuclease treatment", y = "Log2 fold change compared to untreated") +
  scale_fill_paletteer_d("MoMAColors::Lupi") +
  scale_color_paletteer_d("MoMAColors::Lupi") +
  theme_minimal()



############ qPCR plot
# Get qpcr data
setwd("/Users/jduan/bushman/virome_methods/VP protocol test/nuclease_titration")
nuc_qpcr<-read_excel("250516_nuc_titration_qpcr_master.xlsx", sheet="for R analysis")
nuc_qpcr$Nuclease <- factor(nuc_qpcr$Nuclease,
                            levels = c("none", "0.1X", "1X", "10X"))  # or your actual order
nuc_qpcr_df <- nuc_qpcr %>% group_by(Group, Nuclease, Virus) %>% summarize(avg_VCN=mean(VCN)) %>% # get the average across the qpcr values
  rename(reference_genome=Virus)


# Compute log2 fold change relative to untreated
norm_nuc_qpcr <- nuc_qpcr_df %>% group_by(reference_genome) %>%
  mutate(untreated_mean = mean(avg_VCN[Nuclease == "none"])) %>% #calculate the mean of untreated
  ungroup() %>%
  mutate(log2FC = log2((avg_VCN+1e-6)/untreated_mean)) #add pseudo value to avoid zero, and normalize other values to the mean of untreated

# Encode nuclease numerically for regression
norm_nuc_qpcr <- norm_nuc_qpcr %>%
  mutate(Dose = case_when(
    Nuclease == "none" ~ 0,
    Nuclease == "0.1X"     ~ 0.1,
    Nuclease == "1X"       ~ 1,
    Nuclease == "10X"      ~ 10
  ),
  logDose = log10(Dose + 0.01))  # small offset to avoid log10(0)


# Compute mean and standard error per virus and condition
summary_norm_nuc_qpcr <- norm_nuc_qpcr %>%
  group_by(reference_genome, Dose) %>%
  summarise(
    mean_log2FC = mean(log2FC, na.rm = TRUE),
    se_log2FC = sd(log2FC, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(Dose_log = ifelse(Dose == 0, 0.01, Dose))  # tiny value for log scale


# Calculate the overall effect of nuclease for each virus
nuc_VC_qpcr_spearman <- norm_nuc_qpcr %>%
  group_by(reference_genome) %>%
  summarise(
    rho = cor(log2FC, logDose, method = "spearman"),
    p_value = cor.test(log2FC, logDose, method = "spearman")$p.value
  ) %>%
  mutate(star=case_when(p_value < 0.001 ~ "***",
                        p_value < 0.01  ~ "**",
                        p_value < 0.05  ~ "*",
                        TRUE ~ "ns"))

# Make log2 fold-change line plot
nuc_VC_qpcr_plot <- ggplot(summary_norm_nuc_qpcr, aes(x = Dose_log, y = mean_log2FC, color = reference_genome)) +
  geom_line(aes(group = reference_genome), size = 1) +               
  geom_point(size = 2) +                                  
  geom_errorbar(aes(ymin = mean_log2FC - se_log2FC,
                    ymax = mean_log2FC + se_log2FC),
                width = 0.05,                             
                size = 0.5) +                             
  scale_x_log10(
    limits=c(0.01, 20),
    breaks=c(0.01, 0.1, 1, 10),
    labels=c("untreated", "0.1X", "1X", "10X")
  ) +     
  labs(x = "Nuclease treatment", y = "Log2 fold change compared to untreated") +
  scale_fill_paletteer_d("MoMAColors::Lupi") +
  scale_color_paletteer_d("MoMAColors::Lupi") +
  theme_minimal()


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

# Compute log2 fold change relative to untreated
nuc_stool_norm_abundance <- nuc_stool_ref_viral_abundance %>% filter(Mock_Community=="mock community spiked-in") %>%
  group_by(reference_genome) %>%
  mutate(untreated_mean = mean(ref_relative_abundance[Nuclease == "none"])) %>% #calculate the mean of untreated
  ungroup() %>%
  mutate(log2FC = log2((ref_relative_abundance+1e-6)/untreated_mean)) #add pseudo value to avoid zero, and normalize other values to the mean of untreated

# Encode nuclease numerically for regression
nuc_stool_norm_abundance <- nuc_stool_norm_abundance %>%
  mutate(Dose = case_when(
    Nuclease == "none" ~ 0,
    Nuclease == "0.1X"     ~ 0.1,
    Nuclease == "1X"       ~ 1,
    Nuclease == "10X"      ~ 10
  ),
  logDose = log10(Dose + 0.01))  # small offset to avoid log10(0)

# Split between the last digit and the letter of Sample ID to create the sample and replicate columns
nuc_stool_norm_abundance <- nuc_stool_norm_abundance %>% 
  separate(Sample_ID, into=c("Sample", "replicate"), sep="(?<=\\d)(?=[A-Za-z])", remove=FALSE)
nuc_stool_norm_abundance <- nuc_stool_norm_abundance %>% mutate(replicate=factor(replicate)) # make replicate a factor

# Compute mean and standard error per virus and condition
summary_nuc_stool_norm_abundance <- nuc_stool_norm_abundance %>%
  group_by(reference_genome, Dose) %>%
  summarise(
    mean_log2FC = mean(log2FC, na.rm = TRUE),
    se_log2FC = sd(log2FC, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(Dose_log = ifelse(Dose == 0, 0.01, Dose))  # tiny value for log scale


# Calculate the overall effect of nuclease for each virus
nuc_VC_stool_seq_spearman <- nuc_stool_norm_abundance %>%
  group_by(reference_genome) %>%
  summarise(
    rho = cor(log2FC, logDose, method = "spearman"),
    p_value = cor.test(log2FC, logDose, method = "spearman")$p.value
  ) %>%
  mutate(star=case_when(p_value < 0.001 ~ "***",
                        p_value < 0.01  ~ "**",
                        p_value < 0.05  ~ "*",
                        TRUE ~ "ns"))

# Make log2 fold-change line plot
nuc_VC_stool_seq_plot <- ggplot(summary_nuc_stool_norm_abundance, aes(x = Dose_log, y = mean_log2FC, color = reference_genome)) +
  geom_line(aes(group = reference_genome), size = 1) +               
  geom_point(size = 2) +                                  
  geom_errorbar(aes(ymin = mean_log2FC - se_log2FC,
                    ymax = mean_log2FC + se_log2FC),
                width = 0.05,                             
                size = 0.5) +                             
  scale_x_log10(
    limits=c(0.01, 20),
    breaks=c(0.01, 0.1, 1, 10),
    labels=c("untreated", "0.1X", "1X", "10X")
  ) +                                       
  labs(x = "Nuclease treatment", y = "Log2 fold change compared to untreated") +
  scale_fill_paletteer_d("MoMAColors::Lupi") +
  scale_color_paletteer_d("MoMAColors::Lupi") +
  theme_minimal()



####################################### Save plots

# get the plots that are going to be saved
plots <- list(nuc_VC_seq_plot, nuc_VC_qpcr_plot, nuc_VC_stool_seq_plot)

# specify the destination folder
dest_folder <- "/Users/jduan/bushman/virome_methods/paper/manuscript/msystems_submission2/figure"
if(!dir.exists(dest_folder)) dir.create(dest_folder, recursive = TRUE)

# loop over plots and save as PDF
for (i in seq_along(plots)) {
  print(plots[[i]])   # ensures the plot is drawn
  ggsave(
    filename = file.path(dest_folder, paste0("Resub_Fig4_", i, ".pdf")),
    plot = plots[[i]],
    width = 6.5,
    height = 7.5,
    units = "in"
    # no need to specify device; showtext handles font embedding
  )
}


