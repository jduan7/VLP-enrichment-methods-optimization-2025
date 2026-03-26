# This script plots sequencing depth vs average contig length for viral contigs of each sample
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

# Import the contig length data
setwd("/Users/jduan/bushman/virome_methods/paper/manuscript/msystems_submission2/contiglens/viral_contig_lens")
cenote_250221 <- read_excel("250221_miniseq_cenote_R_analysis.xlsx", sheet="cenote_result") %>% 
  select(sample, contig, Mock_Community, virus_seq_length) %>% rename(Sample_ID=sample) %>% mutate(Experiment="Amplification")
cenote_250228 <- read_excel("250228_nextseq_v2_cenote_R_analysis.xlsx", sheet="cenote_result") %>% 
  select(sample, contig, Mock_Community, virus_seq_length) %>% rename(Sample_ID=sample) %>% mutate(Experiment="Amplification")
cenote_250505 <- read_excel("250505_nextseq_cenote_R_analysis.xlsx", sheet="cenote_result") %>% 
  select(Sample_ID, contig, Experiment, Mock_Community, virus_seq_length)
cenote_250617 <- read_excel("250617_nextseq_cenote_R_analysis.xlsx", sheet="cenote_result") %>% 
  select(Sample_ID, contig, Experiment, Mock_Community, virus_seq_length)
cenote_250918 <- read_excel("250918_nextseq_cenote_R_analysis.xlsx", sheet="cenote_result") %>% 
  select(Index, contig, Experiment, Mock_Community, virus_seq_length) %>%
  rename(Sample_ID=Index) %>% filter(Experiment=="base modification")
cenote_260116 <- read_excel("260116_nextseq_cenote_R_analysis.xlsx", sheet="cenote_result") %>% 
  select(Sample_ID, contig, Experiment, Mock_Community, virus_seq_length)

vmock_vcontig_df <- rbind(cenote_250221, cenote_250228, cenote_250505, cenote_250617, cenote_250918, cenote_260116) %>%
  filter(Mock_Community!="no mock community") %>% filter(Mock_Community!="SM buffer")


# Import the raw reads data
setwd("/Users/jduan/Downloads")
reads_df <-read_excel("mSystems00188-26-Table_S1_to_S15.xlsx", sheet="S15_sample_metadata", skip=1) %>% filter(Read_Length_bp =="300")

# Filter for those with spiked-in community
reads_df_vmock <- reads_df %>% filter(Mock_Community!="no mock community") %>% filter(Mock_Community!="SM buffer")

# Assign join_id for later left_join
reads_df_vmock <- reads_df_vmock %>%
  mutate(join_id = if_else(Experiment %in% c("nuclease titration", "stool spike-in"),
                           Sample_ID,
                           Index))

# Match with the raw reads data
depth_vcontig_df <- reads_df_vmock %>% select(join_id, R1_raw_read_count) %>%
  left_join(
    vmock_vcontig_df,
    by = c("join_id"="Sample_ID")
  ) %>% 
  rename(Sample_ID=join_id, seq_depth=R1_raw_read_count)

####### Compare sequencing depth (total reads) vs mean contig length of each sample
summary_depth_vcontig_df <- depth_vcontig_df %>% group_by(Sample_ID, Experiment, seq_depth) %>% 
  summarize(mean_sample_vlen=mean(virus_seq_length, na.rm = TRUE),
            se_sample_vlen=sd(virus_seq_length, na.rm = TRUE)/ sqrt(n()),
            .groups="drop") %>%
  filter(!is.nan(mean_sample_vlen)) # excluding those with no viral contigs
  
  
# Calculate R2 and p value
depth_vcontig_lm_fit <- lm(
  log10(mean_sample_vlen) ~ log10(seq_depth),
  data = summary_depth_vcontig_df
)
depth_vcontig_r2 <- summary(depth_vcontig_lm_fit)$r.squared #R2 value
depth_vcontig_pval <- summary(depth_vcontig_lm_fit)$coefficients["log10(seq_depth)", "Pr(>|t|)"] #p-value testing slope=0

# Plot linear regression
depth_vs_vlen <- ggplot(summary_depth_vcontig_df, aes(seq_depth, y=mean_sample_vlen, color=Experiment)) +
  geom_errorbar(data=summary_depth_vcontig_df,
                aes(x = seq_depth,
                    ymin = mean_sample_vlen - se_sample_vlen, ymax = mean_sample_vlen + se_sample_vlen),
                width = 0.1,
                linewidth = 0.8,
                position = position_dodge(width = 0.75),
                #inherit.aes = FALSE
  ) +
  geom_point(size = 3, shape = 18) +
  geom_smooth(aes(seq_depth, y=mean_sample_vlen),
              method="lm", formula=y~x, se=FALSE, 
              linewidth=1, linetype="dashed", color="black",
              inherit.aes=FALSE) +
  annotate("text", x=Inf, y=Inf, label=paste0("R² = ", round(depth_vcontig_r2, 3), "\np = ", signif(depth_vcontig_pval, digits=3)),
           hjust = 1.1,
           vjust = 1.5,
           size = 4) +
  scale_x_log10() +
  scale_y_log10() +
  labs(x="sample sequencing depth (number of total reads)", y="mean viral contig length of the sample (bp)") +
  theme_bw()


####################################### Save plots
# enable showtext
#library(showtext)
#font_add("Arial", regular = "/System/Library/Fonts/Supplemental/Arial.ttf")  # macOS path
#showtext_auto()  # enables showtext for all plots

# get the plots that are going to be saved
plots <- list(depth_vs_vlen)

# specify the destination folder
dest_folder <- "/Users/jduan/bushman/virome_methods/paper/manuscript/msystems_submission3/figures"
if(!dir.exists(dest_folder)) dir.create(dest_folder, recursive = TRUE)

# loop over plots and save as PDF
for (i in seq_along(plots)) {
  print(plots[[i]])   # ensures the plot is drawn
  ggsave(
    filename = file.path(dest_folder, paste0("Supp_FigS9_", i, ".pdf")),
    plot = plots[[i]],
    width = 6.5,
    height = 7.5,
    units = "in"
    # no need to specify device; showtext handles font embedding
  )
}


################## Export table results as excel
# Create a blank workbook
OUT <- createWorkbook()

# Add some sheets to the workbook
addWorksheet(OUT, "vcontig_depth")
addWorksheet(OUT, "summary_vcontig_depth")

# Write the data to the sheets
writeData(OUT, sheet = "vcontig_depth", x = depth_vcontig_df)
writeData(OUT, sheet = "summary_vcontig_depth", x = summary_depth_vcontig_df)


# Export the file
setwd("/Users/jduan/bushman/virome_methods/paper/manuscript/msystems_submission3/viral_contig_lens")
saveWorkbook(OUT, "depth_viral_contiglen_v3.xlsx")

