# Figure S3. Comparison of methods for purifying viral particles from stool.

# This script combines the reference-based and contig-based results of 250505_nextseq and 250617_nextseq and plots using facet by conditions
# and summarizes the VirMock1 distribution result using PCoA
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
library(tibble)

# Get the metadata
setwd("/Users/jduan/bushman/virome_methods/250429_nextseq")
metadata_250505 <- read_excel("metadata.xlsx")
setwd("/Users/jduan/bushman/virome_methods/250616_nextseq")
metadata_250617 <- read_excel("metadata.xlsx")
metadata <- bind_rows(metadata_250505,metadata_250617)

# Get the output excel file for 250505_nextseq
setwd("/Users/jduan/bushman/virome_methods/250429_nextseq/R analysis")
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
setwd("/Users/jduan/bushman/virome_methods/250616_nextseq/R analysis")
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
cenote_reads_df <- cenote_reads_df %>% mutate(contig_type="viral")
ref_reads_df <- ref_reads_df %>% mutate(contig_type="ref viral") %>% rename(percent_viral_reads=percent_ref_reads, sum_viral_reads=sum_viralref_reads)
viral_reads_abundance <- rbind(cenote_reads_df, ref_reads_df)

# Create a bar plot for the ref_rpkm of each replicate (+ ref rpkm annotated on graph)
i_ref_abundance_plot <- ref_df %>% 
  select(Experiment, Sample_ID, Sample_type, Mock_Community, Treatment, Nuclease, Amplification, VC_aliquot_used, reference_genome, ref_rpkm, mapped_ref_reads, total_reads) %>% 
  left_join(ref_reads_df %>% 
              ungroup() %>% 
              select (Experiment, Sample_ID, sum_replicate_rpkm, percent_viral_reads),
            by = c("Sample_ID" = "Sample_ID", "Experiment"="Experiment")) #match the sum rpkm to their respective Sample_ID so that each input reference in each group has a matching total rpkm
i_ref_abundance_plot <- i_ref_abundance_plot %>% mutate(ref_relative_abundance=ref_rpkm/sum_replicate_rpkm) #calculate the relative abundance of each input reference


############# Ecllipse PCoA plot for Treatment method using Bray-Curtis

stool_ref_abundance_plot <- i_ref_abundance_plot %>% filter(Experiment=="stool spike-in", Mock_Community=="mock community spiked-in")

# Construct a sample-by-species table (row=sample, column=species) with relative abundance
stool_ref_abun_mat <- stool_ref_abundance_plot %>% select(Sample_ID, reference_genome, ref_relative_abundance) %>% 
  pivot_wider(names_from=reference_genome,values_from=ref_relative_abundance)
stool_ref_abun_mat <- data.frame(stool_ref_abun_mat) #convert to dataframe
rownames(stool_ref_abun_mat) <- stool_ref_abun_mat$Sample_ID #make sample names as row names
stool_ref_abun_mat <- stool_ref_abun_mat %>% select(-Sample_ID) #get rid of the sample column
stool_ref_abun_mat[is.na(stool_ref_abun_mat)] <- 0 #replace all NA with 0 before computing bray-curtis
stool_ref_abun_mat <- stool_ref_abun_mat[rowSums(stool_ref_abun_mat) > 0, ]

# Compute Bray-Curtis distance to represent differences in abundance patterns
stool_ref_bray <- vegdist(stool_ref_abun_mat, method="bray")

# Run Prinicipal Coordinate Analysis
stool_ref_pcoa <- cmdscale(stool_ref_bray, k=2, eig=TRUE) # k=2 means 2 dimensions (axes); eig=TRUE returns eigenvalues (indicating how much variance)

# Convert to dataframe and join with the metadata
stool_metadata <- stool_ref_abundance_plot %>%
  distinct(Sample_ID, Treatment, Sample_type, Mock_Community)
#stool_metadata <- metadata %>% filter(Experiment=="stool spike-in", Sample_type=="stool")

stool_ref_pcoa_df <- data.frame(sample=rownames(stool_ref_pcoa$points),
                                Axis1=stool_ref_pcoa$points[,1],
                                Axis2=stool_ref_pcoa$points[,2]
) %>%
  left_join(stool_metadata, by=c("sample"="Sample_ID"))

# Run PERMANOVA by-factor test
stool_metadata <- stool_metadata %>% column_to_rownames("Sample_ID")
mod1 <- adonis2(stool_ref_bray ~ Treatment, data = stool_metadata, by = "margin")
summary_adonis_mod <- data.frame(
  Model=c("Treatment"),
  R2=c(mod1$R2[1]),
  P_value= c(mod1$`Pr(>F)`[1])
)
print(summary_adonis_mod, row.names=FALSE)

adonis_text <- paste0(summary_adonis_mod$Model, ": R² = ", round(summary_adonis_mod$R2, 3),
                      ", p = ", summary_adonis_mod$P_value, collapse = "\n")

# Make the ecllipse plot
ref_pcoa_ecllipse <- ggplot(stool_ref_pcoa_df, aes(x=Axis1, y=Axis2, color=Treatment, shape=Sample_type)) +
  geom_point(size=3, alpha=0.8) +
  stat_ellipse(aes(group = Treatment), level = 0.95, linetype = "dashed") + # draw 95% confidence ellipses for each amplification method, dash lines show approximate spread of replicates within each group
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )+
  coord_equal() +
  labs(
    x = paste0("PCoA1 (", round(stool_ref_pcoa$eig[1] / sum(stool_ref_pcoa$eig) * 100, 1), "%)"),
    y = paste0("PCoA2 (", round(stool_ref_pcoa$eig[2] / sum(stool_ref_pcoa$eig) * 100, 1), "%)"),
    title = "Stool spike-in experiment VirMock1 distribution \nPCoA (Bray-curtis)",
    color="Treatment",
    shape="Sample type"
  ) +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.1,
           label=adonis_text, size=2.5)


####################################### Save plots

# get the plots that are going to be saved
plots <- list(ref_pcoa_ecllipse)

# specify the destination folder
dest_folder <- "/Users/jduan/bushman/virome_methods/paper/manuscript/msystems_submission2/figure"
if(!dir.exists(dest_folder)) dir.create(dest_folder, recursive = TRUE)

# loop over plots and save as PDF
for (i in seq_along(plots)) {
  print(plots[[i]])   # ensures the plot is drawn
  ggsave(
    filename = file.path(dest_folder, paste0("Resub_FigS4_", i, ".pdf")),
    plot = plots[[i]],
    width = 6.5,
    height = 7.5,
    units = "in"
    # no need to specify device; showtext handles font embedding
  )
}

