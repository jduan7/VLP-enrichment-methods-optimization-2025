#Figure 3. Comparison of the effects of different DNA amplification methods on viral genome recovery. 
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
library(tidyr)
library(vegan)

# Get the RPKM and cenote analysis tables
setwd("/Users/jduan/bushman/virome_methods/paper/codes/codes for publication/")
PTA_ref_viral_analysis<-read_excel("250221_miniseq_ref_R_analysis.xlsx", sheet="ref_viral_analysis")
PTA_ref_viral_percentage<-read_excel("250221_miniseq_ref_R_analysis.xlsx", sheet="ref_viral_percentage")
Others_ref_viral_analysis<-read_excel("250228_nextseq_v2_ref_R_analysis.xlsx", sheet="ref_viral_analysis")
Others_ref_viral_percentage<-read_excel("250228_nextseq_v2_ref_R_analysis.xlsx", sheet="ref_viral_percentage")

# Combine the tables
ref_viral_saliva <- rbind(PTA_ref_viral_analysis, Others_ref_viral_analysis)
ref_viral_percentage <- rbind(PTA_ref_viral_percentage, Others_ref_viral_percentage)
ref_viral_percentage <- ref_viral_percentage %>% mutate(contig_type="ref viral") %>% rename(percent_viral_reads=percent_ref_reads, sum_viral_reads=sum_viralref_reads)

# Calculate relative abundance and make bar plots
saliva_ref_abundance_plot <- ref_viral_saliva %>% select(sample, group, amplification, replicate, Treatment, Mock_Community, Sample_type, reference_genome, ref_rpkm, mapped_ref_reads, total_reads)
saliva_ref_abundance_plot <- saliva_ref_abundance_plot %>% 
  left_join(ref_viral_percentage %>% ungroup() %>% 
              select (amplification, replicate, sum_replicate_rpkm, percent_viral_reads),
            by = c("amplification" = "amplification", "replicate" = "replicate")) #match the sum rpkm to their respective group&type so that each input reference in each group has a matching total rpkm
saliva_ref_abundance_plot <- saliva_ref_abundance_plot %>% mutate(ref_relative_abundance=ref_rpkm/sum_replicate_rpkm) #calculate the relative abundance of each input reference

saliva_total_reads_data <- saliva_ref_abundance_plot %>% ungroup %>% group_by(amplification, replicate) %>% summarize(TotalReads=first(total_reads), .groups="drop")
saliva_total_reads_data$TotalReads <- formatC(saliva_total_reads_data$TotalReads, format="e", digits=1) #change to scientific annotation

# Filter for mock community spike in only
saliva_ref_abundance_plot <- saliva_ref_abundance_plot %>% filter(Mock_Community=="mock community spiked-in")

# Arrange in a certain order for later plotting
panel_order <- c("unamplified", "Genomiphi","PTA","WTA2_17cycles","MALBAC_17cycles")
saliva_ref_abundance_plot$amplification <- factor(saliva_ref_abundance_plot$amplification, levels = panel_order)
#viral_reads_abundance_saliva$amplification <- factor(viral_reads_abundance_saliva$amplification, levels = panel_order)

# Making relative abundance plots
p <- ggplot(saliva_ref_abundance_plot, aes(x=factor(replicate), y=ref_relative_abundance, fill=reference_genome)) + 
  geom_bar(stat="identity") + 
  scale_fill_paletteer_d("MoMAColors::Lupi") +
  scale_x_discrete(drop = T) +
  facet_nested(amplification ~ Mock_Community + Treatment + Sample_type, 
               scales = 'free_x', 
               nest_line = element_line(linetype=1)) +
  labs(x="Sample", y="relative abundance based on RPKM", 
       #title = "relative abundance based on reference genomes",
       fill = "reference_genome", color = "reference_genome") +
  theme_minimal(base_family="Helvetica") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background=element_blank(),
        ggh4x.facet.nestline=element_line(color="black"))

############# Ecllipse PCoA plot for amplification methods experiment using Bray-Curtis

# Construct a sample-by-species table (row=sample, column=species) with relative abundance
saliva_ref_abun_mat <- saliva_ref_abundance_plot %>% select(sample, reference_genome, ref_relative_abundance) %>% 
  pivot_wider(names_from=reference_genome,values_from=ref_relative_abundance)
saliva_ref_abun_mat <- data.frame(saliva_ref_abun_mat) #convert to dataframe
rownames(saliva_ref_abun_mat) <- saliva_ref_abun_mat$sample #make sample names as row names
saliva_ref_abun_mat <- saliva_ref_abun_mat %>% select(-sample) #get rid of the sample column

# Compute Bray-Curtis distance to represent differences in abundance patterns
saliva_ref_bray <- vegdist(saliva_ref_abun_mat, method="bray")

# Run Prinicipal Coordinate Analysis
saliva_ref_pcoa <- cmdscale(saliva_ref_bray, k=2, eig=TRUE) # k=2 means 2 dimensions (axes); eig=TRUE returns eigenvalues (indicating how much variance)

# Convert to dataframe and join with the metadata
saliva_metadata <- saliva_ref_abundance_plot %>%
  distinct(sample, amplification, group)
saliva_ref_pcoa_df <- data.frame(sample=rownames(saliva_ref_pcoa$points),
                                 Axis1=saliva_ref_pcoa$points[,1],
                                 Axis2=saliva_ref_pcoa$points[,2]
) %>%
  left_join(saliva_metadata, by="sample")

# Run PERMANOVA
saliva_adonis <- adonis2(saliva_ref_bray ~ amplification * group, data = saliva_metadata) #tests the combined effect of all terms (i.e. amplification, group, and their interaction)
print(saliva_adonis)

# Run by-factor test
mod1 <- adonis2(saliva_ref_bray ~ amplification, data = saliva_metadata, by = "margin")
mod2 <- adonis2(saliva_ref_bray ~ group, data = saliva_metadata, by = "margin")
mod3 <- adonis2(saliva_ref_bray ~ amplification+group, data = saliva_metadata, by = "margin")
mod4 <- adonis2(saliva_ref_bray ~ amplification*group, data = saliva_metadata, by = "margin")
summary_adonis_mod <- data.frame(
  Model=c("Amplification method"
          #"Sample group",
          #"Amplification method + Sample group",
          #"Amplification method x Sample group (interaction)"
  ),
  R2=c(mod1$R2[1]
       #mod2$R2[1],
       #sum(mod3$R2[1:2]),
       #sum(mod4$R2[1:3])
  ),
  P_value= c(mod1$`Pr(>F)`[1]
             #mod2$`Pr(>F)`[1],
             #max(mod3$`Pr(>F)`[1:2], na.rm = TRUE),
             #max(mod4$`Pr(>F)`[1:3], na.rm = TRUE)
  )
)
print(summary_adonis_mod, row.names=FALSE)

adonis_text <- paste0(summary_adonis_mod$Model, ": R² = ", round(summary_adonis_mod$R2, 3),
                      ", p = ", summary_adonis_mod$P_value, collapse = "\n")

# Make the ecllipse plot
pcoa_ecllipse <- ggplot(saliva_ref_pcoa_df, aes(x=Axis1, y=Axis2, color=amplification, shape=group)) +
  geom_point(size=3, alpha=0.8) +
  stat_ellipse(aes(group = amplification), level = 0.95, linetype = "dashed") + # draw 95% confidence ellipses for each amplification method, dash lines show approximate spread of replicates within each group
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )+
  coord_equal() +
  labs(
    x = paste0("PCoA1 (", round(saliva_ref_pcoa$eig[1] / sum(saliva_ref_pcoa$eig) * 100, 1), "%)"),
    y = paste0("PCoA2 (", round(saliva_ref_pcoa$eig[2] / sum(saliva_ref_pcoa$eig) * 100, 1), "%)"),
    #title = "PCoA (Bray-curtis)",
    color="Amplification method",
    shape="Sample group"
  ) +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.1,
           label=adonis_text, size=4)

