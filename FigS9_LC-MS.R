#Figure S3. Relative abundance of each modified and unmodified nucleoside in T4 C, T4 hmC, and T4 ghmC determined by LC-MS are shown in the bar graph. 
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

# Get the LC-MS data
setwd("/Users/jduan/bushman/virome_methods/paper/codes/codes for publication/")
MS<-read_excel("LC-MS_data.xlsx")
MS <- MS %>% select(-Peak_Area) %>% filter(nucleoside !="dT")

# Calculate the normalized signal and relative abundance
MS_summary <- MS %>% group_by(nucleoside) %>% 
  summarize(max_norm_signal=max(Normalized_by_dT_signal))
MS <- MS %>% left_join(MS_summary, by=c("nucleoside")) %>% mutate(relative_abundance=Normalized_by_dT_signal/max_norm_signal*100)

MS_plot <- ggplot(MS, aes(x=factor(nucleoside, levels=c("dC", "hm dC", "ghm dC")),
                          y=relative_abundance,
                          fill=factor(T4_strain, levels=c("T4C", "T4hmC", "T4ghmC")))
) +
  geom_bar(stat="identity", position = position_dodge(width = 0.9), alpha = 0.7) +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Set2") +
  labs(x = "nucleoside", y = "relative abundance (%)", fill = "T4 strain") +
  theme_minimal()
