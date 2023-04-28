## Analysis of Breseq mutation data from Cousin_Individual Appelmans Experiment

# Set working directory
setwd("/Volumes/WhitesonBigMac/adamsed/projects_main/Appelmans/Appelmans2022/breseq/R_files")

# Load packages
library(dplyr)
library(ggplot2) 
library(tidyverse)
library(reshape2)
library(scales)
library(patchwork)
library(grid)

# Colorblind friendly palette
palette <- c(
  "snp_nonsynonymous" = "#88CCEE",
  "snp_intergenic" = "#CC6677",
  "snp_synonymous" = "#DDCC77",
  "snp_nonsense" = "#117733",
  "snp_noncoding" = "#332288",
  "small_indel" = "#AA4499",
  "large_amplification" = "#44AA99",
  "large_deletion" = "red"
)

# Define columns to keep
keep <- c(
  "frequency",
  "gene_product",
  "mutation_category",
  "position",
  "position_end",
  "position_start",
  "seq_id",
  "snp_type",
  "title"
)

# Import TSV files
si <- read_tsv("sister_I.tsv") |>
  filter(mutation_category!="snp_synonymous|nonsynonymous")
s0 <- read_tsv("sister_t0.tsv")
s0$title[s0$title=="Phage_Sisters_rep1"] <- "S1I_T00"
s0$title[s0$title=="Phage_Sisters_rep2"] <- "S2I_T00"
s0$title[s0$title=="Phage_Sisters_rep3"] <- "S3I_T00"

si <- rbind(s0, si)

# Subset necessary columns, split title into replicate and time, replace NA with blank string
si_sub <- si[ ,keep] |>
  separate(col = title, into = c("replicate", "time"))

# Clean up replicate column values
si_sub$replicate <- parse_number(si_sub$replicate)

# Clean up time point column values
si_sub$time <- parse_number(si_sub$time)

# Remove phage name prefix
si_sub$seq_id <- gsub("Enterococcus_phage_", "", as.character(si_sub$seq_id))

# calculate average frequencies from replicates
si_sum <- si_sub |>
  group_by(seq_id, time, position, mutation_category, gene_product) |>
  summarise(freq=mean(frequency))

# plot phage Bob heat map
si_bob_map <- ggplot(
  si_sum[si_sum$seq_id == "Bob", ],
  aes(x = factor(position), y = factor(time))
) +
  geom_tile(
    color="black",
    aes(fill=mutation_category, alpha = freq)
  ) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  theme_classic() +
  scale_fill_manual(values = palette) +
  xlab(NULL) +
  ylab("Cycle Number") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 2),
    axis.ticks.x = element_line(linewidth = 0.1),
    axis.title.y = element_text(size = 8),
    axis.text.y = element_text(size = 7)
  )
ggsave("si_bob_map.png", plot = si_bob_map, width = 6.5, height = 1, units = "in", dpi = 600)

# plot phage Car heat map
si_car_map <- ggplot(
  si_sum[si_sum$seq_id == "Car", ],
  aes(x = factor(position), y = factor(time))
) +
  geom_tile(
    color="black",
    aes(fill=mutation_category, alpha = freq)
  ) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  theme_classic() +
  scale_fill_manual(values = palette) +
  theme(legend.position = "none") +
  xlab(NULL) +
  ylab("Cycle Number") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 2),
    axis.ticks.x = element_line(linewidth = 0.1),
    axis.title.y = element_text(size = 8),
    axis.text.y = element_text(size = 7)
  )
ggsave("si_car_map.png", plot = si_car_map, width = 6.5, height = 1, units = "in", dpi = 600)

# plot all phages to generate complete legend
si_all_map <- ggplot(
  si_sum[si_sum$seq_id != "Bill", ],
  aes(x = factor(position), y = factor(time))
) +
  geom_tile(
    color="black",
    aes(fill=mutation_category, alpha = freq)
  ) +
  facet_grid(seq_id~.) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  theme_classic() +
  scale_fill_manual(
    name = "Mutation Type",
    labels = c(
      "Small Indel",
      "Intergenic SNP",
      "Noncoding SNP",
      "Nonsense SNP",
      "Nonsynonymous SNP",
      "Synonymous SNP"
    ),
    values = palette
  ) +
  xlab("Position") +
  ylab("Experiment Cylce") +
  scale_alpha_continuous(name="Frequency") +
  theme(legend.position = "bottom")

si_legend <- cowplot::get_legend(si_all_map)
si_legend <- ggpubr::as_ggplot(si_legend)
ggsave("si_legend.png", plot = si_legend, width = 6.5, height = 0.5, units = "in", dpi = 600, scale = 2)

### end of script ###