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
sm <- read_tsv("sister_M.tsv") |> 
  select(!"new_copy_number")
s0 <- read_tsv("sister_t0.tsv")
s0$title[s0$title=="Phage_Sisters_rep1"] <- "S1M_T00"
s0$title[s0$title=="Phage_Sisters_rep2"] <- "S2M_T00"
s0$title[s0$title=="Phage_Sisters_rep3"] <- "S3M_T00"

sm <- rbind(s0, sm)

# Subset necessary columns, split title into replicate and time, replace NA with blank string
sm_sub <- sm[ ,keep] |> 
  separate(col = title, into = c("replicate", "time"))

# Clean up replicate column values
sm_sub$replicate <- parse_number(sm_sub$replicate)

# Clean up time point column values
sm_sub$time <- parse_number(sm_sub$time)

# Remove phage name prefix
sm_sub$seq_id <- gsub("Enterococcus_phage_", "", as.character(sm_sub$seq_id))

# calculate average frequencies from replicates
sm_sum <- sm_sub  |>
  group_by(seq_id, time, position, mutation_category, gene_product) |>
  summarise(freq=mean(frequency))

# plot phage Bob heat map
sm_bob_map <- ggplot(
  sm_sum[sm_sum$seq_id == "Bob", ],
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
ggsave("sm_bob_map.png", plot = sm_bob_map, width = 6.5, height = 1, units = "in", dpi = 600)

# plot phage Car heat map
sm_car_map <- ggplot(
  sm_sum[sm_sum$seq_id == "Car", ],
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
    axis.text.x = element_text(size = 2),
    axis.ticks.x = element_line(linewidth = 0.1),
    axis.title.y = element_text(size = 8),
    axis.text.y = element_text(size = 7)
  )
ggsave("sm_car_map.png", plot = sm_car_map, width = 6.5, height = 1, units = "in", dpi = 600)

# plot all phages to generate complete legend
sm_all_map <- ggplot(
  sm_sum[sm_sum$seq_id != "Bill", ], 
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

sm_legend <- cowplot::get_legend(sm_all_map)
sm_legend <- ggpubr::as_ggplot(sm_legend)
ggsave("sm_legend.png", plot = sm_legend, width = 6.5, height = 0.5, units = "in", dpi = 600, scale = 2)

### end of script ###