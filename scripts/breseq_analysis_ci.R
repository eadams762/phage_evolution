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
ci <- read_tsv("cousin_I.tsv") |>
  subset(select = -size) |> 
  filter(mutation_category != "large_deletion")
c0 <- read_tsv("cousin_t0.tsv")
c0$title[c0$title=="Phage_Cousins_rep1"] <- "C1I_T00"
c0$title[c0$title=="Phage_Cousins_rep2"] <- "C2I_T00"
c0$title[c0$title=="Phage_Cousins_rep3"] <- "C3I_T00"

ci <- rbind(c0, ci)

# Subset necessary columns, split title into replicate and time, replace NA with blank string
ci_sub <- ci[ ,keep] |>
  separate(col = title, into = c("replicate", "time"))

# Clean up replicate column values
ci_sub$replicate <- parse_number(ci_sub$replicate)

# Clean up time point column values
ci_sub$time <- parse_number(ci_sub$time)

# Remove phage name prefix
ci_sub$seq_id <- gsub("Enterococcus_phage_", "", as.character(ci_sub$seq_id))

# calculate average frequencies from replicates
ci_sum <- ci_sub |>
  group_by(seq_id, time, position, mutation_category, gene_product) |>
  summarise(freq=mean(frequency))

# plot phage Bob heat map
ci_bob_map <- ggplot(
  ci_sum[ci_sum$seq_id == "Bob", ],
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
ggsave("ci_bob_map.png", plot = ci_bob_map, width = 6.5, height = 1, units = "in", dpi = 600)

# plot phage CCS4 heat map
ci_ccs4_map <- ggplot(
  ci_sum[ci_sum$seq_id == "CCS4", ],
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
ggsave("ci_ccs4_map.png", plot = ci_ccs4_map, width = 6.5, height = 1, units = "in", dpi = 600)

# plot all phages to generate complete legend
ci_all_map <- ggplot(
  ci_sum[ci_sum$seq_id != "SDS1", ],
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
      "Nonsynonymous SNP",
      "Synonymous SNP"
    ),
    values = palette
  ) +
  xlab("Position") +
  ylab("Experiment Cylce") +
  scale_alpha_continuous(name="Frequency") +
  theme(legend.position = "bottom")

ci_legend <- cowplot::get_legend(ci_all_map)
ci_legend <- ggpubr::as_ggplot(ci_legend)
ggsave("ci_legend.png", plot = ci_legend, width = 6.5, height = 0.5, units = "in", dpi = 600, scale = 2)

### end of script ###