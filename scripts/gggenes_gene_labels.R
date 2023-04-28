# Plotting fixed mutations at Appelmans endpoint with gggenes
# Labels only for gene product of mutated genes

setwd("/Volumes/WhitesonBigMac/adamsed/projects_main/Appelmans/Appelmans2022/breseq/R_files")

library(tidyverse)
library(ggplot2)
library(gggenes)
library(ggrepel)

# Define Phage Genomes (annotation tables exported from Geneious)
Bob_genome <- read_tsv("Bob_table.tsv") |> 
  rename(
    "gene"="Name",
    "start"="Minimum",
    "end"="Maximum",
    "strand"="Direction",
    "type"="Track Name"
  ) |>
  add_column(molecule = "Bob")
Bob_genome$strand[Bob_genome$strand == "forward"] <- 1
Bob_genome$strand[Bob_genome$strand == "reverse"] <- 0
Bob_genome$strand <- as.numeric(Bob_genome$strand)

Car_genome <- read_tsv("Car_table.tsv") |> 
  rename(
    "gene"="Name",
    "start"="Minimum",
    "end"="Maximum",
    "strand"="Direction",
    "type"="Track Name"
  ) |>
  add_column(molecule = "Car")
Car_genome$strand[Car_genome$strand == "forward"] <- 1
Car_genome$strand[Car_genome$strand == "reverse"] <- 0
Car_genome$strand <- as.numeric(Car_genome$strand)

CCS4_genome <- read_tsv("CCS4_table.tsv") |> 
  rename(
    "gene"="Name",
    "start"="Minimum",
    "end"="Maximum",
    "strand"="Direction",
    "type"="Track Name"
  ) |>
  add_column(molecule = "CCS4")
CCS4_genome$strand[CCS4_genome$strand == "forward"] <- 1
CCS4_genome$strand[CCS4_genome$strand == "reverse"] <- 0
CCS4_genome$strand <- as.numeric(CCS4_genome$strand)

# read genome length tags file
lengths <- read_csv("genome_lengths.csv")
lengths$name <- sub("^", "           ", lengths$name)

# bind genomes into sister and cousin phage files, transform gene names for plot labels
sisters <- rbind(Bob_genome, Car_genome)
sisters$gene <- gsub(", ",",\n", sisters$gene)
sisters$gene <- gsub("(.{10,}?)\\s", "\\1\n", sisters$gene)
sisters$gene <- gsub("\nI"," I", sisters$gene)

cousins <- rbind(Bob_genome, CCS4_genome)
cousins$gene <- gsub(", ",",\n", cousins$gene)
cousins$gene <- gsub("(.{10,}?)\\s", "\\1\n", cousins$gene)
cousins$gene <- gsub("\nI"," I", cousins$gene)

#read in mutation data, clean and arrange data, generate list of mutated genes
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

# sis_t0 data
s0 <- read_tsv("sister_t0.tsv")
s0$title[s0$title=="Phage_Sisters_rep1"] <- "S1M_T00"
s0$title[s0$title=="Phage_Sisters_rep2"] <- "S2M_T00"
s0$title[s0$title=="Phage_Sisters_rep3"] <- "S3M_T00"

# sis_mix data
sm <- read_tsv("sister_M.tsv") |> 
  select(!"new_copy_number")

sm <- rbind(s0, sm)

sm_sub <- sm[ ,keep] |> 
  separate(col = title, into = c("replicate", "time"))

sm_sub$replicate <- parse_number(sm_sub$replicate)

sm_sub$time <- parse_number(sm_sub$time)

sm_sub$seq_id <- gsub("Enterococcus_phage_", "", as.character(sm_sub$seq_id))

sis_mix <- sm_sub  |>
  filter(seq_id != "Bill") |>
  group_by(seq_id, time, position, mutation_category, gene_product) |>
  summarise(freq=mean(frequency))

# sis_ind  data
si <- read_tsv("sister_I.tsv") |>
  filter(mutation_category!="snp_synonymous|nonsynonymous")

si <- rbind(s0, si)

si_sub <- si[ ,keep] |>
  separate(col = title, into = c("replicate", "time"))

si_sub$replicate <- parse_number(si_sub$replicate)

si_sub$time <- parse_number(si_sub$time)

si_sub$seq_id <- gsub("Enterococcus_phage_", "", as.character(si_sub$seq_id))

sis_ind <- si_sub |>
  filter(seq_id != "Bill") |>
  group_by(seq_id, time, position, mutation_category, gene_product) |>
  summarise(freq=mean(frequency))

# cus_t0 data
c0 <- read_tsv("cousin_t0.tsv")
c0$title[c0$title=="Phage_Cousins_rep1"] <- "C1M_T00"
c0$title[c0$title=="Phage_Cousins_rep2"] <- "C2M_T00"
c0$title[c0$title=="Phage_Cousins_rep3"] <- "C3M_T00"

# cus_mix data
cm <- read_tsv("cousin_M.tsv") |>
  subset(select = -size) |> 
  filter(mutation_category != "large_deletion")

cm <- rbind(c0, cm)

cm_sub <- cm[ ,keep] |>
  separate(col = title, into = c("replicate", "time"))

cm_sub$replicate <- parse_number(cm_sub$replicate)

cm_sub$time <- parse_number(cm_sub$time)

cm_sub$seq_id <- gsub("Enterococcus_phage_", "", as.character(cm_sub$seq_id))

cus_mix <- cm_sub |>
  filter(seq_id != "SDS1") |>
  group_by(seq_id, time, position, mutation_category, gene_product) |>
  summarise(freq=mean(frequency))

# cus_ind data
ci <- read_tsv("cousin_I.tsv") |>
  subset(select = -size) |> 
  filter(mutation_category != "large_deletion")

ci <- rbind(c0, ci)

ci_sub <- ci[ ,keep] |>
  separate(col = title, into = c("replicate", "time"))

ci_sub$replicate <- parse_number(ci_sub$replicate)

ci_sub$time <- parse_number(ci_sub$time)

ci_sub$seq_id <- gsub("Enterococcus_phage_", "", as.character(ci_sub$seq_id))

cus_ind <- ci_sub |>
  filter(seq_id != "SDS1") |>
  group_by(seq_id, time, position, mutation_category, gene_product) |>
  summarise(freq=mean(frequency))

# get logical lists where genes have mutations, subset genome files for mutated non-hypothetical genes
## sisters, individual hosts
sis_ind_mut = list()

for (i in 1:nrow(sis_ind)) {
  tmp = sis_ind$position[i] > sisters$start & 
    sis_ind$position[i] < sisters$end & 
    sis_ind$seq_id[i] == sisters$molecule
  sis_ind_mut[[length(sis_ind_mut)+1]] = tmp
}

sis_ind_logic <- Reduce("|",sis_ind_mut)

sis_ind_genes <- cbind(sisters, sis_ind_logic) |> 
  filter(sis_ind_logic == T)

sis_ind_mut = list()

for (i in 1:nrow(sis_ind)) {
  tmp = sis_ind$position[i] > sisters$start & 
    sis_ind$position[i] < sisters$end & 
    sis_ind$seq_id[i] == sisters$molecule & 
    sis_ind$freq[i] >= 0.1
  sis_ind_mut[[length(sis_ind_mut)+1]] = tmp
}

sis_ind_logic <- Reduce("|",sis_ind_mut)

sis_ind_names <- cbind(sisters, sis_ind_logic) |> 
  filter(sis_ind_logic == T) |>
  filter(type == "Annotated")

## sisters, mixed hosts
sis_mix_mut = list()

for (i in 1:nrow(sis_mix)) {
  tmp = sis_mix$position[i] > sisters$start & 
    sis_mix$position[i] < sisters$end & 
    sis_mix$seq_id[i] == sisters$molecule
  sis_mix_mut[[length(sis_mix_mut)+1]] = tmp
}

sis_mix_logic <- Reduce("|", sis_mix_mut)

sis_mix_genes <- cbind(sisters, sis_mix_logic) |> 
  filter(sis_mix_logic == T)

sis_mix_mut = list()

for (i in 1:nrow(sis_mix)) {
  tmp = sis_mix$position[i] > sisters$start & 
    sis_mix$position[i] < sisters$end & 
    sis_mix$seq_id[i] == sisters$molecule & 
    sis_mix$freq[i] >= 0.1
  sis_mix_mut[[length(sis_mix_mut)+1]] = tmp
}

sis_mix_logic <- Reduce("|", sis_mix_mut)

sis_mix_names <- cbind(sisters, sis_mix_logic) |> 
  filter(sis_mix_logic == T) |>
  filter(type == "Annotated")

# cousins, individual hosts
cus_ind_mut = list()

for (i in 1:nrow(cus_ind)) {
  tmp = cus_ind$position[i] > cousins$start & 
    cus_ind$position[i] < cousins$end & 
    cus_ind$seq_id[i] == cousins$molecule
  cus_ind_mut[[length(cus_ind_mut)+1]] = tmp
}

cus_ind_logic <- Reduce("|",cus_ind_mut)

cus_ind_genes <- cbind(cousins, cus_ind_logic) |> 
  filter(cus_ind_logic == T)

cus_ind_mut = list()

for (i in 1:nrow(cus_ind)) {
  tmp = cus_ind$position[i] > cousins$start & 
    cus_ind$position[i] < cousins$end & 
    cus_ind$seq_id[i] == cousins$molecule & 
    cus_ind$freq[i] >= 0.1
  cus_ind_mut[[length(cus_ind_mut)+1]] = tmp
}

cus_ind_logic <- Reduce("|",cus_ind_mut)

cus_ind_names <- cbind(cousins, cus_ind_logic) |> 
  filter(cus_ind_logic == T) |>
  filter(type == "Annotated")

## cousins, mixed hosts
cus_mix_mut = list()

for (i in 1:nrow(cus_mix)) {
  tmp = cus_mix$position[i] > cousins$start & 
    cus_mix$position[i] < cousins$end & 
    cus_mix$seq_id[i] == cousins$molecule
  cus_mix_mut[[length(cus_mix_mut)+1]] = tmp
}

cus_mix_logic <- Reduce("|",cus_mix_mut)

cus_mix_genes <- cbind(cousins, cus_mix_logic) |> 
  filter(cus_mix_logic == T)

cus_mix_mut = list()

for (i in 1:nrow(cus_mix)) {
  tmp = cus_mix$position[i] > cousins$start & 
    cus_mix$position[i] < cousins$end & 
    cus_mix$seq_id[i] == cousins$molecule & 
    cus_mix$freq[i] >= 0.1
  cus_mix_mut[[length(cus_mix_mut)+1]] = tmp
}

cus_mix_logic <- Reduce("|",cus_mix_mut)

cus_mix_names <- cbind(cousins, cus_mix_logic) |> 
  filter(cus_mix_logic == T) |>
  filter(type == "Annotated")

# remove duplicate gene names
## sis_ind
sis_ind_bob <- sis_ind_names |> 
  filter(molecule == "Bob")
sis_ind_bob <- sis_ind_bob[match(unique(sis_ind_bob$gene), sis_ind_bob$gene),]

sis_ind_car <- sis_ind_names |> 
  filter(molecule == "Car")
sis_ind_car <- sis_ind_car[match(unique(sis_ind_car$gene), sis_ind_car$gene),]

sis_ind_gene_uniq <- rbind(sis_ind_bob, sis_ind_car)

## sis_mix
sis_mix_bob <- sis_mix_names |> 
  filter(molecule == "Bob")
sis_mix_bob <- sis_mix_bob[match(unique(sis_mix_bob$gene), sis_mix_bob$gene),]

sis_mix_car <- sis_mix_names |> 
  filter(molecule == "Car")
sis_mix_car <- sis_mix_car[match(unique(sis_mix_car$gene), sis_mix_car$gene),]

sis_mix_gene_uniq <- rbind(sis_mix_bob, sis_mix_car)

# plot genomes with labels for mutated non-hypothetical genes
## Sisters, individual hosts
si_plot <- ggplot(
  sisters,
  aes(
    xmin = start,
    xmax = end,
    y = molecule,
    fill = gene,
    forward = strand
  )
) +
  geom_gene_arrow(
    inherit.aes=FALSE,
    aes(
      xmin = start,
      xmax = end,
      y = molecule,
      forward = strand
    )
  ) +
  geom_gene_arrow(
    data = sis_ind_genes,
    inherit.aes=FALSE,
    aes(
      xmin = start,
      xmax = end,
      y = molecule,
      fill = type,
      forward = strand
    )
  ) +
  geom_feature(
    data = lengths[lengths$molecule != "CCS4", ],
    feature_height = unit(-10, "mm"),
    aes(x = position, y = molecule)
  ) +
  geom_feature_label(
    data = lengths[lengths$molecule != "CCS4", ], 
    feature_height = unit(-5, "mm"),
    label_height = unit(4, "mm"),
    aes(x = position, y = molecule, label = name),
    size = 15
  ) +
  geom_feature(
    data = sis_ind |> rename("molecule" = "seq_id"),
    feature_height = unit(-4, "mm"),
    aes(x = position, y = molecule)
  ) +
  facet_wrap(~ molecule, scales = "free", ncol = 1) +
  scale_x_discrete(expand=c(0,500)) +
  scale_y_discrete(expand=c(0,0)) +
  ylab(NULL) +
  xlab(NULL) +
  theme_genes() +
  geom_text_repel(
    data = sis_ind_gene_uniq |> mutate(start = (start + end)/2),
    aes(
      x = start,
      y = molecule,
      label = gene,
      segment.alpha = .4),
    size = 4.5,
    inherit.aes = F,
    max.overlaps = Inf,
    box.padding = 1,
    ylim = c(2,NA),
    xlim = c(NA, NA)
  ) +
  scale_fill_manual(values=c('blue','yellow','red')) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 20)
  )
ggsave("si_labels.png", plot = si_plot, width = 3900, height = 1700, units = "px", dpi = 300)

## Sisters, mixed hosts
sm_plot <- ggplot(
  sisters,
  aes(
    xmin = start,
    xmax = end,
    y = molecule,
    fill = gene,
    forward = strand
  )
) +
  geom_gene_arrow(
    inherit.aes=FALSE,
    aes(
      xmin = start,
      xmax = end,
      y = molecule,
      forward = strand
    )
  ) +
  geom_gene_arrow(
    data = sis_mix_genes,
    inherit.aes=FALSE,
    aes(
      xmin = start,
      xmax = end,
      y = molecule,
      fill = type,
      forward = strand
    )
  ) +
  geom_feature(
    data = lengths[lengths$molecule != "CCS4", ],
    feature_height = unit(-10, "mm"),
    aes(x = position, y = molecule)
  ) +
  geom_feature_label(
    data = lengths[lengths$molecule != "CCS4", ],
    feature_height = unit(-5, "mm"),
    label_height = unit(4, "mm"),
    aes(x = position, y = molecule, label = name),
    size = 15
  ) +
  geom_feature(
    data = sis_mix |> rename("molecule" = "seq_id"),
    feature_height = unit(-4, "mm"),
    aes(x = position, y = molecule)
  ) +
  facet_wrap(~ molecule, scales = "free", ncol = 1) +
  scale_x_discrete(expand=c(0,500)) +
  scale_y_discrete(expand=c(0,0)) +
  ylab(NULL) +
  xlab(NULL) +
  theme_genes() +
  geom_text_repel(
    data = sis_mix_gene_uniq |> mutate(start = (start + end)/2),
    aes(x=start,y=molecule,label=gene, segment.alpha=.4),
    size = 4.5,
    inherit.aes = F,
    max.overlaps = Inf,
    box.padding = .75,
    ylim = c(3,NA),
    xlim = c(NA, NA)
  ) +
  scale_fill_manual(values=c('blue','yellow','red')) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 20)
  )
ggsave("sm_labels.png", plot = sm_plot, width = 3900, height = 1900, units = "px", dpi = 300)

## Cousins, individual hosts
ci_plot <- ggplot(
  cousins,
  aes(
    xmin = start,
    xmax = end,
    y = molecule,
    fill = gene,
    forward = strand
  )
) +
  geom_gene_arrow(
    inherit.aes=FALSE,
    aes(
      xmin = start,
      xmax = end,
      y = molecule,
      forward = strand
    )
  ) +
  geom_gene_arrow(
    data = cus_ind_genes,
    inherit.aes=FALSE,
    aes(
      xmin = start,
      xmax = end,
      y = molecule,
      fill = type,
      forward = strand
    )
  ) +
  geom_feature(
    data = lengths[lengths$molecule != "Car", ],
    feature_height = unit(-10, "mm"),
    aes(x = position, y = molecule)
  ) +
  geom_feature_label(
    data = lengths[lengths$molecule != "Car", ],
    feature_height = unit(-5, "mm"),
    label_height = unit(4, "mm"),
    aes(x = position, y = molecule, label = name),
    size = 15
  ) +
  geom_feature(
    data = cus_ind |> rename("molecule" = "seq_id"),
    feature_height = unit(-4, "mm"),
    aes(x = position, y = molecule)
  ) +
  facet_wrap(~ molecule, scales = "free", ncol = 1) +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  ylab(NULL) +
  xlab(NULL) +
  theme_genes() +
  geom_text_repel(
    data = cus_ind_names |> mutate(start = (start + end)/2),
    aes(x=start,y=molecule,label=gene, segment.alpha=.4),
    size = 4.5,
    inherit.aes = F,
    box.padding = 1,
    ylim = c(2,NA),
    xlim = c(NA, NA)
  ) +
  scale_fill_manual(values=c('blue','yellow','red')) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 20)
  )
ggsave("ci_labels.png", plot = ci_plot, width = 3900, height = 1600, units = "px", dpi = 300)

## Cousins, mixed hosts
cm_plot <- ggplot(
  cousins,
  aes(
    xmin = start,
    xmax = end,
    y = molecule,
    fill = gene,
    forward = strand
  )
) +
  geom_gene_arrow(
    inherit.aes=FALSE,
    aes(
      xmin = start,
      xmax = end,
      y = molecule,
      forward = strand
    )
  ) +
  geom_gene_arrow(
    data = cus_mix_genes,
    inherit.aes=FALSE,
    aes(
      xmin = start,
      xmax = end,
      y = molecule,
      fill = type,
      forward = strand
    )
  ) +
  geom_feature(
    data = lengths[lengths$molecule != "Car", ],
    feature_height = unit(-10, "mm"),
    aes(x = position, y = molecule)
  ) +
  geom_feature_label(
    data = lengths[lengths$molecule != "Car", ],
    feature_height = unit(-5, "mm"),
    label_height = unit(4, "mm"),
    aes(x = position, y = molecule, label = name),
    size = 15
  ) +
  geom_feature(
    data = cus_mix |> rename("molecule" = "seq_id"),
    feature_height = unit(-4, "mm"),
    aes(x = position, y = molecule)
  ) +
  facet_wrap(~ molecule, scales = "free", ncol = 1) +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  ylab(NULL) +
  xlab(NULL) +
  theme_genes() +
  geom_text_repel(
    data = cus_mix_names |> mutate(start = (start + end)/2),
    aes(x=start,y=molecule,label=gene, segment.alpha=.4),
    size = 4.5,
    inherit.aes = F,
    box.padding = 1,
    ylim = c(2,NA),
    xlim = c(NA, NA)
  ) +
  scale_fill_manual(values=c('blue','yellow','red')) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 20)
  )
ggsave("cm_labels.png", plot = cm_plot, width = 3900, height = 1600, units = "px", dpi = 300)


## histograms of mutation frequency per 10 Kb
ggplot(unique(sis_ind[,c("seq_id","position")]), aes(x = position)) +
  geom_histogram(breaks = seq(0,150000,10000)) +
  stat_bin(breaks = seq(0,150000,10000), geom="text", aes(label=after_stat(count)), vjust = -0.1) +
  facet_grid(seq_id~.)

ggplot(unique(sis_mix[,c("seq_id","position")]), aes(x = position)) +
  geom_histogram(breaks = seq(0,150000,10000)) +
  stat_bin(breaks = seq(0,150000,10000), geom="text", aes(label=after_stat(count)), vjust = -0.1) +
  facet_grid(seq_id~.)

ggplot(unique(cus_ind[,c("seq_id","position")]), aes(x = position)) +
  geom_histogram(breaks = seq(0,150000,10000)) +
  stat_bin(breaks = seq(0,150000,10000), geom="text", aes(label=after_stat(count)), vjust = -0.1) +
  facet_grid(seq_id~.)

ggplot(unique(cus_mix[,c("seq_id","position")]), aes(x = position)) +
  geom_histogram(breaks = seq(0,150000,10000)) +
  stat_bin(breaks = seq(0,150000,10000), geom="text", aes(label=after_stat(count)), vjust = -0.1) +
  facet_grid(seq_id~.)

### end of script ###