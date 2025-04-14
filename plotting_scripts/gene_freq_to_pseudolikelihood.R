library(openxlsx)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(tidyverse)

# read in gene frequencies and generate dictionary
absolute_gene_IDs <- TRUE
filename <- "gene_freq/S_pneumoniae_gene_frequency_All_absolute_geneIDs.txt"
gene.freq.df <- read.csv(filename, header = 1, sep = "\t")
gene.freq.df$Gene_ID <- as.numeric(gene.freq.df$Gene_ID)

# run for training dataset
filename <- "log_pseudolikelihoods/stratified-shuffled/per-gene/AtB_All_S_pneumoniae_training_N50_stratified_contig_shuf_encoder_decoder_drop0_pseudoL_per_gene.txt"
data <- read.csv(filename, header = 1, sep = "\t")
if (absolute_gene_IDs == TRUE)
{
  data$Gene_ID <- abs(as.numeric(data$Gene_ID))
}

summary.df <- data %>% 
  group_by(Gene_ID) %>% 
  summarize(Avg_log_pseudolikelihood = mean(log_pseudolikelihood), Std_log_pseudolikelihood = sd(log_pseudolikelihood))
summary.df <- summary.df %>% replace(is.na(.), 0)
combined.df <- merge(summary.df, gene.freq.df, by = c("Gene_ID"))

# plot average
p <- ggplot(combined.df, aes(x = Genome_freq, y=Avg_log_pseudolikelihood,)) +
  geom_point(alpha=0.2) +
  geom_smooth(method="loess", colour="blue", fill = "red") +
  theme_light() +
  labs(x = "Gene Frequency", y = "Gene log pseudolikelihood")
p

ggsave(file="Pneumo_log_pseudolikelihood_vs_gene_freq_encoder_decoder_training_data_all_stratified_inference_loess.svg", plot=p, height = 6, width = 8)
ggsave(file="Pneumo_log_pseudolikelihood_vs_gene_freq_encoder_decoder_training_data_all_stratified_inference_loess.png", plot=p, height = 6, width = 8)

breaks <- seq(0, 1, by = 0.1)
p <- ggplot(combined.df, aes(x = cut(Genome_freq, breaks = breaks), y=Avg_log_pseudolikelihood,)) +
  geom_boxplot() +
  theme_light() +
  labs(x = "Gene Frequency", y = "Gene log pseudolikelihood")
p

ggsave(file="Pneumo_log_pseudolikelihood_vs_gene_freq_encoder_decoder_training_data_all_stratified_inference_boxplot.svg", plot=p, height = 6, width = 8)
ggsave(file="Pneumo_log_pseudolikelihood_vs_gene_freq_encoder_decoder_training_data_all_stratified_inference_boxplot.png", plot=p, height = 6, width = 8)

# plot everything
combined.df.all <- merge(data, gene.freq.df, by = c("Gene_ID"))
p <- ggplot(combined.df.all, aes(x = Genome_freq, y=log_pseudolikelihood,)) +
  geom_point(alpha=0.2) +
  theme_light() +
  labs(x = "Gene Frequency", y = "log pseudolikelihood")
p

# run for held-in dataset
filename <- "log_pseudolikelihoods/per_gene/S_pneumoniae_held_in_random_N50_encoder_decoder_drop0_pseudoL_per_gene.txt"
data <- read.csv(filename, header = 1, sep = "\t")
if (absolute_gene_IDs == TRUE)
{
  data$Gene_ID <- abs(as.numeric(data$Gene_ID))
}

summary.df <- data %>% 
  group_by(Gene_ID) %>% 
  summarize(Avg_log_pseudolikelihood = mean(log_pseudolikelihood), Std_log_pseudolikelihood = sd(log_pseudolikelihood))
summary.df <- summary.df %>% replace(is.na(.), 0)
combined.df <- merge(summary.df, gene.freq.df, by = c("Gene_ID"))

# plot average
p <- ggplot(combined.df, aes(x = Genome_freq, y=Avg_log_pseudolikelihood,)) +
  geom_point(alpha=0.2) +
  geom_smooth(method="loess", colour="blue", fill = "red") +
  theme_light() +
  labs(x = "Gene Frequency", y = "Gene log pseudolikelihood")
p

ggsave(file="Pneumo_log_pseudolikelihood_vs_gene_freq_encoder_decoder_novel_trained_data_all_inference_loess.svg", plot=p, height = 6, width = 8)
ggsave(file="Pneumo_log_pseudolikelihood_vs_gene_freq_encoder_decoder_novel_trained_data_all_inference_loess.png", plot=p, height = 6, width = 8)

breaks <- seq(0, 1, by = 0.1)
p <- ggplot(combined.df, aes(x = cut(Genome_freq, breaks = breaks), y=Avg_log_pseudolikelihood,)) +
  geom_boxplot() +
  theme_light() +
  labs(x = "Gene Frequency", y = "Gene log pseudolikelihood")
p

ggsave(file="Pneumo_log_pseudolikelihood_vs_gene_freq_encoder_decoder_novel_trained_data_all_inference_boxplot.svg", plot=p, height = 6, width = 8)
ggsave(file="Pneumo_log_pseudolikelihood_vs_gene_freq_encoder_decoder_novel_trained_data_all_inference_boxplot.png", plot=p, height = 6, width = 8)

# plot everything
combined.df.all <- merge(data, gene.freq.df, by = c("Gene_ID"))
p <- ggplot(combined.df.all, aes(x = Genome_freq, y=log_pseudolikelihood,)) +
  geom_point(alpha=0.2) +
  theme_light() +
  labs(x = "Gene Frequency", y = "log pseudolikelihood")
p

# run for held-out dataset
filename <- "log_pseudolikelihoods/per_gene/S_pneumoniae_held_out_cluster1_N50_encoder_decoder_drop0_pseudoL_per_gene.txt"
data <- read.csv(filename, header = 1, sep = "\t")
if (absolute_gene_IDs == TRUE)
{
  data$Gene_ID <- abs(as.numeric(data$Gene_ID))
}

summary.df <- data %>% 
  group_by(Gene_ID) %>% 
  summarize(Avg_log_pseudolikelihood = mean(log_pseudolikelihood), Std_log_pseudolikelihood = sd(log_pseudolikelihood))
summary.df <- summary.df %>% replace(is.na(.), 0)
combined.df <- merge(summary.df, gene.freq.df, by = c("Gene_ID"))

# plot average
p <- ggplot(combined.df, aes(x = Genome_freq, y=Avg_log_pseudolikelihood,)) +
  geom_point(alpha=0.2) +
  geom_smooth(method="loess", colour="blue", fill = "red") +
  theme_light() +
  labs(x = "Gene Frequency", y = "Gene log pseudolikelihood")
p

ggsave(file="Pneumo_log_pseudolikelihood_vs_gene_freq_encoder_decoder_novel_untrained_data_all_inference_loess.svg", plot=p, height = 6, width = 8)
ggsave(file="Pneumo_log_pseudolikelihood_vs_gene_freq_encoder_decoder_novel_untrained_data_all_inference_loess.png", plot=p, height = 6, width = 8)

breaks <- seq(0, 1, by = 0.1)
p <- ggplot(combined.df, aes(x = cut(Genome_freq, breaks = breaks), y=Avg_log_pseudolikelihood,)) +
  geom_boxplot() +
  theme_light() +
  labs(x = "Gene Frequency", y = "Gene log pseudolikelihood")
p

ggsave(file="Pneumo_log_pseudolikelihood_vs_gene_freq_encoder_decoder_novel_untrained_data_all_inference_boxplot.svg", plot=p, height = 6, width = 8)
ggsave(file="Pneumo_log_pseudolikelihood_vs_gene_freq_encoder_decoder_novel_untrained_data_all_inference_boxplot.png", plot=p, height = 6, width = 8)

# plot everything
combined.df.all <- merge(data, gene.freq.df, by = c("Gene_ID"))
p <- ggplot(combined.df.all, aes(x = Genome_freq, y=log_pseudolikelihood,)) +
  geom_point(alpha=0.2) +
  theme_light() +
  labs(x = "Gene Frequency", y = "log pseudolikelihood")
p

# run for random dataset
filename <- "log_pseudolikelihoods/per_gene/S_pneumoniae_held_in_random_N50_encoder_decoder_drop0_pseudoL_randomise_per_gene.txt"
data <- read.csv(filename, header = 1, sep = "\t")
if (absolute_gene_IDs == TRUE)
{
  data$Gene_ID <- abs(as.numeric(data$Gene_ID))
}

summary.df <- data %>% 
  group_by(Gene_ID) %>% 
  summarize(Avg_log_pseudolikelihood = mean(log_pseudolikelihood), Std_log_pseudolikelihood = sd(log_pseudolikelihood))
summary.df <- summary.df %>% replace(is.na(.), 0)
combined.df <- merge(summary.df, gene.freq.df, by = c("Gene_ID"))

# plot average
p <- ggplot(combined.df, aes(x = Genome_freq, y=Avg_log_pseudolikelihood,)) +
  geom_point(alpha=0.2) +
  geom_smooth(method="loess", colour="blue", fill = "red") +
  theme_light() +
  labs(x = "Gene Frequency", y = "Gene log pseudolikelihood")
p

ggsave(file="Pneumo_log_pseudolikelihood_vs_gene_freq_encoder_decoder_random_data_all_inference_loess.svg", plot=p, height = 6, width = 8)
ggsave(file="Pneumo_log_pseudolikelihood_vs_gene_freq_encoder_decoder_random_data_all_inference_loess.png", plot=p, height = 6, width = 8)

breaks <- seq(0, 1, by = 0.1)
p <- ggplot(combined.df, aes(x = cut(Genome_freq, breaks = breaks), y=Avg_log_pseudolikelihood,)) +
  geom_boxplot() +
  theme_light() +
  labs(x = "Gene Frequency", y = "Gene log pseudolikelihood")
p

ggsave(file="Pneumo_log_pseudolikelihood_vs_gene_freq_encoder_decoder_random_data_all_inference_boxplot.svg", plot=p, height = 6, width = 8)
ggsave(file="Pneumo_log_pseudolikelihood_vs_gene_freq_encoder_decoder_random_data_all_inference_boxplot.png", plot=p, height = 6, width = 8)

# plot everything
combined.df.all <- merge(data, gene.freq.df, by = c("Gene_ID"))
p <- ggplot(combined.df.all, aes(x = Genome_freq, y=log_pseudolikelihood,)) +
  geom_point(alpha=0.2) +
  theme_light() +
  labs(x = "Gene Frequency", y = "log pseudolikelihood")
p

