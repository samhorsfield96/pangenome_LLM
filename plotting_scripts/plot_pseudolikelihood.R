library(openxlsx)
library(ggplot2)
library(ggsci)
library(ggpubr)

filename <- "pseudolikelihood_results_stratified.xlsx"

# data <- read.xlsx(filename, "random",)
# data$Type <- factor(data$Type, levels=c("Trained", "Similar", "Novel"))
# data$Inference <- factor(data$Inference, levels=c("Encoder", "Encoder-decoder"))
# data$Random <- factor(data$Random, levels=c("TRUE", "FALSE"))
# data$pseudolikelihood <- exp(data$log_pseudolikelihood)
# 
# data_subset <- subset(data, Cluster != "Random")
# # Plot two histograms, one for each group
# p <- ggplot(data_subset, aes(x = log_pseudolikelihood, fill = Type)) +
#   geom_histogram(position = "identity", alpha = 0.5, bins = 10) +
#   facet_grid(~Inference, scales="free_y") +
#   theme_minimal() +
#   labs(x = "log pseudolikelihood", y = "Count")
# p
# 
# pairwise_comparisons <- list(c("Trained", "Similar"), c("Trained", "Novel"), c("Similar", "Novel"))
# #pairwise_comparisons <- list(c("Similar", "Novel"))
# p <- ggplot(data_subset, aes(x = Type, y = log_pseudolikelihood, colour = Type)) +
#   geom_boxplot(outlier.shape = NA) + geom_jitter(size=2, alpha=0.9) +
#   facet_grid(Random~Inference, scales="free_y") +
#   guides(colour=guide_legend(title="Alignment")) +
#   theme_light() +
#   theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 11), legend.title=element_text(size=14,face="bold"), legend.text=element_text(size=12), legend.position = "none") +
#   scale_color_npg() + 
#   stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, vjust=-1, size = 5, aes( label=round(after_stat(y), digits=2))) + 
#   stat_summary(fun.y=mean, colour="black", geom="point", shape=18, size=2, show_guide = FALSE) + 
#   stat_compare_means(paired = FALSE, size = 5.5, label.x = 1.0, comparisons = pairwise_comparisons) +
#   labs(y = "log pseudolikelihood", x = "Prompted Genome")
# p
# 
# ggsave(file="Pneumo_log_pseudolikelihood_comparison_encoder_decoder_all_inference.svg", plot=p, height = 6, width = 10)
# ggsave(file="Pneumo_log_pseudolikelihood_comparison_encoder_decoder_all_inference.png", plot=p, height = 6, width = 10)

# pairwise_comparisons <- list(c("Trained", "Similar"), c("Trained", "Novel"), c("Similar", "Novel"))
# #pairwise_comparisons <- list(c("Similar", "Novel"))
# data_subset <- subset(data, Inference == "Encoder-decoder" & Cluster != "Random")
# p <- ggplot(data_subset, aes(x = Type, y = log_pseudolikelihood, colour = Type)) +
#   geom_boxplot(outlier.shape = NA) + geom_jitter(size=2, alpha=0.9) +
#   guides(colour=guide_legend(title="Alignment")) +
#   theme_light() +
#   theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 11), legend.title=element_text(size=14,face="bold"), legend.text=element_text(size=12), legend.position = "none") +
#   scale_color_npg() + 
#   stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, vjust=-1, size = 5, aes( label=round(after_stat(y), digits=2))) + 
#   stat_summary(fun.y=mean, colour="black", geom="point", shape=18, size=2, show_guide = FALSE) + 
#   stat_compare_means(paired = FALSE, size = 5.5, label.x = 1.0, comparisons = pairwise_comparisons) +
#   labs(y = "log pseudolikelihood", x = "Prompted Genome")
# p
# 
# ggsave(file="Pneumo_log_pseudolikelihood_comparison_encoder_decoder_only.svg", plot=p, height = 6, width = 8)
# ggsave(file="Pneumo_log_pseudolikelihood_comparison_encoder_decoder_only.png", plot=p, height = 6, width = 8)
# 
# data_subset <- data
# pairwise_comparisons <- list(c("TRUE", "FALSE"))
# p <- ggplot(data_subset, aes(x = Random, y = log_pseudolikelihood, colour = Type)) +
#   geom_boxplot(outlier.shape = NA) + geom_jitter(size=2, alpha=0.9) +
#   guides(colour=guide_legend(title="Alignment")) +
#   facet_grid(Type~Inference) +
#   theme_light() +
#   theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 11), strip.text.y = element_text(size = 11), legend.title=element_text(size=14,face="bold"), legend.text=element_text(size=12), legend.position = "none") +
#   scale_color_npg() + 
#   stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, vjust=-1, size = 5, aes( label=round(after_stat(y), digits=2))) + 
#   stat_summary(fun.y=mean, colour="black", geom="point", shape=18, size=2, show_guide = FALSE) + 
#   stat_compare_means(paired = FALSE, size = 5.5, vjust = -0.5, comparisons = pairwise_comparisons) +
#   scale_y_continuous(limits = c(NA, 10000)) +
#   labs(y = "log pseudolikelihood", x = "Randomised genome")
# p
# ggsave(file="Pneumo_log_pseudolikelihood_comparison_encoder_decoder_all_inference_random.svg", plot=p, height = 10, width = 10)
# ggsave(file="Pneumo_log_pseudolikelihood_comparison_encoder_decoder_all_inference_random.png", plot=p, height = 10, width = 10)

# plot trained genomes as random
data <- read.xlsx(filename, "All_Pneumo")
data$log_pseudolikelihood <- as.numeric(data$log_pseudolikelihood)
data$Type <- factor(data$Type, levels=c("Random", "Diff. Species", "Trained", "Similar", "Novel"))
data$Inference <- factor(data$Inference, levels=c("Encoder", "Encoder-decoder"))
data$Random <- factor(data$Random, levels=c("TRUE", "FALSE"))
data$pseudolikelihood <- exp(data$log_pseudolikelihood)

pairwise_comparisons <- list(c("Trained", "Similar"), c("Trained", "Novel"), c("Similar", "Novel"), 
                             c("Similar", "Random"), c("Novel", "Random"), c("Trained", "Random"))
#pairwise_comparisons <- list(c("Similar", "Novel"))
data_subset <- subset(data, Inference == "Encoder-decoder")
p <- ggplot(data_subset, aes(x = Type, y = log_pseudolikelihood, colour = Type)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(size=2, alpha=0.9) +
  guides(colour=guide_legend(title="Alignment")) +
  theme_light() +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 11), legend.title=element_text(size=14,face="bold"), legend.text=element_text(size=12), legend.position = "none") +
  scale_color_npg() + 
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, vjust=-1, size = 5, aes( label=round(after_stat(y), digits=2))) + 
  stat_summary(fun.y=mean, colour="black", geom="point", shape=18, size=2, show_guide = FALSE) + 
  stat_compare_means(paired = FALSE, size = 5.5, label.x = 1.0, comparisons = pairwise_comparisons) +
  labs(y = "log pseudolikelihood", x = "Prompted Genome")
p

ggsave(file="All_Pneumo_stratified_log_pseudolikelihood_comparison_encoder_decoder_all_plot.svg", plot=p, height = 6, width = 10)
ggsave(file="All_Pneumo_stratified_log_pseudolikelihood_comparison_encoder_decoder_all_plot.png", plot=p, height = 6, width = 10)

# for presentation
colour_vec <- c("Trained" = "#CD4309", "Similar" = "#2f47ed",  "Novel" = "#fcb43a", "Random" = "#0a895e", "Diff. Species" = "#8e6f48")
pairwise_comparisons <- list(c("Trained", "Similar"), c("Trained", "Novel"), c("Similar", "Novel"))
#pairwise_comparisons <- list(c("Similar", "Novel"))
data_subset <- subset(data, Inference == "Encoder-decoder" & Type != "Random" & Type != "Diff. Species")
p <- ggplot(data_subset, aes(x = Type, y = log_pseudolikelihood, colour = Type)) +
  geom_boxplot(outlier.shape = NA) + 
  guides(colour=guide_legend(title="Alignment")) +
  theme_light() +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 11), legend.title=element_text(size=14,face="bold"), legend.text=element_text(size=12), legend.position = "none") +
  scale_colour_manual(values = colour_vec) + 
  stat_compare_means(paired = FALSE, size = 5.5, label.x = 1.0, label = "p.signif", comparisons = pairwise_comparisons) +
  labs(y = "log pseudolikelihood", x = "Prompted Genome")
p
ggsave(file="All_Pneumo_stratified_log_pseudolikelihood_comparison_presentation_trained.png", plot=p, height = 6, width = 8)

pairwise_comparisons <- list(c("Random", "Similar"), 
                             c("Novel", "Random"), c("Novel", "Similar"))
#pairwise_comparisons <- list(c("Similar", "Novel"))
data_subset <- subset(data, Inference == "Encoder-decoder" & Type != "Trained" & Type != "Diff. Species")
p <- ggplot(data_subset, aes(x = Type, y = log_pseudolikelihood, colour = Type)) +
  geom_boxplot(outlier.shape = NA) + 
  guides(colour=guide_legend(title="Alignment")) +
  theme_light() +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 11), legend.title=element_text(size=14,face="bold"), legend.text=element_text(size=12), legend.position = "none") +
  scale_colour_manual(values = colour_vec) + 
  stat_compare_means(paired = FALSE, size = 5.5, label.x = 1.0, label = "p.signif", comparisons = pairwise_comparisons) +
  labs(y = "log pseudolikelihood", x = "Prompted Genome")
p
ggsave(file="All_Pneumo_stratified_log_pseudolikelihood_comparison_presentation_random.png", plot=p, height = 6, width = 8)

# all comparisons
pairwise_comparisons <- list(c("Diff. Species", "Random"), c("Random", "Similar"), 
                             c("Novel", "Random"))
#pairwise_comparisons <- list(c("Similar", "Novel"))
data_subset <- subset(data, Inference == "Encoder-decoder" & Type != "Trained")
p <- ggplot(data_subset, aes(x = Type, y = log_pseudolikelihood, colour = Type)) +
  geom_boxplot(outlier.shape = NA) + 
  guides(colour=guide_legend(title="Alignment")) +
  theme_light() +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 11), legend.title=element_text(size=14,face="bold"), legend.text=element_text(size=12), legend.position = "none") +
  scale_colour_manual(values = colour_vec) + 
  stat_compare_means(paired = FALSE, size = 5.5, label.x = 1.0, label = "p.signif", comparisons = pairwise_comparisons) +
  labs(y = "log pseudolikelihood", x = "Prompted Genome")
p
ggsave(file="All_Pneumo_log_pseudolikelihood_comparison_presentation_diff_species_random.png", plot=p, height = 6, width = 8)

