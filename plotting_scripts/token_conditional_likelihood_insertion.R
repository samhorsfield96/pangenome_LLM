library(ggplot2)
library(ggpubr)
library(ggsci)
library(tidyverse)
library(stringr)
library(ggbeeswarm)
library(dplyr)
library(optparse)

colours <- c("ST131 B/C (+ve)" = "#696969", 
             "ST131 B/C (-ve)" = "#E64B35FF",
             "ST131 A (+ve)" = "#4DBBD5FF",
             "ST131 A (-ve)" = "#4DBBD5FF",
             "ST95 (+ve)" = "#00A087FF",
             "ST95 (-ve)" = "#00A087FF",
             "ST73 (+ve)" = "#3C5488FF",
             "ST73 (-ve)" = "#3C5488FF")

parse_filename <- function(filename) {
  list(
    bla_pos = if (grepl("\\+ve", filename)) {
      TRUE
    } else {
      FALSE
    },
    ST131_BC = if (grepl("ST131B_C", filename)) {
      TRUE
    } else {
      FALSE
    },
    ST131_A = if (grepl("ST131A", filename)) {
      TRUE
    } else {
      FALSE
    },
    ST73 = if (grepl("CC73", filename)) {
      TRUE
    } else {
      FALSE
    },
    ST95 = if (grepl("CC95", filename)) {
      TRUE
    } else {
      FALSE
    },
    train = if (grepl("train", filename)) {
      TRUE
    } else {
      FALSE
    },
    test = if (grepl("test", filename)) {
      TRUE
    } else {
      FALSE
    },
    val = if (grepl("val", filename)) {
      TRUE
    } else {
      FALSE
    }
  )
}

iqr <- function(x, na.rm = TRUE) {
  qs <- quantile(x, probs = c(0.25, 0.75), na.rm = na.rm)
  (qs[2] - qs[1])
}

option_list <- list(
  make_option("--indir", type="character", default="token_conditonal_likelhoods/", help="Directory containing input files. Default = 'token_conditonal_likelhoods/'"),
  make_option("--outdir", type="character", default="figures/token_conditonal_likelhoods/", help="Output directory for saved plots. Default = 'figures/token_conditonal_likelhoods/'"),
  make_option("--outpref", type="character", default="output", help="Output filename prefix. Default = 'output'"),
  make_option("--labels", type="character", help="Comma-separated list of labels."),
  make_option("--plot-width", type="integer", default=100, help="Width in tokens around highest scoring location. Default = 100")
)
opt <- parse_args(OptionParser(option_list=option_list))

indir <- opt$indir
outdir <- opt$outdir
outpref <- opt$outpref
files_like <- Sys.glob(paste(indir,"*locus_pseudolikelihood.txt", sep = ""))
files_idx <- Sys.glob(paste(indir,"*_file_indices.txt", sep = ""))
files_position <- Sys.glob(paste(indir,"*_absolute_token_locations.csv", sep = ""))

df_idx <- data.frame(
  prefix = str_remove(basename(files_idx), "_N10_file_indices\\.txt$") |>
    str_remove("_rem$"),
  file_indices = files_idx,
  stringsAsFactors = FALSE
)

df_like <- data.frame(
  prefix = str_remove(basename(files_like), "_N10_locus_pseudolikelihood\\.txt$") |>
    str_remove("_rem$"),
  locus_pseudolikelihood = files_like,
  stringsAsFactors = FALSE
)

df_pos <- data.frame(
  prefix = str_remove(basename(files_position), "_absolute_token_locations\\.csv$") |>
    str_remove("_rem$"),
  locus_pseudolikelihood = files_position,
  stringsAsFactors = FALSE
)

merged <- merge(df_idx, df_like, by = "prefix")
merged <- merge(merged, df_pos, by = "prefix")

#final.df <- data.frame(Lineage = c(), Dataset = c(), GeneID = c(), GenomeID = c(), Quantile25 = c(), Quantile50 = c(), Quantile75 = c(), Min = c(), Max = c(), Mean = c())
# non pres
labels <- strsplit(opt$labels, ",")[[1]]
#pres
#labels <- c("ST131 B (-ve)", "ST131 C0 (-ve)", "ST131 C2 (-ve)", "ST131 C2 (+ve)", "ST95 (-ve)", "ST95 (+ve)")
#labels <- c("ST131 A (-ve)", "ST131 A (+ve)", "ST131 C2 (-ve)", "ST131 C2 (+ve)", "ST73 (-ve)")

# width around highest scoring location to plot
plot_width <- opt$`plot-width`

i <- 1
for (i in 1:nrow(merged))
{
  # TODO need to pair indices and locus pseudolikelihood files
  row <- merged[i,]
  
  prefix <- row[[1]]
  filename_idx <- row[[2]]
  filename_like <- row[[3]]
  filename_pos <- row[[4]]
  df.raw <- read.csv(filename_like, sep = "\t")
  df <- df.raw[, -c(1, 2)]
  df[df == 0] <- NA
  
  df.position <- read.csv(filename_pos, sep = ",", header = FALSE)
  colnames(df.position) <- c("Sample", "Target_token", "Position")
  
  sample_ids <- read.csv(filename_idx, sep = "\t", header = FALSE)
  colnames(sample_ids) <- "Sample"
  
  n_positions <- ncol(df)
  colnames(df) <- seq(1, n_positions)
  tmp.df <- cbind(sample_ids, df)
  
  tmp.df <- merge(tmp.df, df.position, by = "Sample", no.dups = TRUE)
  tmp.df <- tmp.df[!duplicated(tmp.df$Sample), ]
  
  parsed_values <- parse_filename(prefix)
  
  lineage <- NA
  lineage <- if (parsed_values$ST131_BC & is.na(lineage)) "ST131 B/C" else lineage
  lineage <- if (parsed_values$ST131_A & is.na(lineage)) "ST131 A" else lineage
  lineage <- if (parsed_values$ST73 & is.na(lineage)) "ST73" else lineage
  lineage <- if (parsed_values$ST95 & is.na(lineage)) "ST95" else lineage
  
  bla_pos <- parsed_values$bla_pos
  bla_pos_str <- if (parsed_values$bla_pos) "+ve" else "-ve"
  
  data_type <- NA
  data_type <- if (parsed_values$train & is.na(data_type)) "Training" else data_type
  data_type <- if (parsed_values$val & is.na(data_type)) "Validation" else data_type
  data_type <- if (parsed_values$test & is.na(data_type)) "Testing" else data_type
  
  max.val <- apply(df,1,max, na.rm= TRUE)
  min.val <- apply(df,1,min, na.rm= TRUE)
  
  max.col <- colnames(df)[apply(df,1,which.max)]
  min.col <- colnames(df)[apply(df,1,which.min)]
  
  # determine if max col matches true positions
  position.matching <- as.data.frame(cbind(as.numeric(max.col), as.numeric(tmp.df$Position)))
  position.matching$distance <- position.matching$V1 - position.matching$V2
  colnames(position.matching) <- c("True", "Predicted", "Difference")
  
  label <- paste0(lineage, " (", bla_pos_str, ")")
  
  to.append <- data.frame(Lineage = lineage, bla_pos = bla_pos, label = label, data_type = data_type, Max_value = max.val, Sample = sample_ids$Sample, True_position = position.matching$True, Prediction_position = position.matching$Predicted, Prediction_distance = position.matching$Difference)
  
  if (i == 1){
    final.df <- to.append
  } else
  {
    final.df <- rbind(final.df, to.append)
  }
}

# ordering
ordering <- c("ST131 B/C (+ve)", "ST131 B/C (-ve)", "ST131 A (+ve)", "ST131 A (-ve)", "ST95 (+ve)", "ST95 (-ve)", "ST73 (+ve)", "ST73 (-ve)")
final.df$label <- factor(final.df$label, levels = ordering)
final.df$data_type <- factor(final.df$data_type, levels = c("Training", "Testing"))

# # negative genomes
# pairwise_comparisons <- list(c("ST131 B/C (+ve)", "ST131 B/C (-ve)"), c("ST131 B/C (+ve)", "ST131 A (+ve)"),
#                              c("ST131 B/C (+ve)", "ST131 A (-ve)"), c("ST131 B/C (+ve)", "ST95 (+ve)"),
#                              c("ST131 B/C (+ve)", "ST95 (-ve)"), c("ST131 B/C (+ve)", "ST73 (+ve)"), c("ST131 B/C (+ve)", "ST73 (-ve)"))


summary.stats.df <- final.df %>%
  group_by(label, data_type) %>%
  summarise(
    IQR = iqr(Prediction_distance),
    .groups = "drop"
  )

summary.stats.df$IQR_text = paste0("IQR: ", round(summary.stats.df$IQR, digits=3))

p.all <- ggplot(final.df, aes(x = label, y = Prediction_distance, colour = label)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_beeswarm(size=2, alpha=0.7) +
  facet_grid(data_type~ .) +
  guides(colour=guide_legend(title="Alignment")) +
  geom_text(data = summary.stats.df, size = 5, colour = "black", aes(x = label, y = 5000, label = IQR_text)) +
  theme_light() +
  theme(axis.text.x = element_text(size = 18, angle = 45, hjust = 1), axis.text.y = element_text(size = 16), axis.title=element_text(size=20,face="bold"), strip.text.x = element_text(size = 16), strip.text.y = element_text(size = 16), 
        legend.title=element_text(size=14,face="bold"), legend.text=element_text(size=12), legend.position = "none") +
  scale_color_manual(values = colours) + 
  scale_y_continuous(breaks = function(z) seq(-5000, range(z)[2], by = 1000)) +
  #stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, vjust=-1, size = 5, aes( label=round(after_stat(y), digits=2))) + 
  #stat_summary(fun.y=mean, colour="black", geom="point", shape=18, size=2, show_guide = FALSE) + 
  #stat_compare_means(paired = FALSE, size = 5.5, label.x = 1.0, comparisons = pairwise_comparisons, label = "p.signif") +
  labs(y = "True-predicted location distance", x = "Lineage")
p.all

ggsave(file=paste0(outdir, outpref, "_max_value_position_comparison.png"), plot=p.all, height = 8, width = 7)
ggsave(file=paste0(outdir, outpref, "_max_value_position_comparison.svg"), plot=p.all, height = 8, width = 10)
