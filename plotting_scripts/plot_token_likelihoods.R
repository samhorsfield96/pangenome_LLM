library(ggplot2)
library(ggpubr)
library(ggsci)

indir <- "token_likelihoods/"
files <- Sys.glob(paste(indir,"*locus_pseudolikelihood.txt", sep = ""))

final.df <- data.frame(Label = c(), Dataset = c(), GeneID = c(), GenomeID = c(), Quantile25 = c(), Quantile50 = c(), Quantile75 = c(), Min = c(), Max = c(), Mean = c())
#labels <- c("ST131 A (-ve)", "ST131 C2 (+ve)", "ST131 C2 (-ve)", "ST73 (-ve)")
labels <- c("ST131 (+ve)", "ST131 (-ve)", "ST73 (-ve)")

i <- 1
for (i in 1:length(files))
{
  df.raw <- read.csv(files[[i]], sep = "\t")
  df.rows <- df.raw[, c(1, 2)]
  df <- df.raw[, -c(1, 2)]
  df[df == 0] <- NA
  
  quantiles.val <- t(apply(df,1,quantile,c(0.25,0.50,0.75), na.rm= TRUE))
  avg.val <- apply(df,1,mean, na.rm= TRUE)
  max.val <- apply(df,1,max, na.rm= TRUE)
  min.val <- apply(df,1,min, na.rm= TRUE)
  
  to.append <- data.frame(Label = labels[[i]], Dataset = files[[i]], GeneID = df.rows$Gene_ID, GenomeID = df.rows$Genome_Index, Quantile25 = quantiles.val[,1], Quantile50 = quantiles.val[,2], Quantile75 = quantiles.val[,3], Min = min.val, Max = max.val, Mean = avg.val)
  final.df <- rbind(final.df, to.append)
}


pairwise_comparisons <- combn(labels, 2, simplify = FALSE)
p <- ggplot(final.df, aes(x = Label, y = Max, colour = Label)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(size=2, alpha=0.9) +
  guides(colour=guide_legend(title="Alignment")) +
  theme_light() +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 11), legend.title=element_text(size=14,face="bold"), legend.text=element_text(size=12), legend.position = "none") +
  scale_color_npg() + 
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, vjust=-1, size = 5, aes( label=round(after_stat(y), digits=2))) + 
  stat_summary(fun.y=mean, colour="black", geom="point", shape=18, size=2, show_guide = FALSE) + 
  stat_compare_means(paired = FALSE, size = 5.5, label.x = 1.0, comparisons = pairwise_comparisons, label = "p.signif") +
  labs(y = "log pseudolikelihood", x = "Lineage")
p

ggsave(file="Ecoli_blaCTX-M_log_pseudolikelihood_comparison_all_signif.svg", plot=p, height = 6, width = 10)
