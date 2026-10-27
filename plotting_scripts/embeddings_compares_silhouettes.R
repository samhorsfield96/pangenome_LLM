library(ggplot2)
library(ggsci)
library(RColorBrewer)
library(ggpubr)
library(reporter)
library(tools)
library(dplyr)
library(forcats)
library(reshape2)
library(tidyr)
library(dplyr)
library(optparse)

parse_filename <- function(filename) {
  list(
    Sampled = if (grepl("_sampled_", filename)) {
      TRUE
    } else {
      FALSE
    },
    Held_out = if (grepl("_held_out_removed", filename)) {
      TRUE
    } else {
      FALSE
    },
    Held_out_only = if (grepl("_held_out_", filename)) {
      TRUE
    } else {
      FALSE
    },
    Novel = if (grepl("_query_", filename)) {
      TRUE
    } else {
      FALSE
    },
    Train = if (grepl("_train_", filename)) {
      TRUE
    } else {
      FALSE
    },
    Test = if (grepl("_test_", filename)) {
      TRUE
    } else {
      FALSE
    },
    Smitis = if (grepl("_mitis_", filename)) {
      TRUE
    } else {
      FALSE
    },
    With_unk = if (!grepl("_without_unk_", filename)) {
      TRUE
    } else {
      FALSE
    }
  )
}

plot_individual <- function(df_paths, min_count, num_bins)
{
  j <- 4
  for (j in 1:length(df_paths))
  {
    # get silhouette values
    filename_silhouette <- df_paths[j]
    df <- read.table(filename_silhouette, sep = "\t", comment.char = "", header = TRUE)
    colnames(df) <- c("Cluster", "Silhouette_width", "Count")
    
    # remove rows with low counts
    df <- df[!df$Count <= min_count, ]
    df <- df[!df$Cluster == "overall", ]
    
    max_count <- max(df$Count)
    df$Proportion <- df$Count / max_count
    
    parsed_filename <- file_path_sans_ext(basename(filename_silhouette))
    df$parsed_filename <- parsed_filename
    
    # add metadata
    parsed_values <- parse_filename(filename_silhouette)
    for (name in names(parsed_values)) {
      df[[name]] <- parsed_values[[name]]
    }
    
    # Use hist() to define breaks and counts
    breaks = seq(0, 1.0, by = 1.0 / num_bins)
    h <- hist(df$Proportion, breaks = breaks, plot = FALSE)
    
    # Use cut() to assign each Count to a bin
    df$bin <- cut(df$Proportion, breaks = h$breaks, include.lowest = TRUE, right = TRUE)
    
    # Extract the bin start and end for each row
    # First, make a data frame mapping factor levels to start/end
    bin_map <- data.frame(
      bin = levels(df$bin),
      bin_start = head(h$breaks, -1),
      bin_end = tail(h$breaks, -1)
    )
    
    # Merge bin information back to df
    df <- merge(df, bin_map, by = "bin", all.x = TRUE)
    
    # Optional: reorder to original order
    df <- df[order(as.numeric(row.names(df))), ]
    
    if (j == 1) {
      df_all <- df
    } else {
      df_all <- rbind(df_all, df)
    }
  }
  
  df_all['Unknown_token'] <- NA
  df_all['Data'] <- NA
  df_all['Subset'] <- NA
  df_all['Query'] <- NA
  
  df_all['Subset'][df_all['Held_out'] != TRUE ] <- "All lineages"
  df_all['Subset'][df_all['Held_out'] == TRUE | df_all['Held_out_only'] == TRUE] <- "Subset lineages"
  
  df_all['Unknown_token'][df_all['With_unk'] == TRUE] <- "With unknown"
  df_all['Unknown_token'][df_all['With_unk'] == FALSE] <- "Without unknown"
  
  df_all['Data'][df_all['Sampled'] == TRUE] <- "10% Data"
  df_all['Data'][df_all['Sampled'] == FALSE] <- "100% Data"
  
  df_all['Query'][df_all['Smitis'] == TRUE] <- "S. mitis"
  df_all['Query'][df_all['Train'] == TRUE] <- "Training"
  df_all['Query'][df_all['Test'] == TRUE] <- "Testing"
  df_all['Query'][df_all['Held_out'] != TRUE & df_all['Held_out_only'] == TRUE] <- "Novel"
  df_all['Query'][df_all['Novel'] == TRUE & df_all['Train'] == TRUE] <- "Novel vs. Training"
  df_all['Query'][df_all['Novel'] == TRUE & df_all['Train'] != TRUE] <- "Novel"
  df_all['Query'][df_all['Test'] == TRUE & df_all['Train'] == TRUE] <- "Testing vs. Training"
  df_all['Query'][df_all['Smitis'] == TRUE & df_all['Train'] == TRUE ] <- "S. mitis vs. Training"
  df_all['Query'][df_all['Smitis'] == TRUE & df_all['Train'] != TRUE] <- "S. mitis"
  
  df_all
}

plot_all <- function(df, outpref, top_K) {
    
    df$Cluster_ori <- df$Cluster
    
    if (top_K > 0) {
      # sort columns
      top_classes <- df %>%
        count(Cluster, sort = TRUE) %>%
        top_n(10, n) %>%
        pull(Cluster)
      
      df <- df %>%
        mutate(Cluster = ifelse(Cluster_ori %in% top_classes, Cluster_ori, "0"))
      
      df <- subset(df, Cluster != "0")
    }
    
    # should be 2 for each
    check <- table(df$Cluster, df$Data)
    unique(check[1])
    unique(check[2])
    unique(check[3])
    
    subset_df_PanBART <- subset(df, Classifier == "PanBART" )
    subset_df_Sketchlib <- subset(df, Classifier == "Sketchlib" )
    
    merged_df <- merge(subset_df_PanBART, subset_df_Sketchlib, by = c("Cluster", "Data", "bin"))
    
    test.stats <- merged_df %>%
      group_by(Data, bin) %>%
      summarise(p_value = wilcox.test(Silhouette_width.x, Silhouette_width.y, paired = TRUE)$p.value,
                counts = n())
    
    p <- ggplot(df, aes(x = Classifier, y = Silhouette_width, fill = Classifier)) + geom_violin() +
      theme_light() +
      facet_grid(bin ~ Data) +
      scale_fill_npg() +
      theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), legend.title=element_text(size=16,face="bold"), legend.text=element_text(size=14), strip.text = element_text(size=12), legend.position = "none") +
      xlab("Tool") + ylab("Silhouette score") +
      stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, vjust=-1, size = 5, aes( label=round(after_stat(y), digits=2))) + 
      stat_summary(fun.y=mean, colour="black", geom="point", shape=18, size=2, show_guide = FALSE) + 
      #stat_compare_means(aes(group = Cluster), paired = TRUE, size = 5.5, label.x = 1.25, label.y = 0.75, method = "wilcox.test")
      geom_text(
        data = test.stats,
        aes(
          x = 1.75, y = Inf,
          label = paste0("p=", formatC(p_value, format = "e", digits = 2), "\n",
                         "N pairs=", counts),
        ),
        inherit.aes = FALSE,
        hjust = 1.1, vjust = 1.1,
        size = 3.5
      )
    p
    
    ggsave(file=paste0(outpref, ".png"), plot=p, height = 6, width = 10)
    
    # plot silhouette vs. frequency
    p <- ggplot(df, aes(x = Count, y = Silhouette_width, colour = Classifier, group= Classifier)) + geom_point(alpha=0.3) +
      geom_smooth(method = "gam") +
      theme_light() +
      facet_grid(Data ~ .) +
      scale_fill_npg() +
      scale_x_log10() +
      theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), legend.title=element_text(size=16,face="bold"), legend.text=element_text(size=14), strip.text = element_text(size=12),) +
      xlab("Cluster size") + ylab("Silhouette score") 
    p
    ggsave(file=paste0(outpref, "_frequency.png"), plot=p, height = 6, width = 10)
    
}

plot_results <- function(df_paths_panbart_silhouette, df_paths_sketchlib_silhouette, outpref, min_count, num_bins)
{
  df_all_panbart <- plot_individual(df_paths_panbart_silhouette, min_count, num_bins)
  df_all_panbart$Classifier <- "PanBART"
  df_all_sketchlib <- plot_individual(df_paths_sketchlib_silhouette, min_count, num_bins)
  df_all_sketchlib$Classifier <- "Sketchlib"
  
  merged_df <- rbind(df_all_sketchlib, df_all_panbart)
  
  # training only, subset
  subset_df <- subset(merged_df, Query == "Training" & Unknown_token == "With unknown" & Subset == "Subset lineages")
  outpref_plot <- paste0(outpref, "_training_subset_lineages")
  plot_all(subset_df, outpref_plot, 0)
  
  # training only, all
  subset_df <- subset(merged_df, Query == "Training" & Unknown_token == "With unknown" & Subset != "Subset lineages")
  outpref_plot <- paste0(outpref, "_training_all_lineages")
  plot_all(subset_df, outpref_plot, 0)
  
  # testing only, subset
  subset_df <- subset(merged_df, Query == "Testing" & Unknown_token == "With unknown" & Subset == "Subset lineages")
  outpref_plot <- paste0(outpref, "_testing_subset_lineages")
  plot_all(subset_df, outpref_plot, 0)
  
  # testing only, all
  subset_df <- subset(merged_df, Query == "Testing" & Unknown_token == "With unknown" & Subset != "Subset lineages")
  outpref_plot <- paste0(outpref, "_testing_all_lineages")
  plot_all(subset_df, outpref_plot, 0)
  
  # held only, subset TODO, looks wrong
  subset_df <- subset(merged_df, Query == "Novel" & Unknown_token == "With unknown" & Subset == "Subset lineages")
  outpref_plot <- paste0(outpref, "_novel_subset_lineages_with_unknown")
  plot_all(subset_df, outpref_plot, 0)
  
  # held only, subset TODO, looks wrong
  subset_df <- subset(merged_df, Query == "Novel" & Unknown_token == "With unknown" & Subset != "Subset lineages")
  outpref_plot <- paste0(outpref, "_novel_all_lineages_with_unknown")
  plot_all(subset_df, outpref_plot, 0)
}

option_list <- list(
  make_option("--indir", type="character", help="Input directory containing panbart/ and sketchlib/ subdirectories with *_silhouette.tsv files."),
  make_option("--outpref", type="character", default="output", help="Output prefix for saved plots. Default = 'output'"),
  make_option("--min-count", type="integer", default=1, help="Minimum cluster count to include. Default = 1"),
  make_option("--num-bins", type="integer", default=3, help="Number of bins for proportion histogram. Default = 3")
)
opt <- parse_args(OptionParser(option_list=option_list))

indir <- opt$indir
min_count <- opt$`min-count`
num_bins <- opt$`num-bins`

df_paths_panbart_silhouette <- list.files(path = paste0(indir, "/panbart/"), pattern = "*_silhouette.tsv", full.names = TRUE)
df_paths_sketchlib_silhouette <- list.files(path = paste0(indir, "/sketchlib/"), pattern = "*_silhouette.tsv", full.names = TRUE)
outpref <- opt$outpref

plot_results(df_paths_panbart_silhouette, df_paths_sketchlib_silhouette, outpref, min_count, num_bins)
