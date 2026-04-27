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
library(stringr)
library(pdfCluster)
library(aricode)
library(optparse)

entropy <- function(x) {
  p <- table(x) / length(x)
  -sum(p * log(p))
}

parse_filename <- function(filename) {
  list(
    Sampled = if (grepl("sampled", filename)) {
      TRUE
    } else {
      FALSE
    },
    Held_out = if (grepl("held_out_removed", filename)) {
      TRUE
    } else {
      FALSE
    },
    Held_out_only = if (grepl("held_out", filename)) {
      TRUE
    } else {
      FALSE
    },
    Novel = if (grepl("query", filename)) {
      TRUE
    } else {
      FALSE
    },
    Train = if (grepl("train", filename)) {
      TRUE
    } else {
      FALSE
    },
    Test = if (grepl("test", filename)) {
      TRUE
    } else {
      FALSE
    },
    Smitis = if (grepl("mitis", filename)) {
      TRUE
    } else {
      FALSE
    },
    With_unk = if (!grepl("without_unk", filename)) {
      TRUE
    } else {
      FALSE
    },
    Species = if (grepl("E_coli", filename)) {
      "E_coli"
    } else {
      "S_pneumoniae"
    }
  )
}

assign_database_name <- function(df_all){
  
  df_all_sample <- df_all[, c("Species", "Data", "Subset")]
  
  df_all_sample$Species_parsed <- case_when(
    df_all_sample$Species == "E_coli"       ~ "Ecoli",
    df_all_sample$Species == "S_pneumoniae" ~ "Pneumo",
    .default = NA_character_
  )
  df_all_sample$Data_parsed <- gsub("% Data", "", df_all_sample$Data)
  df_all_sample$Subset_parsed <- substr(df_all_sample$Subset, 1, 3)
  
  DatasetName <- paste0(df_all_sample$Species_parsed, df_all_sample$Data_parsed, df_all_sample$Subset_parsed)
  DatasetName <- factor(DatasetName, levels <- c("Pneumo10Sub", "Pneumo100Sub", "Pneumo10All", "Pneumo100All", "Ecoli10Sub", "Ecoli100Sub", "Ecoli10All", "Ecoli100All"))
  DatasetName
}

plot_individual <- function(df_paths, df_paths_embeddings, df_paths_accuracy, outpref)
{
  
  # match up files
  key1 <- basename(df_paths) |>
    str_remove("_leiden_.*\\.tsv$")
  key1_species <- ifelse(grepl("E_coli", df_paths), "E_coli", "S_pneumoniae")
  key1 <- paste(key1, key1_species, sep = "_")
  # key2 <- basename(df_paths_embeddings) |>
  #   str_remove("_leiden_.*\\.csv$")
  
  key2 <- basename(df_paths_embeddings) |>
    str_remove("_UMAP\\.csv$")
  key2_species <- ifelse(grepl("E_coli", df_paths_embeddings), "E_coli", "S_pneumoniae")
  key2 <- paste(key2, key2_species, sep = "_")
  
  key3 <- basename(df_paths_accuracy) |>
    str_remove("_leiden_.*\\.tsv$")
  key3_species <- ifelse(grepl("E_coli", df_paths_accuracy), "E_coli", "S_pneumoniae")
  key3 <- paste(key3, key3_species, sep = "_")
  
  matched <- merge(
    data.frame(key = key1, leiden = df_paths),
    data.frame(key = key2, umap = df_paths_embeddings),
    by = "key"
  )
  
  matched <- merge(
    matched,
    data.frame(key = key3, accuracy = df_paths_accuracy),
    by = "key"
  )
  
  j <- 10
  for (j in 1:nrow(matched))
  {
    row <- matched[j,]
    leiden.path <- row$leiden
    umap.path <- row$umap
    accuracy.path <- row$accuracy
    key <- row$key
    
    df.leiden <- read.table(leiden.path, sep = "\t", comment.char = "", header = TRUE)
    df.leiden$Taxon <- gsub("\\..*", "", df.leiden$Taxon)
    df.umap <- read.table(umap.path, sep = ",", comment.char = "", header = TRUE)
    df.umap$Sample <- gsub("\\..*", "", df.umap$Sample)
    
    df.merged <- merge(df.leiden, df.umap, by.x = "Taxon", by.y = "Sample")
    
    df.accuracy <- read.table(accuracy.path, sep = "\t", comment.char = "", header = TRUE)
    max.row <- df.accuracy[which.max(df.accuracy$AMI),]
    ARI <- max.row$ARI
    AMI <- max.row$AMI
    
    # hold original cluster
    df.merged$predicted_label_ori <- df.merged$predicted_label
    df.merged$parsed_filename <- key
    
    # sort columns, keep top N clusters coloured
    total_colours <- 10
    top_classes <- df.merged %>%
      count(predicted_label, sort = TRUE) %>%
      slice_head(n = total_colours) %>%  # safer than top_n which is deprecated
      pull(predicted_label)
    
    mapping <- setNames(seq_len(total_colours), as.character(top_classes))
    
    df.merged <- df.merged %>%
      mutate(
        predicted_label = ifelse(predicted_label %in% top_classes, mapping[as.character(predicted_label)], "0"),
        predicted_label = factor(predicted_label, levels = c(as.character(seq_len(total_colours)), "0"))
      )
    
    
    # Put "0" last instead of first
    df.merged$predicted_label <- fct_relevel(df.merged$predicted_label, "0", after = Inf)
    
    #Sort rows so "0" comes last (plotted on top)
    df.merged <- df.merged[order(df.merged$predicted_label != "0"), ]
    
    # Build a named colour vector:
    n_top <- nlevels(df.merged$predicted_label) - 1
    # brewer colours
    my_cols <- c(
      setNames(brewer.pal(min(12, n_top), "Set3")[1:n_top], levels(df.merged$predicted_label)[-nlevels(df.merged$predicted_label)]),
      "0" = "gray30"
    )
    
    #ARI <- adj.rand.index(df.merged$predicted_label_ori, df.merged$true_label)
    #AMI <- AMI(factor(as.character(df.merged$predicted_label_ori)), factor(as.character(df.merged$true_label)))
    
    df.merged$ARI <- ARI
    df.merged$AMI <- AMI
    
    #npg colours
    # my_cols <- c(setNames(pal_npg("nrc")(n_top), levels(df$Cluster)[-nlevels(df$Cluster)]),
    #              "0" = "grey")
    
    p <- ggplot(df.merged, aes(x = UMAP1, y = UMAP2, color = predicted_label)) + geom_point(size=0.01, alpha =0.5) +
      theme_light() +
      scale_color_manual(values=my_cols) +
      theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), legend.title=element_text(size=16,face="bold"), legend.text=element_text(size=14), strip.text = element_text(size=12), legend.position = "none") +
      annotate(
        "label",
        x = Inf, y = Inf,   # top right corner
        label = paste0("ARI = ", round(ARI, 3), "\n",
                       "AMI = ", round(AMI, 3)),
        hjust = 1.1, vjust = 1.1,  # nudges inward from the edge
        label.r = unit(0.15, "lines"), # corner rounding
        fill = "white",            # background color
        color = "black",            # text + border color
        size = 3
      ) +
      theme(
        panel.background = element_rect(fill = "black"), # plot background
        #plot.background  = element_rect(fill = "black"), # outer background
        panel.grid.major = element_line(color = "gray40"), # major grid
        panel.grid.minor = element_line(color = "gray20")  # minor grid
      )
    p
    
    # add metadata
    parsed_values <- parse_filename(key)
    for (name in names(parsed_values)) {
      df.merged[[name]] <- parsed_values[[name]]
    }
    
    ggsave(file=paste0(outpref, "_", key, ".png"), plot=p, height = 6, width = 6)
    
    if (j == 1) {
      df_all = df.merged
    } else {
      df_all <- rbind(df_all, df.merged)
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
  
  df_all$DatasetName <- assign_database_name(df_all)
  
  # Get one row per unique 'group'
  test_df <- df_all %>%
    distinct(parsed_filename, .keep_all = TRUE)
  
  df_all
}

plot_all <- function(df, outpref, plot_height, plot_width, svg=FALSE) {
  # top 10 cluster sizes per dataset  
  n_top <- 10
  
  # # sort columns
  # top_classes <- df %>%
  #   group_by(DatasetName) %>%
  #   count(predicted_label, sort = TRUE) %>%
  #   top_n(n_top, n) %>%
  #   pull(predicted_label)
  # 
  # df <- df %>%
  #   group_by(DatasetName) %>%
  #   mutate(predicted_label = ifelse(predicted_label_ori %in% top_classes, predicted_label, "0"))
  
  # sort columns
  # top_classes <- df %>%
  #   count(Cluster, sort = TRUE) %>%
  #   top_n(10, n) %>%
  #   pull(Cluster)
  
  # df <- df %>%
  #   mutate(Cluster = ifelse(Cluster_ori %in% top_classes, Cluster_ori, "0"))
  
  # df <- df %>%
  #   mutate(Cluster = ifelse(Cluster_ori %in% top_classes, Cluster_ori, "0"),
  #          Cluster = factor(Cluster_ori))
  # Put "0" last instead of first
  df$predicted_label <- fct_relevel(df$predicted_label, "0", after = Inf)
  
  #Sort rows so "0" comes last (plotted on top)
  df <- df[order(df$predicted_label != "0"), ]
  
  # brewer colours
  df_ranked <- df %>%
    filter(predicted_label != "0") %>%   # <-- add this
    group_by(Data, Classifier, Subset, Species, predicted_label) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(Data, Classifier, Subset, Species) %>%
    slice_max(count, n = n_top) %>%
    mutate(rank = row_number()) %>%
    ungroup() %>%
    select(Data, Classifier, Subset, Species, predicted_label, rank)
  
  # Join rank back onto main df
  df <- df %>%
    left_join(df_ranked, by = c("Data", "Classifier", "Subset", "Species", "predicted_label")) %>%
    mutate(rank = ifelse(is.na(rank), 0, rank))  # non-top clusters get rank 0
  
  # "0" first = drawn first = appears on bottom
  df <- df %>% arrange(desc(rank == 0), rank)
  
  set3_cols <- brewer.pal(12, "Set3")
  
  rank_cols <- c(
    setNames(
      set3_cols[(seq_along(1:n_top) - 1) %% 12 + 1],  # cycles through 12 colours
      as.character(1:n_top)
    ),
    "0" = "gray30"
  )
  
  # Build a named colour vector:
  # n_top <- nlevels(df$predicted_label) - 1
  # 
  # # brewer colours
  # my_cols <- c(
  #   setNames(brewer.pal(min(12, n_top), "Set3")[1:n_top], levels(df$predicted_label)[-nlevels(df$predicted_label)]),
  #   "0" = "gray30"
  # )
  
  # npg colours
  #my_cols <- c(setNames(pal_npg("nrc")(n_top), levels(df$Cluster)[-nlevels(df$Cluster)]),
  #             "0" = "grey")
  
  annotations <- df %>% 
    group_by(DatasetName, Classifier) %>%
    summarise(AMI = unique(AMI), ARI = unique(ARI))
  
  p <- ggplot(df, aes(x = UMAP1, y = UMAP2, color = as.character(rank))) + geom_point(size=0.001, alpha =0.2) +
    theme_light() +
    facet_grid(DatasetName ~ Classifier) +
    scale_color_manual(values=rank_cols) +
    theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), legend.title=element_text(size=16,face="bold"), legend.text=element_text(size=14), strip.text = element_text(size=12), legend.position = "none") +
    geom_label(
      data = annotations,
      aes(
        x = Inf, y = Inf,
        label = paste0("ARI = ", round(ARI, 3), "\n",
                       "AMI = ", round(AMI, 3)),
      ),
      inherit.aes = FALSE,
      hjust = 1.1, vjust = 1.1,
      fill = "white",
      color = "black",
      size = 3
    ) +
    theme(
      panel.background = element_rect(fill = "black"), # plot background
      #plot.background  = element_rect(fill = "black"), # outer background
      panel.grid.major = element_line(color = "gray40"), # major grid
      panel.grid.minor = element_line(color = "gray20")  # minor grid
    )
  p
  
  ggsave(file=paste0(outpref, ".png"), plot=p, height = plot_height, width = plot_width)
  if (svg)
  {
    ggsave(file=paste0(outpref, ".svg"), plot=p, height = plot_height, width = plot_width)
  }
}

plot_results <- function(df_paths_panbart, df_paths_sketchlib, df_paths_panbart_embeddings, df_paths_sketchlib_embeddings, df_paths_panbart_accuracy, df_paths_sketchlib_accuracy, outpref)
{
  df_all_panbart <- plot_individual(df_paths_panbart, df_paths_panbart_embeddings, df_paths_panbart_accuracy, paste0(outpref, "panbart"))
  df_all_panbart$Classifier <- "PanBART"
  df_all_sketchlib <- plot_individual(df_paths_sketchlib, df_paths_sketchlib_embeddings, df_paths_sketchlib_accuracy, paste0(outpref, "sketchlib"))
  df_all_sketchlib$Classifier <- "Sketchlib"
  
  merged_df <- rbind(df_all_sketchlib, df_all_panbart)
  
  # training only, subset - fine
  subset_df <- subset(merged_df, Query == "Training" & Unknown_token == "With unknown" & Subset == "Subset lineages")
  outpref_plot <- paste0(outpref, "leiden_training_subset_lineages")
  plot_all(subset_df, outpref_plot, 14, 10)
  
  # # training only, all, 100% data
  # subset_df <- subset(merged_df, Query == "Training" & Unknown_token == "With unknown" & Subset != "Subset lineages" & Data == "100% Data")
  # outpref_plot <- paste0(outpref, "leiden_training_all_lineages_1.0_data")
  # plot_all(subset_df, outpref_plot)
  # 
  # training only, all - fine
  subset_df <- subset(merged_df, Query == "Training" & Unknown_token == "With unknown" & Subset != "Subset lineages")
  outpref_plot <- paste0(outpref, "leiden_training_all_lineages")
  plot_all(subset_df, outpref_plot, 14, 10)
  
  # testing only, subset - fine
  subset_df <- subset(merged_df, Query == "Testing" & Unknown_token == "With unknown" & Subset == "Subset lineages")
  outpref_plot <- paste0(outpref, "leiden_testing_subset_lineages")
  plot_all(subset_df, outpref_plot, 14, 10)
  
  # testing only, all - fine
  subset_df <- subset(merged_df, Query == "Testing" & Unknown_token == "With unknown" & Subset != "Subset lineages")
  outpref_plot <- paste0(outpref, "leiden_testing_all_lineages")
  plot_all(subset_df, outpref_plot, 14, 10)
  
  # held-out only, with unknown
  subset_df <- subset(merged_df, Query == "Novel" & Unknown_token == "With unknown" & Subset == "Subset lineages")
  outpref_plot <- paste0(outpref, "leiden_novel_subset_lineages_with_unknown")
  plot_all(subset_df, outpref_plot, 14, 10)
  
  # held only, without unknown
  # subset_df <- subset(merged_df, (Query == "Novel" & Unknown_token != "With unknown" & Subset == "Subset lineages") | (Query == "Novel" & Classifier == "Sketchlib" & Subset == "Subset lineages"))
  # outpref_plot <- paste0(outpref, "leiden_novel_subset_lineages_without_unknown")
  # plot_all(subset_df, outpref_plot)
}

option_list <- list(
  make_option("--indir", type="character", help="Input directory containing panbart/ and sketchlib/ subdirectories with species subdirectories."),
  make_option("--outdir", type="character", help="Output directory for saved plots."),
  make_option("--species", type="character", default="S_pneumoniae,E_coli", help="Comma-separated species subdirectory names. Default = 'S_pneumoniae,E_coli'")
)
opt <- parse_args(OptionParser(option_list=option_list))

indir <- opt$indir
outdir <- opt$outdir
species_list <- strsplit(opt$species, ",")[[1]]


outpref <- outdir

# # get highest scoring ARI with K, leiden resolution
# j <- 1
# for (j in 1:length(df_paths_panbart))
# {
#   filename <- df_paths_panbart[j]
#   base <- gsub(indir, "", filename)
#   df <- read.csv(filename, sep="\t")
#   max.row <- df[which.max(df$AMI), ]
#   
#   max.row$file <- base
#   if (j == 1)
#   {
#     total.df <- max.row
#   } else {
#     total.df <- rbind(total.df, max.row)
#   }
# }
# 
# for (j in 1:length(df_paths_sketchlib))
# {
#   filename <- df_paths_sketchlib[j]
#   base <- gsub(indir, "", filename)
#   df <- read.csv(filename, sep="\t")
#   max.row <- df[which.max(df$AMI), ]
#   
#   max.row$file <- base
#   
#   total.df <- rbind(total.df, max.row)
# }

# df_paths_panbart <- list.files(path = paste0(indir, "panbart/", species), pattern = "*_true.tsv", full.names = TRUE)
# df_paths_sketchlib <- list.files(path = paste0(indir, "sketchlib/", species), pattern = "*_true.tsv", full.names = TRUE)
# df_paths_panbart_embeddings <- list.files(path = paste0(indir, "panbart/", species, "/embeddings"), pattern = "*_UMAP.csv", full.names = TRUE)
# df_paths_sketchlib_embeddings <- list.files(path = paste0(indir, "sketchlib/", species, "/embeddings"), pattern = "*_UMAP.csv", full.names = TRUE)
#df_paths_panbart_accuracy <- list.files(path = paste0(indir, "panbart/", species), pattern = "*_per_iter_accuracy.tsv", full.names = TRUE)
#df_paths_sketchlib_accuracy <- list.files(path = paste0(indir, "sketchlib/", species), pattern = "*_per_iter_accuracy.tsv", full.names = TRUE)

df_paths_panbart <- unlist(lapply(species_list, function(sp) {
  list.files(
    path = paste0(indir, "panbart/", sp),
    pattern = "*_true.tsv",
    full.names = TRUE
  )
}))
df_paths_sketchlib <- unlist(lapply(species_list, function(sp) {
  list.files(
    path = paste0(indir, "sketchlib/", sp),
    pattern = "*_true.tsv",
    full.names = TRUE
  )
}))
df_paths_panbart_embeddings <- unlist(lapply(species_list, function(sp) {
  list.files(
    path = paste0(indir, "panbart/", sp, "/embeddings"),
    pattern = "*_UMAP.csv",
    full.names = TRUE
  )
}))
df_paths_sketchlib_embeddings <- unlist(lapply(species_list, function(sp) {
  list.files(
    path = paste0(indir, "sketchlib/", sp, "/embeddings"),
    pattern = "*_UMAP.csv",
    full.names = TRUE
  )
}))
df_paths_panbart_accuracy <- unlist(lapply(species_list, function(sp) {
  list.files(
    path = paste0(indir, "panbart/", sp),
    pattern = "*_per_iter_accuracy.tsv",
    full.names = TRUE
  )
}))
df_paths_sketchlib_accuracy <- unlist(lapply(species_list, function(sp) {
  list.files(
    path = paste0(indir, "sketchlib/", sp),
    pattern = "*_per_iter_accuracy.tsv",
    full.names = TRUE
  )
}))

#TESTING
df_paths <- df_paths_panbart
df_paths_embeddings <- df_paths_panbart_embeddings
df_paths_accuracy <- df_paths_panbart_accuracy

plot_results(df_paths_panbart, df_paths_sketchlib, df_paths_panbart_embeddings, df_paths_sketchlib_embeddings, df_paths_panbart_accuracy, df_paths_sketchlib_accuracy, outpref)
