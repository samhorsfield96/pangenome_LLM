library(ggplot2)
library(ggsci)
library(ggpubr)
library(rlang)
library(tools)
library(optparse)

search_string <- function(filename, string){
  search_rg <- paste0(".*", string, "_(-?[0-9.]+(?:[eE][-+]?[0-9]+)?)(?:_|\\.|$).*")
  label = if (grepl(string, filename)) {
    as.numeric(sub(search_rg, "\\1", filename))
  } else {
    NA
  }
  return(label)
}

parse_filename <- function(filename) {
  list(
    Sampled = if (grepl("_sampled_", filename)) {
      TRUE
    } else {
      FALSE
    },
    Held_out = if (grepl("_held_out_removed_", filename)) {
      TRUE
    } else {
      FALSE
    },
    Held_out_only = if (grepl("_held_out_only_", filename)) {
      TRUE
    } else {
      FALSE
    },
    With_unk = if (grepl("with_unk", filename)) {
      TRUE
    } else {
      FALSE
    },
    Sketchlib = if (grepl("sketchlib", filename)) {
      TRUE
    } else {
      FALSE
    }
  )
}

parse_df <- function(df_paths) {
  
  j <- 27
  for (j in 1:length(df_paths))
  {
    filename <- df_paths[j]
    df <- read.table(filename, sep = "\t", comment.char = "", header = TRUE)
    
    colnames(df)[which(names(df) == "Label")] <- "Cluster"
    df <- df[!df$Cluster == "overall", ]
    df <- df[!df$Cluster == "overall_balanced", ]
    
    df$filename <- parsed_filename <- file_path_sans_ext(basename(filename))
    
    # add metadata
    parsed_values <- parse_filename(filename)
    for (name in names(parsed_values)) {
      df[[name]] <- parsed_values[[name]]
    }
    
    if (j == 1)
    {
      df_all <- df
    } else {
      df_all <- rbind(df_all, df)
    }
  }
  
  df_all['Data'] <- NA
  df_all['Lineages'] <- NA
  df_all['Unknown_token'] <- NA
  df_all['Classifier'] <- NA
  
  df_all['Data'][df_all['Sampled'] == TRUE] <- "10%"
  df_all['Data'][df_all['Sampled'] == FALSE] <- "100%"
  
  df_all['Lineages'][df_all['Held_out'] == TRUE] <- "Subset lineages"
  df_all['Lineages'][df_all['Held_out'] == FALSE] <- "All lineages"
  df_all['Lineages'][df_all['Held_out_only'] == TRUE] <- "Held-out lineages"
  
  df_all['Unknown_token'][df_all['With_unk'] == TRUE] <- "With unknown"
  df_all['Unknown_token'][df_all['With_unk'] == FALSE] <- "Without unknown"
  
  df_all['Classifier'][df_all['Sketchlib'] == TRUE] <- "Sketchlib"
  df_all['Classifier'][df_all['Sketchlib'] == FALSE] <- "PanBART"
  
  
  df_all$Data <- factor(df_all$Data, levels = c("10%", "100%"))
  df_all$Lineages <- factor(df_all$Lineages, levels = c("All lineages", "Subset lineages", "Held-out lineages"))
  df_all$Classifier <- factor(df_all$Classifier, levels = c("PanBART", "Sketchlib"))
  
  df_all
}

generate_plots <- function(data_subset, outpref)
{
  p <- ggplot(data_subset, aes(x = Test_count, y=Accuracy, colour = Classifier)) +
    #geom_line(linewidth=1.5) +
    geom_point(alpha=0.1) +
    geom_smooth(method = "gam")+
    scale_colour_npg() +
    theme_light() +
    facet_grid(K ~ Data) +
    scale_x_log10() +
    theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), legend.title=element_text(size=16,face="bold"), legend.text=element_text(size=14), strip.text = element_text(size=12)) +
    #guides(colour=guide_legend(title=legend_title)) +
    labs(x = "Cluster size", y = "Accuracy") +
    geom_hline(yintercept = 1, linetype="dashed")
  p
  ggsave(file=paste(outpref, "_frequency.png", sep = ""), plot=p, height = 10, width = 8)
}

plot_results <- function(df_paths_panbart, df_paths_sketchlib, outpref) {
  df_panbart <- parse_df(df_paths_panbart)
  df_sketchlib <- parse_df(df_paths_sketchlib)
  
  df_merged <- rbind(df_panbart, df_sketchlib)
  
  df_test <- subset(df_merged, Classifier == "Sketchlib")
  
  # all data, with unknown
  data_subset = subset(df_merged, Lineages == "All lineages" & ((Classifier == "PanBART" & Unknown_token != "Without unknown") | (Classifier == "Sketchlib" & Unknown_token == "Without unknown")))
  outpref_plot <- paste0(outpref, "_all_lineages_w_unknown")
  generate_plots(data_subset, outpref_plot)
  
  # Subset lineages with unknown
  data_subset = subset(df_merged, Lineages == "Subset lineages" & ((Classifier == "PanBART" & Unknown_token != "Without unknown") | (Classifier == "Sketchlib" & Unknown_token == "Without unknown")))
  outpref_plot <- paste0(outpref, "_held_out_removed_w_unknown")
  generate_plots(data_subset, outpref_plot)
}

option_list <- list(
  make_option("--indir", type="character", help="Input directory containing panbart/ and sketchlib/ subdirectories with *per_label_accuracy.tsv files."),
  make_option("--outpref", type="character", default="output", help="Output prefix for saved plots. Default = 'output'")
)
opt <- parse_args(OptionParser(option_list=option_list))

indir <- opt$indir

{
  df_paths_panbart <- list.files(path = paste0(indir, "/panbart/"), pattern = "*per_label_accuracy.tsv", full.names = TRUE)
  df_paths_sketchlib <- list.files(path = paste0(indir, "/sketchlib/"), pattern = "*per_label_accuracy.tsv", full.names = TRUE)
  outpref <- opt$outpref
  plot_results(df_paths_panbart, df_paths_sketchlib, outpref)
}
