library(stringr)
library(dplyr)
library(tidyr)
library(dplyr)
library(optparse)

parse_filename <- function(filename) {
  list(
    Sampled = if (grepl("sampled", filename)) {
      TRUE
    } else {
      FALSE
    },
    Held_out = if (grepl("held_out", filename)) {
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
    })
}

generate_table <- function(paths) {
  j <- 2
  for (j in 1:length(paths)) {
    filename <- paths[j]
    df <- read.csv(filename, sep = "\t")
    
    if (nrow(df) == 0)
    {
      next
    }
    
    # remove duplicate rows
    df <- df[!duplicated(df), ]
    
    base <- basename(filename)
    target_gene <- str_match(base, "_([^_]+)_N")[,2]
    
    df$target_gene <- target_gene
    df$filename <- base
    
    parsed_values <- parse_filename(base)
    for (name in names(parsed_values)) {
      df[[name]] <- parsed_values[[name]]
    }
    
    if (j == 1) {
      df_all <- df
    } else {
      df_all <- rbind(df_all, df)
    }
  }
  
  df_all['Data'] <- NA
  df_all['Subset'] <- NA
  df_all['Query'] <- NA
  
  df_all['Subset'][df_all['Held_out'] != TRUE ] <- "All lineages"
  df_all['Subset'][df_all['Held_out'] == TRUE] <- "Subset lineages"
  
  df_all['Data'][df_all['Sampled'] == TRUE] <- "10% Data"
  df_all['Data'][df_all['Sampled'] == FALSE] <- "100% Data"
  
  df_all['Query'][df_all['Train'] == TRUE] <- "Training"
  df_all['Query'][df_all['Test'] == TRUE] <- "Testing"
  
  df_all <- df_all %>%
    separate_rows(UniRef, sep = ",") %>%
    mutate(
      level = str_extract(UniRef, "^UniRef[0-9]+")
    ) %>%
    filter(is.na(level) | level != "") %>%
    pivot_wider(
      names_from  = level,
      values_from = UniRef,
      values_fill = NA
    )
  
  df_all
}

option_list <- list(
  make_option("--species", type="character", help="Species name (subdirectory within --indir-pos and --indir-neg)."),
  make_option("--indir-pos", type="character", help="Directory containing species subdirectory with positive SHAP annotation .tsv files."),
  make_option("--indir-neg", type="character", help="Directory containing species subdirectory with negative SHAP annotation .tsv files."),
  make_option("--outdir", type="character", help="Output directory for saved tables.")
)
opt <- parse_args(OptionParser(option_list=option_list))

species <- opt$species
indir_pos <- opt$`indir-pos`
paths_pos <- list.files(path = paste0(indir_pos, species), pattern = paste0("*hits_annotated_pos.tsv"), full.names = TRUE)

indir_neg <- opt$`indir-neg`
paths_neg <- list.files(path = paste0(indir_neg, species), pattern = paste0("*hits_annotated_neg.tsv"), full.names = TRUE)

outdir <- opt$outdir

pos_df <- generate_table(paths_pos)
neg_df <- generate_table(paths_neg)

write.csv(pos_df, paste0(outdir, species, "_pos_df_all.csv"))
write.csv(neg_df, paste0(outdir, species, "_neg_df_all.csv"))

# remove matching tokens
pos_df_filtered <- pos_df %>%
  anti_join(neg_df, by = "Token")

neg_df_filtered <- neg_df %>%
  anti_join(pos_df, by = "Token")

write.csv(pos_df_filtered, paste0(outdir, species, "_pos_df_filtered_all.csv"))
write.csv(neg_df_filtered, paste0(outdir, species, "_neg_df_filtered_all.csv"))

# write all individual tables
{
  #df_all <- pos_df
  #type <- "pos"
  
  #df_all <- neg_df
  #type <- "neg"
  
  #df_all <- pos_df_filtered
  #type <- "pos_filtered"
  
  df_all <- neg_df_filtered
  type <- "neg_filtered"
  
  # tokens and genes
  count_df_token_genes <- df_all %>%
    group_by(target_gene, Token, Gene) %>%
    summarise(count = n()) %>%   # count rows per group/item combo
    arrange(target_gene, desc(count))
  
  write.csv(count_df_token_genes, paste0(outdir, species, "_", type, "_count_df_token_genes.csv"))
  
  # genes only
  count_df_genes <- df_all %>%
    group_by(target_gene, Gene) %>%
    summarise(count = n()) %>%   # count rows per group/item combo
    arrange(target_gene, desc(count))
  
  write.csv(count_df_genes, paste0(outdir, species, "_", type, "_count_df_genes.csv"))
  
  count_df_product <- df_all %>%
    group_by(target_gene, Product) %>%
    summarise(count = n()) %>%   # count rows per group/item combo
    arrange(target_gene, desc(count))
  
  write.csv(count_df_product, paste0(outdir, species, "_", type, "_count_df_product.csv"))
  
  count_df_token_product <- df_all %>%
    group_by(target_gene, Token, Product) %>%
    summarise(count = n()) %>%   # count rows per group/item combo
    arrange(target_gene, desc(count))
  
  write.csv(count_df_token_product, paste0(outdir, species, "_", type, "_count_df_token_product.csv"))
  
  count_df_uniref50 <- df_all %>%
    group_by(target_gene, UniRef50) %>%
    summarise(count = n()) %>%   # count rows per group/item combo
    arrange(target_gene, desc(count))
  
  write.csv(count_df_uniref50, paste0(outdir, species, "_", type, "_count_df_uniref50.csv"))
  
  count_df_token_uniref50 <- df_all %>%
    group_by(target_gene, Token, UniRef50) %>%
    summarise(count = n()) %>%   # count rows per group/item combo
    arrange(target_gene, desc(count))
  
  write.csv(count_df_token_uniref50, paste0(outdir, species, "_", type, "_count_df_token_uniref50.csv"))
  
  count_df_product_uniref50 <- df_all %>%
    group_by(target_gene, Product, UniRef50) %>%
    summarise(count = n()) %>%   # count rows per group/item combo
    arrange(target_gene, desc(count))
  
  write.csv(count_df_product_uniref50, paste0(outdir, species, "_", type, "_count_df_product_uniref50.csv"))
  
  count_df_uniref90 <- df_all %>%
    group_by(target_gene, UniRef90) %>%
    summarise(count = n()) %>%   # count rows per group/item combo
    arrange(target_gene, desc(count))
  
  write.csv(count_df_uniref90, paste0(outdir, species, "_", type, "_count_df_uniref90.csv"))
  
  count_df_token_uniref90 <- df_all %>%
    group_by(target_gene, Token, UniRef90) %>%
    summarise(count = n()) %>%   # count rows per group/item combo
    arrange(target_gene, desc(count))
  
  write.csv(count_df_token_uniref90, paste0(outdir, species, "_", type, "_count_df_token_uniref90.csv"))
  
  count_df_product_uniref90 <- df_all %>%
    group_by(target_gene, Product, UniRef90) %>%
    summarise(count = n()) %>%   # count rows per group/item combo
    arrange(target_gene, desc(count))
  
  write.csv(count_df_product_uniref90, paste0(outdir, species, "_", type, "_count_df_product_uniref90.csv"))
  
  count_df_uniref100 <- df_all %>%
    group_by(target_gene, UniRef100) %>%
    summarise(count = n()) %>%   # count rows per group/item combo
    arrange(target_gene, desc(count))
  
  write.csv(count_df_uniref100, paste0(outdir, species, "_", type, "_count_df_uniref100.csv"))
  
  count_df_token_uniref100 <- df_all %>%
    group_by(target_gene, Token, UniRef90) %>%
    summarise(count = n()) %>%   # count rows per group/item combo
    arrange(target_gene, desc(count))
  
  write.csv(count_df_token_uniref100, paste0(outdir, species, "_", type, "_count_df_token_uniref100.csv"))
  
  count_df_product_uniref100 <- df_all %>%
    group_by(target_gene, Product, UniRef100) %>%
    summarise(count = n()) %>%   # count rows per group/item combo
    arrange(target_gene, desc(count))
  
  write.csv(count_df_product_uniref100, paste0(outdir, species, "_", type, "_count_df_product_uniref100.csv"))
}
