library(readr)
library(ggplot2)
library(randomForest)

df<-read.csv("synteny_grid_summary.txt", sep = "\t", header = 0)
colnames(df) <- c("param_set_ID", "batch_size", "block_size", "n_layer", "n_head", "n_embd", "dropout", "learning_rate", "max_epochs", 
                  "weight_decay", "beta1", "beta2", "grad_clip", "decay_lr", "warmup_epochs", "lr_decay_epochs", 
                  "min_lr", "best_train_loss", "best_train_perp", "best_val_loss", "best_val_perp")

p <- ggplot(df, aes(x = param_set_ID, y = best_val_loss)) + geom_bar(stat = "identity", fill = "skyblue") + ylab("Validation loss") + xlab("Parameter set ID")
p

# Or using random forests (requires the randomForest package)
df_training <- df[, -which(names(df) %in% c("param_set_ID", "best_train_loss", "best_train_perp", "best_val_perp"))]
rf_model <- randomForest(best_val_loss ~ ., data = df_training)
importance(rf_model)

best_model <- df[which.min(df$best_val_loss),]
