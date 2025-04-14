library(ggplot2)
library(ggsci)

df<-read.csv("AtB_All_S_pneumoniae_held_in_model_mask0.15_lr0.00001_drop0.4_parsed_training_log.txt", sep = "\t", header = 1)

data_subset = subset(df, Type != "Test")
p <- ggplot(data_subset, aes(x = Epoch, y=Loss, colour = Type)) +
  geom_line(linewidth=1.5) +
  scale_colour_npg() +
  theme_light() +
  theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), legend.title=element_text(size=16,face="bold"), legend.text=element_text(size=14)) +
  guides(colour=guide_legend(title="Data type")) +
  scale_y_log10() +
  labs(x = "Epoch", y = "Cross Entropy Loss")
p
ggsave(file="All_Pneumo_PanBART_loss.svg", plot=p, height = 6, width = 8)
ggsave(file="All_Pneumo_PanBART_loss.png", plot=p, height = 6, width = 8)

p <- ggplot(data_subset, aes(x = Epoch, y=Precision, colour = Type)) +
  geom_line(linewidth=1.5) +
  scale_colour_npg() +
  theme_light() +
  theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), legend.title=element_text(size=16,face="bold"), legend.text=element_text(size=14), legend.position = "none") +
  guides(colour=guide_legend(title="Data type")) +
  labs(x = "Epoch", y = "Precision")
p
ggsave(file="All_Pneumo_PanBART_precision.svg", plot=p, height = 6, width = 8)
ggsave(file="All_Pneumo_PanBART_precision.png", plot=p, height = 6, width = 8)

p <- ggplot(data_subset, aes(x = Epoch, y=Recall, colour = Type)) +
  geom_line(linewidth=1.5) +
  scale_colour_npg() +
  theme_light() +
  theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), legend.title=element_text(size=16,face="bold"), legend.text=element_text(size=14), legend.position = "none") +
  guides(colour=guide_legend(title="Data type")) +
  labs(x = "Epoch", y = "Recall")
p
ggsave(file="All_Pneumo_PanBART_recall.svg", plot=p, height = 6, width = 8)
ggsave(file="All_Pneumo_PanBART_recall.png", plot=p, height = 6, width = 8)

p <- ggplot(data_subset, aes(x = Epoch, y=Kappa, colour = Type)) +
  geom_line(linewidth=1.5) +
  scale_colour_npg() +
  theme_light() +
  theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), legend.title=element_text(size=16,face="bold"), legend.text=element_text(size=14), legend.position = "none") +
  guides(colour=guide_legend(title="Data type")) +
  labs(x = "Epoch", y = "Cohen's Kappa")
p
ggsave(file="All_Pneumo_PanBART_kappa.svg", plot=p, height = 6, width = 8)
ggsave(file="All_Pneumo_PanBART_kappa.png", plot=p, height = 6, width = 8)
