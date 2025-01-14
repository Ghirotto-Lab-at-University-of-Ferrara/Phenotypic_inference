
library(ggplot2)
library(cowplot)
library(ggrepel)

args<-commandArgs(trailingOnly=TRUE)

sample_id<-as.character(args[1])
sample_csv <- read.csv(paste(sample_id,"_phenotypicPrediction.csv",sep=""),header=T, sep=";") 

# Conteggio delle previsioni dei colori degli occhi
eye_colour_counts <- table(sample_csv$predicted_eye_colour)

# Creazione del dataframe per il grafico a torta degli occhi
eye_colour_df <- data.frame(
  Eye_color = names(eye_colour_counts),
  Count = as.numeric(eye_colour_counts)
)

eye_colour_level <- c("Brown", "Intermediate", "Blue", "Not predicted")

# Rinomina la colonna nel dataframe eye_colour_df
names(eye_colour_df)[1] <- "Eye_colour"

# Esegui la merge dei dataframe
eye_colour_df <- merge(eye_colour_df, data.frame(Eye_colour = eye_colour_level), 
                       by = "Eye_colour", all = TRUE)

eye_colour_df$Count[is.na(eye_colour_df$Count)] <- 0  # Imposta il conteggio a zero per le nuove categorie

eye_colour_df$Percent <- eye_colour_df$Count / sum(eye_colour_df$Count) * 100  # Calcola la percentuale

eye_colour_df$Eye_colour <- factor(eye_colour_df$Eye_colour, levels = eye_colour_level)


# Plot a torta per gli occhi
eye_pie <- ggplot(eye_colour_df, aes(x = "", y = Count, fill = Eye_colour)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0) +
  theme_void() + 
  labs(title = "Eye colour phenotype") +
  scale_fill_manual(values = c("Blue" = "#3D85C7", 
                               "Intermediate" = "#597133", 
                               "Brown" = "#713012", 
                               "Not predicted" = "#BFBFBF"),
                    name = "Eye colour") +  
  geom_label_repel(data = subset(eye_colour_df, Percent > 0),
                   aes(label = paste0(round(Percent, 1), "%")),
                   position = position_stack(vjust = 0.5),
                   colour = "white", box.padding = 0.2, force = 0.01,
                   segment.color = "white", segment.size = 1,
                   point.padding = 0.5, show.legend = FALSE) +
  theme(legend.key.width = unit(0.3, "cm"), legend.key.height = unit(0.3, "cm"),
        plot.title = element_text(hjust = 0.5, color = "black", face = "bold", size = 12),
        legend.title = element_text(color = "black", face = "bold", size = 10),
        legend.text = element_text(size = 8))



# Conteggio delle previsioni dei colori dei capelli
hair_colour_counts <- table(sample_csv$predicted_hair_colour)

# Creazione del dataframe per il grafico a torta degi capelli
hair_colour_df <- data.frame(
  Hair_color = names(hair_colour_counts),
  Count = as.numeric(hair_colour_counts)
)

hair_colour_level <- c("Black", "Dark brown/black", "Brown/dark brown", "Brown", "Dark blond", "Blond", "Red", "Not predicted")

# Rinomina la colonna nel dataframe hair_colour_df
names(hair_colour_df)[1] <- "Hair_colour"

# Esegui la merge dei dataframe
hair_colour_df <- merge(hair_colour_df, data.frame(Hair_colour = hair_colour_level), 
                        by = "Hair_colour", all = TRUE)

hair_colour_df$Count[is.na(hair_colour_df$Count)] <- 0  # Imposta il conteggio a zero per le nuove categorie

hair_colour_df$Percent <- hair_colour_df$Count / sum(hair_colour_df$Count) * 100  # Calcola la percentuale

hair_colour_df$Hair_colour <- factor(hair_colour_df$Hair_colour, levels = hair_colour_level)


# Plot a torta per gli occhi
hair_pie <- ggplot(hair_colour_df, aes(x = "", y = Count, fill = Hair_colour)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0) +
  theme_void() + 
  labs(title = "Hair colour phenotype") +
  scale_fill_manual(values = c("Black" = "#000000",
                               "Dark brown/black" = "#1c0e04",
                               "Brown/dark brown" = "#381c08",
                               "Brown" = "#6f370f",
                               "Dark blond" = "#7e5414",
                               "Blond" = "#dcba7f",
                               "Red" = "#cf5f17",
                               "Not predicted" = "#BFBFBF"),
                    name = "Hair colour") +  
  geom_label_repel(data = subset(hair_colour_df, Percent > 0),
                   aes(label = paste0(round(Percent, 1), "%")),
                   position = position_stack(vjust = 0.5),
                   colour = "white", box.padding = 0.2, force = 0.01,
                   segment.color = "white", segment.size = 1,
                   point.padding = 0.5, show.legend = FALSE) +
  theme(legend.key.width = unit(0.3, "cm"), legend.key.height = unit(0.3, "cm"),
        plot.title = element_text(hjust = 0.5, color = "black", face = "bold", size = 12),
        legend.title = element_text(color = "black", face = "bold", size = 10),
        legend.text = element_text(size = 8))



# Conteggio delle previsioni dei colori degli occhi
skin_colour_counts <- table(sample_csv$predicted_skin_colour)

# Creazione del dataframe per il grafico a torta degli occhi
skin_colour_df <- data.frame(
  Skin_color = names(skin_colour_counts),
  Count = as.numeric(skin_colour_counts)
)

skin_colour_level <- c("Dark to black", "Dark", "Intermediate", "Pale", "Very Pale", "Not predicted")

# Rinomina la colonna nel dataframe skin_colour_df
names(skin_colour_df)[1] <- "Skin_colour"

# Esegui la merge dei dataframe
skin_colour_df <- merge(skin_colour_df, data.frame(Skin_colour = skin_colour_level), 
                        by = "Skin_colour", all = TRUE)

skin_colour_df$Count[is.na(skin_colour_df$Count)] <- 0  # Imposta il conteggio a zero per le nuove categorie

skin_colour_df$Percent <- skin_colour_df$Count / sum(skin_colour_df$Count) * 100  # Calcola la percentuale

skin_colour_df$Skin_colour <- factor(skin_colour_df$Skin_colour, levels = skin_colour_level)


# Plot a torta per gli occhi
skin_pie <- ggplot(skin_colour_df, aes(x = "", y = Count, fill = Skin_colour)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0) +
  theme_void() + 
  labs(title = "Skin colour phenotype") +
  scale_fill_manual(values = c("Dark to black" = "#452F17",
                               "Dark" = "#81582B",
                               "Intermediate" = "#C18747",
                               "Pale" = "#D5AD81",
                               "Very Pale" = "#EAD6C0",
                               "Not predicted" = "#BFBFBF"),
                    name = "Skin colour") +  
  geom_label_repel(data = subset(skin_colour_df, Percent > 0),
                   aes(label = paste0(round(Percent, 1), "%")),
                   position = position_stack(vjust = 0.5),
                   colour = "white", box.padding = 0.2, force = 0.01,
                   segment.color = "white", segment.size = 1,
                   point.padding = 0.5, show.legend = FALSE) +
  theme(legend.key.width = unit(0.3, "cm"), legend.key.height = unit(0.3, "cm"),
        plot.title = element_text(hjust = 0.5, color = "black", face = "bold", size = 12),
        legend.title = element_text(color = "black", face = "bold", size = 10),
        legend.text = element_text(size = 8))




# Creazione dei singoli insiemi di grafici (Plot 1)

plot_title <- ggdraw() + draw_label(paste(sample_id, "Phenotypes",sep=" "),
                                    size = 18, fontface = "bold") + 
  theme(plot.margin = margin(0, 0, 0, 0))

plot_combined <- plot_grid(eye_pie, hair_pie, skin_pie,
                           ncol = 3, align = "hv")

# Costruzione della griglia combinata
final_plot <- plot_grid(plot_title, plot_combined, ncol = 1, rel_heights = c(0.1, 0.5))



ggsave(filename=paste(sample_id,"_phenotypicPredictionGL.pdf",sep=""), 
       plot = final_plot,
       width = 297, 
       height = 210, 
       units = "mm")

ggsave(paste(sample_id,"_phenotypicPredictionGL_eye.jpg",sep=""), plot = eye_pie, width = 3.5, height = 3, dpi = 300)
ggsave(paste(sample_id,"_phenotypicPredictionGL_hair.jpg",sep=""), plot = hair_pie, width = 3.5, height = 3, dpi = 300)
ggsave(paste(sample_id,"_phenotypicPredictionGL_skin.jpg",sep=""), plot = skin_pie, width = 3.5, height = 3, dpi = 300)

