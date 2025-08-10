# ==== Библиотеки ====
library(pheatmap)
library(dplyr)
library(RColorBrewer)
library(cluster)
library(stringr)
library(tidyr)
library(tibble)

# ==== Путь к данным/картинке ====
infile  <- "/Users/sofakonopleva/Downloads/concat_table.tsv"

# ==== Загрузка и предобработка ====
data <- read.table(infile, header = TRUE, sep = "\t", stringsAsFactors = FALSE)[-1]

data <- data %>%
  mutate(
    locus = str_extract(locus, "^[^:]+")
  )

# сводим до уникальных комбинаций Species+Pathway и агрегируем по статусу
sp_path <- data %>%
  group_by(Species, Pathway, status) %>%
  summarise(med_log_tpm = median(log_tpm, na.rm = TRUE), .groups = "drop")

# делаем уникальный id строки "Species_Pathway"
data_matrix <- sp_path %>%
  mutate(row_id = paste(Species, Pathway, sep = "_")) %>%
  select(row_id, status, med_log_tpm) %>%
  tidyr::pivot_wider(names_from = status, values_from = med_log_tpm, values_fill = NA) %>%
  tibble::column_to_rownames("row_id") %>%
  as.matrix()

data_matrix <- as.data.frame(data_matrix) %>%
  mutate(across(everything(), ~replace(.x,
                                       !is.finite(.x),
                                       median(.x[is.finite(.x)],
                                              na.rm=TRUE)))) %>%
  as.matrix()

# ==== Центрирование по строкам (даёт синее/красное) ====
data_centered <- sweep(data_matrix, 1, rowMeans(data_matrix, na.rm = TRUE), "-")

row_annotation <- sp_path %>%
  distinct(Species, Pathway) %>%
  mutate(row_id = paste(Species, Pathway, sep = "_")) %>%
  select(row_id, Pathway) %>%
  column_to_rownames("row_id") %>%
  .[rownames(data_matrix), , drop = FALSE]

# ---- Цвета для Pathway ----

pathway_colors <- c(
  "TAL" = "#00BFFF",  # оранжевый
  "TK"  = "#FFD700",  # жёлтый
  "EMP" = "#FF7F50"   # тёмно-синий
)

# формируем список цветов для аннотаций
annotation_colors <- list(Pathway = pathway_colors)

# Пределы шкалы и палитра (сужаем до 2–98 перцентилей, симметрично)
finite_vals <- as.numeric(data_centered[is.finite(data_centered)])
q <- stats::quantile(finite_vals, c(0.02, 0.98), na.rm = TRUE)
lim <- max(abs(q))
breaks <- seq(-lim, lim, length.out = 101)
cols <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)

# ==== Кластеризация по центрированной матрице ====
dist_matrix   <- dist(data_centered)
hclust_result <- hclust(dist_matrix, method = "ward.D2")
clusters      <- cutree(hclust_result, k = 3)
row_annotation <- row_annotation[rownames(data_centered), , drop = FALSE]

# Добавляем кластеры в аннотацию + их цвета
row_annotation$Pathway <- trimws(row_annotation$Pathway)
row_annotation$Cluster <- factor(
  paste("Cluster", clusters[rownames(row_annotation)]),
  levels = paste("Cluster", 1:3)
)
annotation_colors$Cluster <- setNames(brewer.pal(3, "Set1"), paste("Cluster", 1:3))

labels_row <- rownames(data_centered)

# ==== Рисуем ====
png("~/Desktop/heatmap_species_pathway.png", width = 12, height = 10, units = "in", res = 300, bg = "transparent")
pheatmap(data_centered,
         labels_row = rownames(data_matrix),
         annotation_row  = row_annotation,
         annotation_colors = annotation_colors,
         cluster_rows    = hclust_result,
         cluster_cols    = TRUE,
         color           = cols,
         breaks          = breaks,
         na_col          = "grey90",
         fontsize        = 12,
         fontsize_row    = 12,
         fontsize_col    = 15,
         main            = 'Hitmap log_TPM',
         cutree_rows     = 3,
         border_color    = NA)
dev.off()



sp_path$Cluster <- clst


pathway_counts <- sp_path %>%
  group_by(Cluster, Pathway, status) %>%
  summarise(count = n(), .groups = 'drop')


status_cluster_counts <- sp_path %>%
  group_by(Cluster, status, Pathway) %>%
  summarise(count = n(), .groups = 'drop')

# Подсчитываем общее количество в каждом кластере и статусе
status_totals <- status_cluster_counts %>%
  group_by(Cluster, status) %>%
  summarise(total = sum(count), .groups = 'drop')

# Добавляем проценты
status_cluster_percent <- status_cluster_counts %>%
  left_join(status_totals, by = c("Cluster", "status")) %>%
  mutate(percentage = round(count / total * 100, 1))

# Создаем палитру цветов для путей
pathway_colors <- c("EMP" = "#FF6B6B", "TAL" = "#4ECDC4", "TK" = "#45B7D1")

# Стековая гистограмма
p1 <- ggplot(status_cluster_percent, aes(x = status, y = percentage, fill = Pathway)) +
  geom_bar(stat = "identity", position = "stack", alpha = 0.8, color = "white", size = 0.3) +
  geom_text(aes(label = ifelse(percentage >= 5, paste0(percentage, "%"), "")), 
            position = position_stack(vjust = 0.5), 
            color = "white", 
            size =5, 
            fontface = "bold") +
  facet_wrap(~Cluster, scales = "fixed") +
  scale_fill_manual(values = pathway_colors) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 20, face = "bold"), 
    strip.text = element_text(size = 20, face = "bold"),
    legend.position = "top",
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 20)
  ) +
  labs(
    x = "Status (HC - healthy, IBS - sick)",
    y = "Percentage distribution",
    fill = "Pathway"
  )
png("~/Desktop/hist_pathway.png", 
    width = 12, height = 10, units = "in", res = 300,
    bg = "transparent")  # Прозрачный фон
p1
dev.off()


