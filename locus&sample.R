library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(tibble)

# Читаем данные
data <- read.table("concat_table.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

data <- data %>%
  mutate(
    taxon = str_extract(Lineage, "p__[^;_]+"),
    locus = str_extract(locus, "^[^:]+")
  )

# Создаем матрицу pivot с log_tpm значениями
heatmap_matrix <- data %>%
  select(locus, Run, log_tpm) %>%
  pivot_wider(names_from = Run, values_from = log_tpm, values_fill = 0) %>%
  column_to_rownames("locus") %>%
  as.matrix()

# Подготавливаем аннотации для строк (locus - Pathway)
row_annotation <- data %>%
  select(locus, Pathway) %>%
  distinct() %>%
  column_to_rownames("locus")

# Подготавливаем аннотации для колонок (Run - status)
col_annotation <- data %>%
  select(Run, status) %>%
  distinct() %>%
  column_to_rownames("Run")

# Создаем цветовые палитры для аннотаций
# Для Pathway (строки)
pathway_colors <- rainbow(length(unique(data$Pathway)))
names(pathway_colors) <- unique(data$Pathway)

# Для status (колонки)
status_colors <- brewer.pal(min(8, length(unique(data$status))), "Set3")
names(status_colors) <- unique(data$status)

# Объединяем цвета для аннотаций
ann_colors <- list(
  Pathway = pathway_colors,
  status = status_colors
)

# Создаем heatmap с log_tpm значениями
pheatmap(
  heatmap_matrix,
  annotation_row = row_annotation,      # Аннотация строк (locus -> Pathway)
  annotation_col = col_annotation,      # Аннотация колонок (Run -> status)
  annotation_colors = ann_colors,       # Цвета для аннотаций
  cluster_rows = TRUE,                  # Кластеризация строк
  cluster_cols = TRUE,                  # Кластеризация колонок
  show_rownames = TRUE,                 # Показать названия строк
  show_colnames = TRUE,                 # Показать названия колонок
  scale = "row",                        # Масштабирование по строкам (Z-score)
  color = colorRampPalette(c("blue", "white", "red"))(100),  # Цветовая палитра
  fontsize = 10,                        # Размер шрифта
  fontsize_row = 8,                     # Размер шрифта для названий строк
  fontsize_col = 8,                     # Размер шрифта для названий колонок
  main = "Экспрессия генов (log_tpm): локусы vs сэмплы"  # Заголовок
)

