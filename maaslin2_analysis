library(readr)
library(dplyr)
library(tidyr)
library(Maaslin2)
library(tibble)

data1 <- data %>%
    select(assemble, locus, tpm, Host_age, host_sex, gastrointest_disord)
table <- data1 %>%
    select(assemble, locus, tpm) %>%
    pivot_wider(names_from = locus, values_from = tpm)
metadata <- data1 %>%
    select(assemble, Host_age, host_sex, gastrointest_disord) %>%
    distinct()  
    
table_no_na <- table
table_no_na[is.na(table_no_na)] <- 0    
rownames(table_no_na) <- table_no_na$assemble
rownames(table_no_na_df) <- table_no_na_df$assemble
metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata$assemble

features_only <- table_no_na %>%
    select(-assemble)
features_matrix <- t(as.matrix(features_only))
colnames(features_matrix) <- table_no_na$assemble
input <- as.data.frame(features_matrix)
input_log <- log(input + 1)


fit <- Maaslin2(
    input_data = input,
    input_metadata = metadata,
    output = "maaslin2_output",
    analysis_method = 'LM',
    fixed_effects = c("gastrointest_disord"),
    reference = c("gastrointest_disord:Control"),
    normalization = "NONE",
    transform = "NONE",
    standardize = TRUE,
    cores = 1,
    min_abundance = 1,    
    min_prevalence = 0.1  
    )
