suppressMessages(library(tidyverse))

## load in arguments
args <- commandArgs(trailingOnly = TRUE)

## throw error if not enough arguments
if (length(args) != 4) {

    stop("Usage: add_end.R {path/original/gff} {#skip} {path/lengths} {path/new/gff}")
}

## store arguments as variables
original_gff_path <- args[1]
skip <- args[2]
chromosome_length_path <- args[3]
new_gff_path <- args[4]

## load in files
original_gff <- read_tsv(original_gff_path, skip = as.numeric(skip), col_names = FALSE) %>%
    rename(seqid = 1, source = 2, type = 3, start = 4, end = 5, score = 6, strand = 7, phase = 8, attributes = 9) %>%
    suppressMessages()
chromosome_length <- read_tsv(chromosome_length_path, col_names = FALSE) %>%
    rename(seqid = 1, length = 2) %>%
    suppressMessages()

## combine data frames and add chromosome length
new_gff <- original_gff %>%
    inner_join(chromosome_length, by = "seqid") %>%
    mutate(end = case_when(source == "RepeatMasker" ~ length,
                           TRUE ~ end))

## select only the gff columns
new_gff <- new_gff[ , 1:9]

## save new gff
write_tsv(new_gff, new_gff_path, escape = "none")