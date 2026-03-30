## load packages
req_packages <- c("janitor", "optparse", "tidyverse")
invisible(suppressWarnings(suppressMessages(
    lapply(req_packages, require, character.only = TRUE)
)))

## store commandline arguments with optparse package
# option_list = list(
#     make_option(c("-v", "--vcf"), type="character", default=NULL, 
#               help="Variant Call File (VCF) path", metavar="character"),
#     make_option(c("-x", "--skip"), type="numeric", default=NULL, 
#               help="Header length of input VCF", metavar="numeric"),
#     make_option(c("-g", "--gtf"), type="character", default=NULL, 
#               help="General Feature Format (GFF) file path", metavar="character"),
#     make_option(c("-l", "--snplocation"), type="character", default="snpsloc.txt", 
#               help="Output SNP location file path [default= %default]", metavar="character"),
#     make_option(c("-s", "--snpgenotype"), type="character", default="SNP.txt", 
#               help="Output SNP genotype file path [default= %default]", metavar="character"),
#     make_option(c("-h", "--genelocation"), type="character", default="geneloc.txt", 
#               help="Output gene location file path [default= %default]", metavar="character"),
# ); 

# opt_parser = OptionParser(option_list=option_list);
# opt = parse_args(opt_parser);

## store commandline arguments with base R
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 8) {

    opt <- list()
    opt["vcf"] <- args[1]
    opt["skip"] <- as.numeric(args[2])
    opt["expression"] <- args[3]
    opt["gtf"] <- args[4]
    opt["snplocation"] <- args[5]
    opt["snpgenotype"] <- args[6]
    opt["genelocation"] <- args[7]
    opt["reads"] <- args[8]

} else {

    error_message <- "R_matrixeqtl_input_zwa.R: Create inputs for MatrixEQTL from VCF and GTF \nUsage: R_matrixeqtl_input_zwa.R $path/to/vcf $number.of.vcf.lines.to.skip $path/to/gene/expression/input $path/to/gtf $path/to/snploc/output $path/to/snp/output $path/to/geneloc/output $path/to/read/output"
    stop(error_message)

}


## load in VCF
vcf <- read_tsv(opt$vcf, skip = opt$skip) %>%
    janitor::clean_names()

## get location of SNPs
snpsloc <- vcf %>%
    rownames_to_column("snp") %>%
    mutate(snp = paste("Snp", snp, sep = "_")) %>%
    select(snp, number_chrom, pos) %>%
    rename(chr = number_chrom) %>%
    unique()

## load in gene expression dataframe
gene <- read_tsv(opt$expression)

## get genotype of SNPs in samples
snp <- vcf[,11:ncol(vcf)] %>%
    mutate(across(everything(), ~ str_remove(., ":.*"))) %>%
    mutate(across(everything(), ~ case_when(. == "0/0" ~ 0,
                                            . == "0/1" | . == "1/0" ~ 1,
                                            . ==  "1/1" ~ 2,
                                            TRUE ~ NA))) %>%
    rownames_to_column("id") %>%
    mutate(id = paste("Snp", id, sep = "_")) %>%
    filter(id %in% snpsloc$snp)

## reformat names to match gene expression data
colnames(snp) <- str_to_upper(colnames(snp))
colnames(snp)[1] <- "id"

## add replicates to dataframe and only keep columns with gene expression data
snp <- snp %>%
  pivot_longer(!id, names_to = "strain", values_to = "snp") %>%
  crossing(rep = paste0("rep", 1:3)) %>%
  mutate(strain_rep = paste0(strain, "_", rep)) %>%
  select(-c(strain, rep)) %>%
  pivot_wider(names_from = strain_rep, values_from = snp) %>%
  select(id, any_of(colnames(gene)))

reads <- gene %>%
    select(gene, any_of(colnames(snp)))

## load in GFF
gtf <- read_tsv(opt$gtf, col_names = FALSE) %>%
    janitor::clean_names()

## pull gene name and location
geneloc <- gtf %>%
    rename(type = 3, format = 9, chr = 1, s1 = 4, s2 = 5) %>%
    filter(type == "exon") %>%
    select(format, chr, s1, s2) %>%
    separate_wider_delim(format, delim = ";", names = c("transcript_id", "geneid", "gene_name", "drop")) %>%
    mutate(gene_name = str_remove_all(gene_name, " gene_name "),
           gene_name = str_remove_all(gene_name, "\"")) %>%
    select(gene_name, chr, s1, s2)

## save files
write_tsv(snpsloc, opt$snplocation)
write_tsv(snp, opt$snpgenotype)
write_tsv(geneloc, opt$genelocation)
write_tsv(reads, opt$reads)