#!/usr/bin/env Rscript

# collect_gffs.R
#
# Author: daniel.lundin@dbb.su.se

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
#suppressPackageStartupMessages(library(dtplyr))
#suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE))
#suppressPackageStartupMessages(library(tidyr))
#suppressPackageStartupMessages(library(stringr))
#suppressPackageStartupMessages(library(magrittr))
#suppressPackageStartupMessages(library(readr))

write("collect_gffs.R starting", stderr())

# opt <- list(options = list(overlap = 1e4), args = c('GCA_001244405.1_Eubacterium_massiliense_genomic.gff.gz', 'unique.accnos'))
option_list = list(
)
opt = parse_args(
  OptionParser(
    usage = "%prog [options]\n\n\tCollects a number of gff tsv files (*.tsv in working directory) output by subset_gff.R into a single one: genomes.tsv.gz.",
    option_list = option_list
  ), 
  positional_arguments = TRUE
)

data <- data.table(
  strand = character(), seqname = character(), source = character(), feature = character(),
  start = integer(), end = integer(), score = character(), frame = character(), 
  attribute = character(), accno = character()
)

n <- 0
for ( f in Sys.glob("*.tsv.gz") ) {
  data <- funion(data, fread(f))
  n <- n + 1
}

write(sprintf("Read %d files, writing to genomes.tsv.gz", n), stderr())

fwrite(data, "genomes.tsv.gz", sep = "\t")

write("Done", stderr())
