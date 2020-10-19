#!/usr/bin/env Rscript

# subset_gff.R
#
# Author: daniel.lundin@dbb.su.se

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dtplyr))
suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(readr))

# Range to include. Entries starting in the interval from each target gene minus OVERLAP and
# plus OVERLAP, will be included in the output.
OVERLAP <- 1e5

# Get arguments
# opt <- list(options = list(overlap = 1e4), args = c('d', 'd/unique.accnos'))
option_list = list(
  make_option(
    "--overlap", type = 'integer', default = OVERLAP, action = 'store', help = 'Size of overlap, default %default.'
  )
)
opt = parse_args(
  OptionParser(
    usage = "%prog [options] genome.gff.gz acc_list\n\n\tSubsets genome.gff.gz to only contain genes nearby (see --overlap) accessions in acc_list.",
    option_list = option_list
  ), 
  positional_arguments = TRUE
)

write(sprintf("Subsetting gffs in %s with %s", opt$args[1], opt$args[2]), stderr())

accnos <- fread(opt$args[2], col.names = c('accno'))
matches <- data.table(
  strand = character(), seqname = character(), source = character(), feature = character(),
  start = integer(), end = integer(), score = character(), frame = character(),
  attribute = character(), accno = character()
)

for ( f in Sys.glob(sprintf("%s/*.gff.gz", opt$args[1])) ) {
  print(sprintf("--> %s <--", f))
  gff <- fread(
    cmd = sprintf("gunzip -c %s | grep -v '#' | grep '\t'", f),
    col.names = c('seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute')
  ) %>%
    lazy_dt() %>%
    mutate(accno = str_extract(attribute, "ID=[^;]+") %>% str_remove("ID=")) %>%
    as.data.table()
  setkey(gff, strand, start, end)

  # Find matches
  m <- lazy_dt(accnos) %>% inner_join(lazy_dt(gff), by = 'accno') %>% as.data.table()
  m <- m[, .(strand, intvstart = start - opt$options$overlap, intvend = end + opt$options$overlap)]

  setkey(m, strand, intvstart, intvend)

  m <- foverlaps(gff, m, type = 'within', by.x = c('strand', 'start', 'end'), by.y = c('strand', 'intvstart', 'intvend'))[!is.na(intvstart)] %>%
    lazy_dt() %>% select(-intvstart, -intvend) %>%
    as.data.table()

  matches <- funion(matches, m)
}

fwrite(matches, "genomes.tsv.gz")

write("Done", stderr())
