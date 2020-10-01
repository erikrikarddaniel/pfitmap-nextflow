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

# Range to include. Entries starting in the interval from each target gene minus OVERLAP and
# plus OVERLAP, will be included in the output.
OVERLAP <- 1e5

# Get arguments
# opt <- list(options = list(overlap = 1e4), args = c('GCA_001244405.1_Eubacterium_massiliense_genomic.gff.gz', 'unique.accnos'))
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

write(sprintf("Subsetting %s with %s", opt$args[1], opt$args[2]), stderr())

gff <- fread(
  cmd = sprintf("gunzip -c %s | grep -v '#' | grep '\t'", opt$args[1]),
  col.names = c('seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute')
) %>%
  lazy_dt() %>%
  mutate(accno = str_extract(attribute, "ID=[^;]+") %>% str_remove("ID=")) %>%
  as.data.table()
setkey(gff, strand, start, end)
accnos <- fread(opt$args[2], col.names = c('accno'))

# Find matches
matches <- lazy_dt(accnos) %>% inner_join(lazy_dt(gff), by = 'accno') %>% 
  as.data.table()
matches$start <- matches$start - opt$options$overlap
matches$end   <- matches$end   + opt$options$overlap
matches <- matches[, .(strand, start, end)]
setkey(matches, strand, start, end)

foverlaps(gff, matches, type = 'within', by.x = c('strand', 'start', 'end'), by.y = c('strand', 'start', 'end'))[!is.na(start)] -> t

write("Done", stderr())
