#!/usr/bin/env nextflow

genomes   = Channel.fromPath(params.genomes)
hmm_files = Channel.fromPath(params.hmms)
results = params.results
     
process singleFaa {
  input: 
  file("*.faa.gz") from genomes.collect()

  output:
  file 'all_genomes.faa' into all_genomes_ch

  shell:
  """
  zcat *.faa.gz > all_genomes.faa
  
  """
}

process hmmSearch {
  publishDir results, mode: 'copy'

  input:
  file genome from all_genomes_ch
  file hmm     from hmm_files
  
  output:
  file "*.tblout" into tblout_ch
  file "*.domtblout" into domtblout_ch

  shell:
  """
  hmmsearch --tblout "${hmm.baseName}.tblout" --domtblout "${hmm.baseName}.domtblout" $hmm $genome
  """
}

process getMetadata {
  
  output:
  file 'gtdb_metadata.tsv' into gtdbmetadata_ch

  shell:
  """
  wget https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/ar122_metadata.tsv
  wget https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/bac120_metadata.tsv
  cat ar122_metadata.tsv > gtdb_metadata.tsv
  cat bac120_metadata.tsv |grep -v ^accession >> gtdb_metadata.tsv
  """
}

