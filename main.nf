#!/usr/bin/env nextflow

// This should be moved to a config file for tests
params.genomes = "testdata/genomes/*.faa.gz"
params.hmms    = "testdata/profiles/*.hmm"

genomes   = Channel.fromPath(params.genomes)
hmm_files = Channel.fromPath(params.hmms)
      
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

process hmmsearch {
  publishDir "/home/ghada/dev/pfitmap_nextflow/testdata", mode: 'copy'

  input:
  file genomes from all_genomes_ch
  file hmm     from hmm_files
  
  output:
  file "${hmm_aln.baseName}.hmm" into hmm_profiles_ch
  file "*.tblout" into tblout_ch
  file "*.domtblout" into domtblout_ch

  shell:
  """
  hmmsearch --tblout "${hmm.basename}.tblout" --domtblout "${hmm.basename}.domtblout" $hmm $genomes
  """
}
