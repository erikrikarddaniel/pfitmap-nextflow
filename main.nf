#!/usr/bin/env nextflow

// This should be moved to a config file for tests
params.genomes = "testdata/genomes/*.faa.gz"

genomes = Channel.fromPath(params.genomes)
      
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
