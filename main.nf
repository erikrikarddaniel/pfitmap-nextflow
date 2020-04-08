#!/usr/bin/env nextflow

params.input= "testdata/genomes/*.faa.gz"

genomes = Channel.fromPath(params.input)
      
process singleFaa {

  input: 
  file(fastas) from genomes.collect()

  output:
  file 'all_genomes.faa' into hmmSearch_ch

  shell:
  """
  zcat $fastas > all_genomes.faa
  """
}
