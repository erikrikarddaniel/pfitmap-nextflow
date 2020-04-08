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

}

process hmmsearch {
  publishDir "/home/ghada/dev/pfitmap_nextflow/testdata", mode: 'copy'

  input:
  file "HMMER_DB_FAA" from hmmSearch_ch
  file "hmm_aln" from Channel.fromPath("/home/ghada/dev/pfitmap-nextflow/testdata/profiles/*alnfaa")
  
  output:
  file "${hmm_aln.baseName}.hmm" into hmm_profiles_ch
  file "*.tblout" into tblout_ch
  file "*.domtblout" into domtblout_ch

  shell:
  """
  hmmbuild $hmm_aln $hmm_profiles_ch
  hmmsearch --tblout "${hmm_profiles.basename}.tblout" --domtblout "${hmm_profiles.basename}.domtblout" $hmm_profiles_ch $HMMER_DB_FAA
  """
}
