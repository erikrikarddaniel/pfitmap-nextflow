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
