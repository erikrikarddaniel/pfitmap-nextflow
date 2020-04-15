#!/usr/bin/env nextflow

genomes   = Channel.fromPath(params.genomes)
hmm_files = Channel.fromPath(params.hmms)
profiles_hierarchy = Channel.fromPath(params.profiles_hierarchy)
DBSOURCE = Channel.value(params.DBSOURCE)
hmm_mincov = Channel.value(params.hmm_mincov)
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

process pfClassify {
 
  input:
  val hmm_mincov from hmm_mincov
  val DBSOURCE form DBSOURCE
  file "hmm_profile_hierarchy.tsv" from profiles_hierarchy
  file "gtdb_metadata.tsv" from gtdbmetadata_ch
  file ("*.tblout") from tblout_ch
  file ("*.domtblout") from domtblout_ch
  
  output:
  file 'gtdb.tsv.gz' into gtdb.tsv_ch
  file 'gtdb.pf.db' into gtdb.pf.db_ch

  shell:
  """
  pf-classify.r --verbose --hmm_mincov=${hmm_mincov} --dbsource=${DBSOURCE} --gtdbmetadata=gtdb_metadata.tsv --profilehierarchies=hmm_profile_hierarchy.tsv --singletable=gtdb.tsv.gz --sqlitedb=gtdb.pf.db  *.tblout *.domtblout
  """
}
