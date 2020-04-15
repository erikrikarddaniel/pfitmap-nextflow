/**
 * main.nf: Nextflow workflow for GTDB
 *
 * The GTDB genomes are downloaded and annotated
 * The workflow starts from a set of annotated genomes in the format of faa.gz files
 * Requirements: directory with all hmm profiles to be run 
 *
 *   Concatenate all faa.gz files into a single onee
 *   Performs an hmm_search of all hmm profiles on all the genomes
 *   Downloads the metadata files for Archaea and Bacterial genomes from gtdb latest version repository and concatenates them into a single metadata file
 *   Classify the found 
 *
 * daniel.lundin@lnu.se
 */

#!/usr/bin/env nextflow
params.help = false
genomes   = Channel.fromPath(params.inputgenomes)
hmm_files = Channel.fromPath(params.hmms)
profiles_hierarchy = Channel.fromPath(params.profiles_hierarchy)
DBSOURCE = Channel.value(params.DBSOURCE)
hmm_mincov = Channel.value(params.hmm_mincov)
results = params.outputdir
     
def helpMessage() {
  log.info """

  Usage:

  The typical command for running the pipeline is as follows:

  nextflow run main.nf --inputgenomes path/to/genomes --outputdir path/to/results --hmm_mincov value --DBSOURCE NCBI:NR:*date*

  Mandatory arguments:
  --inputgenomes path/to/genomes		Path to annotated genomes in the format faa.gz 
  --outputdir path/to/results 			Path to the results directory
  --hmm_mincov value 				Set a value for the threshold of coverage hmm_profile/querry (default = 0.9)
  --DBSOURCE GTDB:release			Set the database source in the format GTDB:release, where 'release' mention which GTDB release was used
  """.stripIndent()
}

// Show help message
if (params.help) {
  helpMessage()
  exit 0
}

process singleFaa {
  input: 
  file("*.faa.gz") from genomes.collect()

  output:
  file 'all_genomes.faa' into all_genomes_ch

  shell:
  """
  find . -name “*.faa.gz” | xargs gunzip -c | gzip -c > all_genomes.faa.gz
  
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
