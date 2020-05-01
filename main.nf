#!/usr/bin/env nextflow

/**
 * main.nf: Nextflow workflow for GTDB
 *
 * The GTDB genomes are downloaded and annotated
 * The workflow starts from a set of annotated genomes in the format of faa.gz files
 * Requirements: directory with all hmm profiles to be run 
 *
 *   Concatenate all faa.gz files into a single one
 *   Performs an hmm_search of all hmm profiles on all the genomes
 *   Downloads the metadata files for Archaea and Bacterial genomes from gtdb latest version repository and concatenates them into a single metadata file
 *   Classify the found 
 *
 * ghada.nouraia@dbb.su.se daniel.lundin@dbb.su.se
 */

// Parameters
params.help                     = false
params.inputgenomes             = ''
params.profiles_hierarchy       = ''
params.dbsource                 = ''
params.hmm_mincov               = 0.9
params.gtdb_arc_metadata        = ''
params.gtdb_bac_metadata        = ''
     
def helpMessage() {
  log.info """

  Usage:

  The typical command for running the pipeline is as follows:

  nextflow run main.nf --inputgenomes path/to/genomes --outputdir path/to/results --hmm_mincov value --dbsource NCBI:NR:*date*

  Mandatory arguments:
  --inputgenomes path/to/genomes		Path to annotated genomes in the format faa.gz 
  --outputdir path/to/results 			Path to the results directory
  --hmm_mincov value 				Set a value for the threshold of coverage hmm_profile/query (default = 0.9)
  --dbsource GTDB:release			Set the database source in the format GTDB:release, where 'release' mention which GTDB release was used
  """.stripIndent()
}

// Show help message
if (params.help) {
  helpMessage()
  exit 0
}

// Create channels to start processing
System.err.println(sprintf("genomes: %s", params.inputgenomes))
genomes                 = Channel.fromPath(params.inputgenomes)
System.err.println("hmms")
hmm_files               = Channel.fromPath(params.hmms)
System.err.println("profiles")
profiles_hierarchy      = Channel.fromPath(params.profiles_hierarchy)
dbsource                = Channel.value(params.dbsource)
hmm_mincov              = Channel.value(params.hmm_mincov)
System.err.println("arc")
gtdb_arc_metadata       = Channel.fromPath(params.gtdb_arc_metadata)
System.err.println("bac")
gtdb_bac_metadata       = Channel.fromPath(params.gtdb_bac_metadata)
results                 = params.outputdir

process singleFaa {
  input: 
  file("*.faa.gz") from genomes.collect()

  output:
  file 'all_genomes.faa' into all_genomes_hmmsearch_ch
  file 'all_genomes.faa' into all_genomes_classify_ch

  shell:
  """
  for f in \$(find . -name '*.faa.gz'|xargs readlink); do
    a=\$(basename \$f | sed 's/\\..*//');
    gunzip -c \$f | sed "/^>/s/\$/ [\$a]/";
  done > all_genomes.faa
  """
}

process hmmSearch {
  publishDir results, mode: 'copy'

  input:
  file genome from all_genomes_hmmsearch_ch
  file hmm     from hmm_files
  
  output:
  file ("*.tblout") into tblout_ch
  file ("*.domtblout") into domtblout_ch

  shell:
  """
  hmmsearch --tblout "${hmm.baseName}.tblout" --domtblout "${hmm.baseName}.domtblout" $hmm $genome
  """
}

process getMetadata {
  publishDir results, mode: 'copy'

  input:
  file arc_metadata from gtdb_arc_metadata
  file bac_metadata from gtdb_bac_metadata

  output:
  file 'gtdb_metadata.tsv' into gtdbmetadata_ch

  shell:
  """
  cat $arc_metadata > gtdb_metadata.tsv
  cat $bac_metadata |grep -v ^accession >> gtdb_metadata.tsv
  """
}

process pfClassify {
  publishDir results, mode: 'copy'

  input:
  val hmm_mincov from hmm_mincov
  val dbsource from dbsource
  file "hmm_profile_hierarchy.tsv" from profiles_hierarchy
  file "gtdb_metadata.tsv" from gtdbmetadata_ch
  file tblout from tblout_ch
  file domtblout from domtblout_ch
  file genomes from all_genomes_classify_ch

  output:
  file "gtdb.pf.db" into gtdb_pf_db_ch
  file "gtdb.tsv.gz" into gtdb_tsv_ch
  file "gtdb.pf-classify.warnings.txt" into gtdb_classify_warnings_ch

  shell:
  """
  pf-classify.r --hmm_mincov=${hmm_mincov} --dbsource=${dbsource} --gtdbmetadata=gtdb_metadata.tsv --profilehierarchies=hmm_profile_hierarchy.tsv --singletable=gtdb.tsv.gz --seqfaa=${genomes} --sqlitedb=gtdb.pf.db  $tblout $domtblout > gtdb.pf-classify.warnings.txt 2>&1
  """
}

process formatSequences {
  publishDir results, mode: 'copy'
  input:
  file "all_genomes.faa" from all_genomes_ch

  output:
  file "all_genomes.c.faa" into formatted_genomes_ch
  
  shell:
  """
  cat "all_genomes.faa" | awk '/^>/ {printf("\\n%s\\n",\$\$0); next;} {printf("%s",\$\$0);} END {printf("\\n");}' | sed '/^\$\$/d' | sed '/^>/s/[ \\t]\\+/_/g' | sed '/^>/s/[-(),:;+|]//g' | sed '/^>/!s/ //g' > all_genomes.c.faa
  """
}

