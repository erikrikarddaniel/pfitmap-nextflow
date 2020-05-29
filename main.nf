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
params.inputgenomes             = null
params.profiles_hierarchy       = null
params.dbsource                 = 'GTDB:GTDB:latest'
params.hmm_mincov               = 0.9
params.gtdb_arc_metadata        = null
params.gtdb_bac_metadata        = null

params.max_cpus = 2
params.max_time = "240.h"

def helpMessage() {
  log.info """

  Usage:

  The typical command for running the pipeline is as follows:

  nextflow run main.nf --inputgenomes path/to/genomes --outputdir path/to/results --hmm_mincov value --dbsource GTDB:GTDB:release

  Mandatory arguments:
  --inputgenomes path/to/genomes_directory		Path and name of the directory containing annotated genomes in the format faa.gz 
  --gtdb_bac_metadata path/to/file			Path and name of tsv file including the metadata for bacterial genomes
  --gtdb_arc_metadata path/to/file 			Path and name of tsv file including the metadata for archaeal genomes
  --hmm							Path to the HMM profile files 
  --profiles_hierarchy	path/to/file			Path and name of tsv file including hmm profile names and information (See README file for more details)		
  --hmm_mincov value					Set a value for the threshold of coverage hmm_profile/querry (default = 0.9)
  --dbsource db:db:release				Set the database source in the format db:db:release, where [db] is the name of the database and [release] mentions 
									the release number/name (default = GTDB:GTDB:latest)
  --outputdir path/to/results				Path to the results directory

  Non Mandatory parameters:
  --max_cpus						Maximum number of CPU cores to be used (default = 2)
  --max_time						Maximum time per process (default = 10 days)
  
  """.stripIndent()
}

// Show help message
if (params.help) {
  helpMessage()
  exit 0
}

if (! ( params.dbsource =~ /.+:.+:.+/ ) ) { error "Error in dbsource format\ndbsource should be in the format db:db:release\nSee more using --help"}

if( !params.inputgenomes ) {
  error "Missing inputgenomes parameter\n[Parameter error] Please specify the parameter --inputgenomes\nSee more using --help" 
}

if( !params.profiles_hierarchy ) {
  error "Missing profiles_hierarchy parameter\n[Parameter error] Please specify the parameter --profiles_hierarchy\nSee more using --help"
}

if( !params.gtdb_arc_metadata ) {
  error "Missing gtdb_arc_metadata parameter\n[Parameter error] Please specify the parameter --gtdb_arc_metadata\nSee more using --help"
}

if( !params.gtdb_bac_metadata ) {
  error "Missing gtdb_bac_metadata parameter\n[Parameter error] Please specify the parameter --gtdb_bac_metadata\nSee more using --help"
}

// Create channels to start processing

genomes   = Channel.fromPath(params.inputgenomes, checkIfExists : true)
hmm_files = Channel.fromPath("$params.hmms/*.hmm") //, checkIfExists : true)
profiles_hierarchy = Channel.fromPath(params.profiles_hierarchy, checkIfExists : true)
dbsource = Channel.value(params.dbsource)
hmm_mincov = Channel.value(params.hmm_mincov)
gtdb_arc_metadata = Channel.fromPath(params.gtdb_arc_metadata, checkIfExists : true)
gtdb_bac_metadata = Channel.fromPath(params.gtdb_bac_metadata, checkIfExists : true)
results = params.outputdir

// Return personalized error when one of the files is missing
workflow.onError {
  filename = workflow.errorReport
  println ""
  println "---------------- Error Message -----------------"
  println " The file ${filename} is missing"  
  println " Please verify that the file exists" 
  println "------------------------------------------------"
  println ""
}

process singleFaa {
  publishDir "$results/genomes", mode: "copy"

  input: 
  file genome_dir from genomes

  output:
  file 'all_genomes.faa' into all_genomes_hmmsearch_ch
  file 'all_genomes.faa' into all_genomes_classify_ch
  file 'processed_genomes.txt' into processed_genomes_ch

  shell:
  """
  for f in \$(find ${genome_dir}/ -name '*.faa.gz'); do
    a=\$(basename \$f | sed 's/\\..*//')
    echo "\$a: \$f" >> processed_genomes.txt
    gunzip -c \$f | sed "/^>/s/\$/ [\$a]/" >> all_genomes.faa
  done
  """
}

process hmmSearch {
  publishDir "$results/hmmsearch", mode: 'copy'
  cpus params.max_cpus

  input:
  file genome  from all_genomes_hmmsearch_ch
  file hmms from hmm_files.collect()
  
  output:
  file ("*.tblout") into tblout_ch
  file ("*.domtblout") into domtblout_ch
  file ("*.hmmout.gz") into hmmout_ch

  shell:
  """
  for h in *.hmm; do
    hmmsearch --cpu ${task.cpus} --tblout \$(basename \$h .hmm).tblout --domtblout \$(basename \$h .hmm).domtblout \$h $genome | gzip -c > \$(basename \$h .hmm).hmmout.gz
  done
  """
}

process getMetadata {
  publishDir "$results/metadata", mode: 'copy'

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
  publishDir "$results/classification", mode: 'copy'

  input:
  val hmm_mincov from hmm_mincov
  val dbsource from dbsource
  file profiles_hierarchy from profiles_hierarchy
  file "gtdb_metadata.tsv" from gtdbmetadata_ch
  file tblouts from tblout_ch
  file domtblouts from domtblout_ch
  file genomes from all_genomes_classify_ch

  output:
  file "gtdb.pf.db" into gtdb_pf_db_ch
  file "gtdb.tsv.gz" into gtdb_tsv_ch
  file "gtdb.pf-classify.warnings.txt" into gtdb_classify_warnings_ch

  shell:
  """
  pf-classify.r --hmm_mincov=${hmm_mincov} --dbsource=${dbsource} --gtdbmetadata=gtdb_metadata.tsv --profilehierarchies=$profiles_hierarchy --singletable=gtdb.tsv.gz --seqfaa=${genomes} --sqlitedb=gtdb.pf.db  *.tblout *.domtblout > gtdb.pf-classify.warnings.txt 2>&1
  """
}

process db2Feather {
  publishDir "$results/feather", mode: 'copy'

  input:
  file "gtdb.pf.db" from gtdb_pf_db_ch

  output:
  file "gtdb.pf-db2feather.warnings.txt" into gtdb_db2feather_warnings_ch
  file "*.feather" into feather_files

  shell:
  """
  pf-db2feather.r --gtdb --prefix=pfitmap gtdb.pf.db > gtdb.pf-db2feather.warnings.txt  2>&1
  """
}

workflow.onComplete {
  println
  println "*** Workflow finished, status: ${ workflow.success ? 'OK' : 'failed' }. Make sure you check $results/classification/gtdb.pf-classify.warnings.txt. ***"
  println
}
