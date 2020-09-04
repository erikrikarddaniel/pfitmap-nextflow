#!/usr/bin/env nextflow

/**
 * main.nf: Nextflow workflow for Pfitmap/GTDB
 *
 * The GTDB genomes are expected to be downloaded and annotated.
 *
 * The workflow starts from a set of annotated genomes in the format of faa.gz files (--inputfaas) 
 * and, optionally, gff.gz files (--inputgffs) plus a set of hmm profiles (--hmms). The protein 
 * sequences will be searched with HMMER using the hmm files and subsequently classified into which 
 * profile it fits best into. The latter uses a table describing the hierarchy of hmm profiles
 * (--profiles_hierarchy; see --help).
 *
 * Requirements: 
 *   directory with faa.gz files
 *   directory with all hmm profiles to be run 
 *   file describing the hmm profile hierarchy
 *
 * Processing steps:
 *   Concatenate all faa.gz files into a single one
 *   Optionally, concatenate all gff.gz files into a single one
 *   Perform an hmmsearch of all hmm profiles on all the proteomes
 *   Download the metadata files for archaeal and bacterial genomes from gtdb latest version 
 *     repository and concatenate them into a single metadata file
 *   Classify the hits 
 *
 * ghada.nouraia@dbb.su.se daniel.lundin@dbb.su.se
 */

// Parameters
params.help                     = false
params.inputfaas                = null
params.inputgffs                = null
params.hmms                     = null
params.profiles_hierarchy       = null
params.dbsource                 = 'GTDB:GTDB:latest'
params.hmm_mincov               = 0.7
params.gtdb_arc_metadata        = null
params.gtdb_bac_metadata        = null
params.featherprefix           = 'pfitmap-gtdb'

params.max_cpus = 2
params.max_time = "240.h"

def helpMessage() {
  log.info """

  Usage:

  The typical command for running the pipeline is as follows:

  nextflow run main.nf --inputfaas path/to/genomes.faa.gzs [--inputgffs path/to/genomes.gff.gzs] --outputdir path/to/results --hmm_mincov value --dbsource GTDB:GTDB:release

  Mandatory arguments:
    --inputfaas path/to/genomes.faa.gzs  		Path of directory containing annotated genomes in the format faa.gz 
    --gtdb_bac_metadata path/to/file			Path of tsv file including the metadata for bacterial genomes
    --gtdb_arc_metadata path/to/file 			Path of tsv file including the metadata for archaeal genomes
    --hmms path/to/hmm_directory                        Path of directory with HMM profile files 
    --profiles_hierarchy path/to/file			Path of tsv file including hmm profile names and information (See README.md file for more details)		
    --hmm_mincov value					Set a value for the threshold of coverage hmm_profile/querry (default = 0.7)
    --dbsource db:db:release				Set the database source in the format db:db:release, where [db] is the name of the database and [release] mentions 
                                                          the release number/name (default = GTDB:GTDB:latest)
    --outputdir path/to/results				Path to the results directory
    --featherprefix prefix                             Prefix for generated feather files (default "pfitmap-gtdb").

  Non Mandatory parameters:
    --inputgffs path/to/genomes.gff.gzs  		Path of directory containing annotated genomes in the format gff.gz 
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

if( !params.inputfaas ) {
  error "Missing inputfaas parameter\n[Parameter error] Please specify the parameter --inputfaas\nSee more using --help" 
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
genome_faas        = Channel.fromPath(params.inputfaas, checkIfExists : true)
if ( params.inputgffs ) { 
  genome_gffs        = Channel.fromPath(params.inputgffs, checkIfExists : true) 
}
else {
  genome_gffs = Channel.empty()
}
hmm_files          = Channel.fromPath("$params.hmms/*.hmm")
profiles_hierarchy = Channel.fromPath(params.profiles_hierarchy, checkIfExists : true)
dbsource           = Channel.value(params.dbsource)
hmm_mincov         = Channel.value(params.hmm_mincov)
featherprefix      = Channel.value(params.featherprefix)
gtdb_arc_metadata  = Channel.fromPath(params.gtdb_arc_metadata, checkIfExists : true)
gtdb_bac_metadata  = Channel.fromPath(params.gtdb_bac_metadata, checkIfExists : true)
results            = params.outputdir

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

process join_faas {
  publishDir "$results/genomes", mode: "copy"

  input: 
  file genome_dir from genome_faas

  output:
  file 'all_genomes.faa' into all_genome_faas_hmmsearch_ch
  file 'all_genomes.faa' into all_genome_faas_classify_ch
  file 'processed_genome_faas.txt' into processed_genome_faas_ch

  shell:
  """
  for f in \$(find ${genome_dir}/ -name '*.faa.gz'); do
    a=\$(basename \$f | sed 's/\\..*//')
    echo "\$a: \$f" >> processed_genome_faas.txt
    gunzip -c \$f | sed "/^>/s/\$/ [\$a]/" >> all_genomes.faa
  done
  """
}

process join_gffs {
  publishDir "$results/genomes", mode: "copy"

  input: 
  file genome_dir from genome_gffs

  output:
  file 'all_genomes.gff.gz' into all_genome_gffs_ch

  shell:
  """
  for f in \$(find ${genome_dir}/ -name '*.gff.gz'); do
    a=\$(basename \$f | sed 's/\\..*//')
    gunzip -c \$f | sed "/^>/s/\$/ [\$a]/" >> all_genomes.gff
  done
  gzip all_genomes.gff
  """
}

/**
 * Run the hmmsearches.
 *
 * This is now called *once* but should be called once per hmm file. I haven't yet found
 * a way of duplicating the channel with the faa file though.
 */
process hmmsearch {
  publishDir "$results/hmmsearch", mode: 'copy'
  cpus params.max_cpus

  input:
  file genome  from all_genome_faas_hmmsearch_ch
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
    val featherprefix from featherprefix
    val dbsource from dbsource
    file profiles_hierarchy from profiles_hierarchy
    file "gtdb_metadata.tsv" from gtdbmetadata_ch
    file tblouts from tblout_ch
    file domtblouts from domtblout_ch
    file genome_faas from all_genome_faas_classify_ch

  output:
    file "gtdb.tsv.gz" into gtdb_tsv_ch
    file "gtdb.pf-classify.warnings.txt" into gtdb_classify_warnings_ch
    file "missing_genomes.txt" into gtdb_classify_missing_ch
    file "*.feather" into feather_files

  shell:
  """
  pf-classify.r --hmm_mincov=${hmm_mincov} --dbsource=${dbsource} --gtdbmetadata=gtdb_metadata.tsv --profilehierarchies=$profiles_hierarchy --singletable=gtdb.tsv.gz --seqfaa=${genome_faas} --featherprefix=${featherprefix}  --missing=missing_genomes.txt *.tblout *.domtblout > gtdb.pf-classify.warnings.txt 2>&1
  """
}

workflow.onComplete {
  println
  println "*** Workflow finished, status: ${ workflow.success ? 'OK' : 'failed' }. Make sure you check $results/classification/gtdb.pf-classify.warnings.txt. ***"
  println
}
