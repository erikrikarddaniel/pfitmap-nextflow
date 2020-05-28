/**
 * main.nf: Nextflow workflow for TARA searches.
 *
 * The workflow starts from a set of downloaded paired fastq files and:
 *
 *   Cuts adaptors (trim_galore)
 *   Quality trims (sickle)
 *   Runs FastQC/MultiQC on the cleaned reads
 *   Runs Diamond against:
 *     Andy's Scripps Pier bloom transcriptome (bloom_transcriptome)
 *     Selected nitrogenase gene products (genes)
 *     MMETSP proteins (mmetsp)
 *   Maps reads to a selection of genomes (genomes)
 *   Runs LAST against:
 *     bloom_transcriptome
 *     genes
 *     mmetsp
 *   Summarises by:
 *     producing an html report
 *
 * daniel.lundin@lnu.se
 */

params.help = false
params.inputdir = "batch/"
params.outputdir = "results/"
params.genomeindex = ""
params.geneslastindex = ""
params.genesdiamondindex = ""
params.bloom_transcriptomelastindex = ""
params.bloom_transcriptomediamondindex = ""
params.mmetsplastindex = ""
params.mmetspdiamondindex = ""
params.project = "snic2019-8-203"
params.samplefeather = ""       // Path to a table in feather format containing sample information à la "simple_sample_table.feather"
params.genesfeather = ""        // Path to a table in feather format containing gene information à la "reference/genes/genes_seq_summary.feather" (see Makefile)
params.mmetspfeather = ""       // Path to a table in feather format containing gene information à la "reference/mmetsp/mmetsp_seq_summary.feather" (see Makefile)

params.max_cpus = 1
params.last_cpus = 1
params.max_time = "240.h"

params.run_last2genes = false
params.run_last2mmetsp = false
params.run_last2bloom_transcriptome = false

//println "--> run_last: ${params.run_last} <--"

// Result files will contain a timestamp to differentiate between different runs
now = new Date()
timestamp = now.format("yyyyMMdd_HHmmss", TimeZone.getTimeZone('CET'))
println "MAX_TIME: ${params.max_time}"

def helpMessage() {
  log.info """

  Usage:

  The typical command for running the pipeline is as follows:

  nextflow run main.nf --inputdir batch --outputdir results --genomeindex PATH

  Mandatory arguments:
  --genomeindex PATH			Path to a bowtie2index for genomes to map to.
  --geneslastindex PATH			Path to LAST database with genes.
  --genesdiamondindex PATH		Path to Diamond database with genes.
  --bloom_transcriptomelastindex PATH	Path to LAST database with bloom_transcriptome.
  --bloom_transcriptomediamondindex PATH	Path to Diamond database with bloom_transcriptome.
  --mmetsplastindex PATH		Path to LAST database with mmetsp.
  --mmetspdiamondindex PATH		Path to Diamond database with mmetsp.
  --samplefeather PATH          	Path to a feather file with sample information (à la "simple_sample_table.feather")

  Mandatory arguments with defaults:
  --inputdir PATH			Directory with raw sequence read files in gzipped fastq format
  --outputdir PATH			Directory in which output files will be placed

  Non-mandatory arguments:
  --max_cpus				Maximum number of cpu cores, default 1
  --last_cpus				Number of cpu cores to allocate to LAST, default 1
  --max_time				Maximum time per process, default 10 days
  --run_last2genes			Run LAST on the genes database, default false
  --run_last2mmetsp			Run LAST on the mmetsp database, default false
  --run_last2bloom_transcriptome	Run LAST on the bloom_transcriptome database, default false
  """.stripIndent()
}

// Show help message
if (params.help) {
  helpMessage()
  exit 0
}

// Handle index paths
if ( params.genomeindex == "" ) {
  exit 1, "You need to specify the path to a bowtie2 index for genomes to map to (--genomeindex)"
}
genomeindex_dir = file(params.genomeindex).getParent()
genomeindex_name = file(params.genomeindex).getName()

if ( params.geneslastindex == "" ) {
  exit 1, "You need to specify the path to a LAST index for genes to align to (--geneslastindex)"
}
geneslastindex_dir = file(params.geneslastindex).getParent()
geneslastindex_name = file(params.geneslastindex).getName()

if ( params.genesdiamondindex == "" ) {
  exit 1, "You need to specify the path to a Diamond index for genes to align to (--genesdiamondindex)"
}
genesdiamondindex_dir = file(params.genesdiamondindex).getParent()
genesdiamondindex_name = file(params.genesdiamondindex).getName()

if ( params.bloom_transcriptomelastindex == "" ) {
  exit 1, "You need to specify the path to a LAST index for bloom_transcriptome to align to (--bloom_transcriptomelastindex)"
}
bloom_transcriptomelastindex_dir = file(params.bloom_transcriptomelastindex).getParent()
bloom_transcriptomelastindex_name = file(params.bloom_transcriptomelastindex).getName()

if ( params.bloom_transcriptomediamondindex == "" ) {
  exit 1, "You need to specify the path to a Diamond index for bloom_transcriptome to align to (--bloom_transcriptomediamondindex)"
}
bloom_transcriptomediamondindex_dir = file(params.bloom_transcriptomediamondindex).getParent()
bloom_transcriptomediamondindex_name = file(params.bloom_transcriptomediamondindex).getName()

if ( params.mmetsplastindex == "" ) {
  exit 1, "You need to specify the path to a LAST index for mmetsp to align to (--mmetsplastindex)"
}
mmetsplastindex_dir = file(params.mmetsplastindex).getParent()
mmetsplastindex_name = file(params.mmetsplastindex).getName()

if ( params.mmetspdiamondindex == "" ) {
  exit 1, "You need to specify the path to a Diamond index for mmetsp to align to (--mmetspdiamondindex)"
}
mmetspdiamondindex_dir = file(params.mmetspdiamondindex).getParent()
mmetspdiamondindex_name = file(params.mmetspdiamondindex).getName()

// Summary input file params
if ( params.samplefeather == "" ) {
  exit 1, "You need to specify the path to a sample table for the summary html (--samplefeather)"
}
if ( params.genesfeather == "" ) {
  exit 1, "You need to specify the path to a gene info table for the summary html (--genesfeather)"
}
if ( params.mmetspfeather == "" ) {
  exit 1, "You need to specify the path to a MMETSP info table for the summary html (--mmetspfeather)"
}

lastalopts = "" // Place holder

// Creating a channel for input fastq files: kicks off the processing.
inputfiles = Channel.fromFilePairs("$params.inputdir/*{1,2}.fastq.gz")

// And for files for the html summary
rmarkdown 	= Channel.fromPath("assets/summary.Rmd")
samplefeather 	= Channel.fromPath(params.samplefeather)
genesfeather    = Channel.fromPath(params.genesfeather)
mmetspfeather   = Channel.fromPath(params.mmetspfeather)

/**
 * 1. Trim and remove adapters.
 */
process trim_galore {
  cpus 1
  time params.max_time

  input:
  set name, file(reads) from inputfiles

  output:
  set val(name), file("*fq.gz") into (
    trimgalore_reads_map2genomes,
    trimgalore_reads_last2genes, trimgalore_reads_diamond2genes,
    trimgalore_reads_last2bloom_transcriptome, trimgalore_reads_diamond2bloom_transcriptome,
    trimgalore_reads_last2mmetsp, trimgalore_reads_diamond2mmetsp
  )
  set val(name), file("*.trim_galore.log") into ch_trimming_logs
  set val(name), file("*_trimming_report.txt") into ch_trimming_logs_rcounts
  file "*_fastqc.{zip,html}" into trimgalore_fastqc_reports

  """
  trim_galore --paired --fastqc --gzip --quality 20 $reads 2>&1 > ${name}.trim_galore.log
  """
}

/**
 * Count the number of reads in each sample. The next process collects that
 * into a common file for all samples.
 */
process sample_count {
  cpus 1
  time params.max_time
  executor 'local'

  input:
  set name, file(trimlogs) from ch_trimming_logs_rcounts

  output:
  file("*.samplecounts.tsv") into ch_sample_counts

  """
  echo \"${name}	\$(grep 'Reads written' ${trimlogs[0]} | sed 's/.*: *\\([0-9,]*\\) (.*/\\1/' | sed 's/,//g')\" > ${name}.samplecounts.tsv
  """
  //echo \"${name}	${trimlogs}:`grep 'Reads written' ${trimlogs[0]}`\" > ${name}.samplecounts.tsv
}

process sample_counts {
  cpus 1
  time params.max_time
  executor 'local'

  publishDir "${params.outputdir}", mode: "copy", pattern: "*.tsv"

  input:
  file(counts) from ch_sample_counts.collect()

  output:
  file("samplecounts.${timestamp}.tsv") into ch_sample_counts_out
  file("samplecounts.${timestamp}.tsv") into ch_sample_counts_rmd
  file("samplecounts.md5") into ch_sample_counts_out_md5

  """
  cat ${counts} > samplecounts.${timestamp}.tsv
  md5sum samplecounts.${timestamp}.tsv > samplecounts.md5
  """
}

/**
 * MultiQC on output from trim_galore.
 */
process multiqc {
  cpus 1
  time params.max_time

  publishDir "${params.outputdir}/multiqc", mode: "copy", pattern: "*.html"

  input:
  file fastqc_reports from trimgalore_fastqc_reports.collect()

  output:
  file "multiqc_${timestamp}_multiqc_report.html" into multiqc_report
  file "multiqc.md5" into multiqc_report_md5
  file "*_data"

  """
  multiqc -i "multiqc_${timestamp}" .
  md5sum multiqc_${timestamp}_multiqc_report.html | sed 's: \\+:&multiqc/:' > multiqc.md5
  """
}

/**
 * Map pairs to genome index.
 */
process map2genomes {
  cpus params.max_cpus
  time params.max_time

  input:
  set name, file(reads) from trimgalore_reads_map2genomes

  output:
  file "*.bam" into (
    map2genomes_result, 	// Consumed by tar process
    ch_map2genomes_idxstats
  )
  file "*.log" into map2genomes_log

  """
  bowtie2 --threads ${task.cpus} -x ${genomeindex_dir}/${genomeindex_name} -1 ${reads[0]} -2 ${reads[1]} 2>${name}.genomes.bowtie2.log | samtools view -Sb | samtools sort > ${name}.genomes.sorted.bam
  """
}

/**
 * Write idxstats files for bams.
 */
process idxstats {
  cpus 1
  time params.max_time

  input:
  file bam from ch_map2genomes_idxstats

  output:
  file("*.bam.bai") into ch_bam_bai
  file("*.idxstats.tsv.gz") into ch_idxstats_results
  file("*.idxstats.tsv.gz") into ch_idxstats_results_rmd

  """
  samtools index ${bam} 
  samtools idxstats $bam | grep -Pv "\t0\t0\$" | gzip -c > \$(basename ${bam}).idxstats.tsv.gz
  """
}

/**
 * Tar up bam files from map2genomes.
 */
process tar_map2genomes {
  cpus 1
  time params.max_time

  publishDir "${params.outputdir}/nucleotides/genomes", mode: "copy"

  input:
  file bam_files from map2genomes_result.collect()
  file log_files from map2genomes_log.collect()
  file bam_bai_files from ch_bam_bai.collect()
  file idxstats from ch_idxstats_results.collect()

  output:
  file "map2genomes.${timestamp}.tar" into tar_map2genomes_result
  file "map2genomes_tar.md5" into tar_map2genomes_result_md5

  """
  tar --dereference -cf map2genomes.${timestamp}.tar ${bam_files} ${log_files} ${bam_bai_files} ${idxstats}
  md5sum map2genomes.${timestamp}.tar | sed 's: \\+:&nucleotides/genomes/:' > map2genomes_tar.md5
  """
}

/**
 * LAST pairs to genes.
 */
process last2genes {
  cpus params.last_cpus
  time params.max_time

  input:
  set name, file(reads) from trimgalore_reads_last2genes

  output:
  set val(name), file("*.maf.gz") into last2genes_result1 // Consumed by last2genestsv process
  file "*.maf.gz" into last2genes_result2 		  // Consumed by tar process

  when:
  params.run_last2genes

  """
  unpigz -c -p ${task.cpus} ${reads} | lastal -P ${task.cpus} -Q0 ${lastalopts} ${geneslastindex_dir}/${geneslastindex_name}.lastdb | pigz -c -p ${task.cpus} > ${name}.genes.maf.gz
  """
}

process last2genestsv {
  cpus 1
  time params.max_time

  input:
  set name, file(maf_file) from last2genes_result1

  output:
  file "*.tsv.gz" into last2genestsv_result1	// tar_last2genes
  file "*.tsv.gz" into last2genestsv_result_rmd	// Rmd

  """
  unpigz -c -p ${task.cpus} ${maf_file} | maf-convert blasttab | pigz -c -p ${task.cpus} > ${name}.genes.last.tsv.gz
  """
}

/**
 * Tar up bam files from last2genes.
 */
process tar_last2genes {
  cpus 1
  time params.max_time

  publishDir "${params.outputdir}/nucleotides/genes", mode: "copy", pattern: "*.tar"

  input:
  file maf_files from last2genes_result2.collect()
  file tsv_files from last2genestsv_result1.collect()

  output:
  file "last2genes.${timestamp}.tar" into tar_last2genes_result
  file "last2genes.md5" into tar_last2genes_result_md5

  """
  tar --dereference -cf last2genes.${timestamp}.tar ${maf_files} ${tsv_files}
  md5sum last2genes.${timestamp}.tar | sed 's: \\+:&nucleotides/genes/:' > last2genes.md5
  """
}

/**
 * Diamond pairs to genes.
 */
process diamond2genes {
  cpus params.last_cpus
  time params.max_time

  input:
  set name, file(reads) from trimgalore_reads_diamond2genes

  output:
  set val(name), file("*.daa") into diamond2genes_result1	// diamond2genestsv process
  file "*.daa" into diamond2genes_result2 	// Consumed by tar process
  file "*.daa.log.gz" into diamond2genes_log

  """
  diamond blastx --outfmt 100 --threads ${task.cpus} --query ${reads[0]} --db ${genesdiamondindex_dir}/${genesdiamondindex_name}.dmnd --out ${name}_1.genes.daa 2>&1 | gzip -c > ${name}_1.genes.daa.log.gz
  diamond blastx --outfmt 100 --threads ${task.cpus} --query ${reads[1]} --db ${genesdiamondindex_dir}/${genesdiamondindex_name}.dmnd --out ${name}_2.genes.daa 2>&1 | gzip -c > ${name}_2.genes.daa.log.gz
  """
}

process diamond2genestsv {
  cpus params.last_cpus
  time params.max_time

  input:
  set name, file(daa_files) from diamond2genes_result1

  output:
  file "*.tsv.gz" into diamond2genestsv_result0	// md5_diamond2genes
  file "*.tsv.gz" into diamond2genestsv_result1	// tar_diamond2genes
  file "*.tsv.gz" into diamond2genestsv_result_rmd

  """
  diamond view -p ${task.cpus} --daa  ${daa_files[0]} -f 6 | pigz -c -p ${task.cpus} > ${name}_1.genes.dmnd.tsv.gz
  diamond view -p ${task.cpus} --daa  ${daa_files[1]} -f 6 | pigz -c -p ${task.cpus} > ${name}_2.genes.dmnd.tsv.gz
  """
}

/**
 * Tar up bam files from diamond2genes.
 */
process tar_diamond2genes {
  cpus 1
  time params.max_time

  publishDir "${params.outputdir}/proteins/genes", mode: "copy", pattern: "*.tar"

  input:
  file daa_files from diamond2genes_result2.collect()
  file tsv_files from diamond2genestsv_result1.collect()
  file log_files from diamond2genes_log.collect()

  output:
  file "diamond2genes.${timestamp}.tar" into tar_diamond2genes_result
  file "diamond2genes.md5" into tar_diamond2genes_result_md5

  """
  tar --dereference -cf diamond2genes.${timestamp}.tar ${daa_files} ${tsv_files} ${log_files}
  md5sum diamond2genes.${timestamp}.tar | sed 's: \\+:&proteins/genes/:' > diamond2genes.md5
  """
}

/**
 * LAST pairs to bloom_transcriptome.
 */
process last2bloom_transcriptome {
  cpus params.last_cpus
  time params.max_time

  input:
  set name, file(reads) from trimgalore_reads_last2bloom_transcriptome

  output:
  set val(name), file("*.maf.gz") into last2bloom_transcriptome_result1  // Consumed by last2bloom_transcriptometsv process
  file "*.maf.gz" into last2bloom_transcriptome_result2 	// Consumed by tar process

  when:
  params.run_last2bloom_transcriptome

  """
  unpigz -c -p ${task.cpus} ${reads} | lastal -P ${task.cpus} -Q0 ${lastalopts} ${bloom_transcriptomelastindex_dir}/${bloom_transcriptomelastindex_name}.lastdb | pigz -c -p ${task.cpus} > ${name}.bloom_transcriptome.maf.gz
  """
}

process last2bloom_transcriptometsv {
  cpus 1
  time params.max_time

  input:
  set name, file(maf_file) from last2bloom_transcriptome_result1

  output:
  file "*.tsv.gz" into last2bloom_transcriptometsv_result1	// tar_last2bloom_transcriptome
  file "*.tsv.gz" into last2bloom_transcriptometsv_result_rmd

  """
  unpigz -c -p ${task.cpus} ${maf_file} | maf-convert blasttab | pigz -c -p ${task.cpus} > ${name}.bloom_transcriptome.last.tsv.gz
  """
}

/**
 * Tar up bam files from last2bloom_transcriptome.
 */
process tar_last2bloom_transcriptome {
  cpus 1
  time params.max_time

  publishDir "${params.outputdir}/nucleotides/bloom_transcriptome", mode: "copy", pattern: "*.tar"

  input:
  file maf_files from last2bloom_transcriptome_result2.collect()
  file tsv_files from last2bloom_transcriptometsv_result1.collect()

  output:
  file "last2bloom_transcriptome.${timestamp}.tar" into tar_last2bloom_transcriptome_result
  file "last2bloom_transcriptome.md5" into tar_last2bloom_transcriptome_result_md5

  """
  tar --dereference -cf last2bloom_transcriptome.${timestamp}.tar ${maf_files} ${tsv_files}
  md5sum last2bloom_transcriptome.${timestamp}.tar | sed 's: \\+:&nucleotides/bloom_transcriptome/:' > last2bloom_transcriptome.md5
  """
}

/**
 * LAST pairs to bloom_transcriptome.
 */
process diamond2bloom_transcriptome {
  cpus params.last_cpus
  time params.max_time

  input:
  set name, file(reads) from trimgalore_reads_diamond2bloom_transcriptome

  output:
  set val(name), file("*.daa") into diamond2bloom_transcriptome_result1	// diamond2bloom_transcriptometsv process
  file "*.daa" into diamond2bloom_transcriptome_result2 	// Consumed by tar process
  file "*.daa.log.gz" into diamond2bloom_transcriptome_log

  """
  diamond blastx --outfmt 100 --threads ${task.cpus} --query ${reads[0]} --db ${bloom_transcriptomediamondindex_dir}/${bloom_transcriptomediamondindex_name}.dmnd --out ${name}_1.bloom_transcriptome.daa 2>&1 | gzip -c > ${name}_1.bloom_transcriptome.daa.log.gz
  diamond blastx --outfmt 100 --threads ${task.cpus} --query ${reads[1]} --db ${bloom_transcriptomediamondindex_dir}/${bloom_transcriptomediamondindex_name}.dmnd --out ${name}_2.bloom_transcriptome.daa 2>&1 | gzip -c > ${name}_2.bloom_transcriptome.daa.log.gz
  """
}

process diamond2bloom_transcriptometsv {
  cpus params.last_cpus
  time params.max_time

  input:
  set name, file(daa_files) from diamond2bloom_transcriptome_result1

  output:
  file "*.tsv.gz" into diamond2bloom_transcriptometsv_result1	// tar_diamond2bloom_transcriptome
  file "*.tsv.gz" into diamond2bloom_transcriptometsv_result_rmd

  """
  diamond view -p ${task.cpus} --daa  ${daa_files[0]} -f 6 | pigz -c -p ${task.cpus} > ${name}_1.bloom_transcriptome.dmnd.tsv.gz
  diamond view -p ${task.cpus} --daa  ${daa_files[1]} -f 6 | pigz -c -p ${task.cpus} > ${name}_2.bloom_transcriptome.dmnd.tsv.gz
  """
}

/**
 * Tar up bam files from diamond2bloom_transcriptome.
 */
process tar_diamond2bloom_transcriptome {
  cpus 1
  time params.max_time

  publishDir "${params.outputdir}/proteins/bloom_transcriptome", mode: "copy", pattern: "*.tar"

  input:
  file daa_files from diamond2bloom_transcriptome_result2.collect()
  file tsv_files from diamond2bloom_transcriptometsv_result1.collect()
  file log_files from diamond2bloom_transcriptome_log.collect()

  output:
  file "diamond2bloom_transcriptome.${timestamp}.tar" into tar_diamond2bloom_transcriptome_result
  file "diamond2bloom_transcriptome.md5" into tar_diamond2bloom_transcriptome_result_md5

  """
  tar --dereference -cf diamond2bloom_transcriptome.${timestamp}.tar ${daa_files} ${tsv_files} ${log_files}
  md5sum diamond2bloom_transcriptome.${timestamp}.tar | sed 's: \\+:&proteins/bloom_transcriptome/:' > diamond2bloom_transcriptome.md5
  """
}

/**
 * LAST pairs to mmetsp.
 */
process last2mmetsp {
  cpus params.last_cpus
  time params.max_time

  input:
  set name, file(reads) from trimgalore_reads_last2mmetsp

  output:
  set val(name), file("*.maf.gz") into last2mmetsp_result1  // Consumed by last2mmetsptsv process
  file "*.maf.gz" into last2mmetsp_result2 	// Consumed by tar process

  when:
  params.run_last2mmetsp

  """
  unpigz -c -p ${task.cpus} ${reads} | lastal -P ${task.cpus} -Q0 ${lastalopts} ${mmetsplastindex_dir}/${mmetsplastindex_name}.lastdb | pigz -c -p ${task.cpus} > ${name}.mmetsp.maf.gz
  """
}

process last2mmetsptsv {
  cpus 1
  time params.max_time

  input:
  set name, file(maf_file) from last2mmetsp_result1

  output:
  file "*.tsv.gz" into last2mmetsptsv_result1	// tar_last2mmetsp
  file "*.tsv.gz" into last2mmetsptsv_result_rmd

  """
  unpigz -c -p ${task.cpus} ${maf_file} | maf-convert blasttab | pigz -c -p ${task.cpus} > ${name}.mmetsp.last.tsv.gz
  """
}

/**
 * Tar up bam files from last2mmetsp.
 */
process tar_last2mmetsp {
  cpus 1
  time params.max_time

  publishDir "${params.outputdir}/nucleotides/mmetsp", mode: "copy", pattern: "*.tar"

  input:
  file maf_files from last2mmetsp_result2.collect()
  file tsv_files from last2mmetsptsv_result1.collect()

  output:
  file "last2mmetsp.${timestamp}.tar" into tar_last2mmetsp_result
  file "last2mmetsp.md5" into tar_last2mmetsp_result_md5

  """
  tar --dereference -cf last2mmetsp.${timestamp}.tar ${maf_files} ${tsv_files}
  md5sum last2mmetsp.${timestamp}.tar | sed 's: \\+:&nucleotides/mmetsp/:' > last2mmetsp.md5
  """
}

/**
 * Diamond pairs to mmetsp.
 */
process diamond2mmetsp {
  cpus params.last_cpus
  time params.max_time

  input:
  set name, file(reads) from trimgalore_reads_diamond2mmetsp

  output:
  set val(name), file("*.daa") into diamond2mmetsp_result1	// diamond2mmetsptsv process
  file "*.daa" into diamond2mmetsp_result2 	// Consumed by tar process
  file "*.daa.log.gz" into diamond2mmetsp_log

  """
  diamond blastx --outfmt 100 --threads ${task.cpus} --query ${reads[0]} --db ${mmetspdiamondindex_dir}/${mmetspdiamondindex_name}.dmnd --out ${name}_1.mmetsp.daa 2>&1 | gzip -c > ${name}_1.mmetsp.daa.log.gz
  diamond blastx --outfmt 100 --threads ${task.cpus} --query ${reads[1]} --db ${mmetspdiamondindex_dir}/${mmetspdiamondindex_name}.dmnd --out ${name}_2.mmetsp.daa 2>&1 | gzip -c > ${name}_2.mmetsp.daa.log.gz
  """
}

process diamond2mmetsptsv {
  cpus params.last_cpus
  time params.max_time

  input:
  set name, file(daa_files) from diamond2mmetsp_result1

  output:
  file "*.tsv.gz" into diamond2mmetsptsv_result1	// tar_diamond2mmetsp
  file "*.tsv.gz" into diamond2mmetsptsv_result_rmd

  """
  diamond view -p ${task.cpus} --daa  ${daa_files[0]} -f 6 | pigz -c -p ${task.cpus} > ${name}_1.mmetsp.dmnd.tsv.gz
  diamond view -p ${task.cpus} --daa  ${daa_files[1]} -f 6 | pigz -c -p ${task.cpus} > ${name}_2.mmetsp.dmnd.tsv.gz
  """
}

/**
 * Tar up bam files from diamond2mmetsp.
 */
process tar_diamond2mmetsp {
  cpus 1
  time params.max_time

  publishDir "${params.outputdir}/proteins/mmetsp", mode: "copy", pattern: "*.tar"

  input:
  file daa_files from diamond2mmetsp_result2.collect()
  file tsv_files from diamond2mmetsptsv_result1.collect()
  file log_files from diamond2mmetsp_log.collect()

  output:
  file "diamond2mmetsp.${timestamp}.tar" into tar_diamond2mmetsp_result
  file "diamond2mmetsp.md5" into tar_diamond2mmetsp_result_md5

  """
  tar --dereference -cf diamond2mmetsp.${timestamp}.tar ${daa_files} ${tsv_files} ${log_files}
  md5sum diamond2mmetsp.${timestamp}.tar | sed 's: \\+:&proteins/mmetsp/:' > diamond2mmetsp.md5
  """
}

/**
 * Generate some summary stats by knitting an Rmarkdown doc.
 */
process rmarkdown {
  //cpus params.max_cpus
  cpus 4
  time params.max_time
  
  publishDir "${params.outputdir}/summary/", mode: "copy"

  input:
  file rmd	from rmarkdown
  file scounts	from ch_sample_counts_rmd
  file idxstats	from ch_idxstats_results_rmd.collect()
  file lgenes	from last2genestsv_result_rmd.ifEmpty(['empty',[]]).collect() 
  file dgenes	from diamond2genestsv_result_rmd.collect()
  file lbloom	from last2bloom_transcriptometsv_result_rmd.ifEmpty(['empty',[]]).collect() 
  file dbloom	from diamond2bloom_transcriptometsv_result_rmd.collect()
  file lmmetsp	from last2mmetsptsv_result_rmd.ifEmpty(['empty',[]]).collect() 
  file dmmetsp	from diamond2mmetsptsv_result_rmd.collect()
  file samples	from samplefeather
  file genes    from genesfeather
  file mmetsp   from mmetspfeather

  output:
  file "summary.${timestamp}.html"
  file "*.${timestamp}.feather"
  file "summary.md5" into summary_md5

  """
  sed 's/__DATE__/${timestamp}/' ${rmd} | sed 's/__THREADS__/${task.cpus}/g' > summary.${timestamp}.Rmd
  Rscript --vanilla -e "rmarkdown::render('summary.${timestamp}.Rmd')"
  md5sum *.${timestamp}.feather summary.${timestamp}.html | sed 's: \\+:&summary/:' > summary.md5
  """
}

/**
 * Calculate the md5 sum file.
 */
process md5sums {
  cpus 1
  time params.max_time
  executor "local"

  publishDir "${params.outputdir}/", mode: "copy"

  input:
  file(scount)         from ch_sample_counts_out_md5
  file(multiqc)        from multiqc_report_md5
  file(map2genomes)    from tar_map2genomes_result_md5
  file(last2genes)     from tar_last2genes_result_md5.ifEmpty(['empty',[]])
  file(diamond2genes)  from tar_diamond2genes_result_md5
  file(last2bloom)     from tar_last2bloom_transcriptome_result_md5.ifEmpty(['empty',[]])
  file(diamond2bloom)  from tar_diamond2bloom_transcriptome_result_md5
  file(last2mmetsp)    from tar_last2mmetsp_result_md5.ifEmpty(['empty',[]])
  file(diamond2mmetsp) from tar_diamond2mmetsp_result_md5
  file(summary)	       from summary_md5

  output:
  file "nextflow_results.${timestamp}.md5" into ch_md5sums

  """
  cat *.md5 > nextflow_results.${timestamp}.md5
  """
}
