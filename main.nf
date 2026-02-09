nextflow.enable.dsl=2

workflow {
  // Sicherheits-Check für EPI2ME
  if( !params.reads || !params.db_fasta ) { return }

  Channel
    .fromPath("${params.reads}/*/*.fastq.gz", checkIfExists: true)
    .map { f -> tuple(f.parent.name, f) }
    .groupTuple()
    .map { sample, files -> tuple(sample, files.sort()) }
    .set { ch_samples }

  merged   = MERGE_FASTQ(ch_samples)
  qc_raw   = QC_SEQKIT(merged)

  // Trimmen nur wenn Primer angegeben sind
  trimmed  = (params.primer_fwd?.trim() && params.primer_rev?.trim()) ? TRIM_CUTADAPT(merged) : merged
  filtered = FILTER_NANOFILT(trimmed)

  fasta     = FASTQ_TO_FASTA(filtered)
  clustered = CLUSTER_VSEARCH(fasta)
  kept      = FILTER_CLUSTERS(clustered)

  // Datenbank vorbereiten
  ref_fa   = Channel.value( file(params.db_fasta) )
  dbdir_ch = MAKEBLASTDB(ref_fa)

  kept_with_db = kept.combine(dbdir_ch)
  blasted      = BLASTN(kept_with_db)

  tax = JOIN_COUNTS_BLAST(blasted)
  
  // Zusammenfassen
  tax_tables_list = tax.map { sample, taxfile -> taxfile }.collect() 
  summary = AGGREGATE_RESULTS(tax_tables_list)

  REPORT_HTML(summary)
}

process MERGE_FASTQ {
  tag "$sample"
  publishDir "${params.out_dir}/per_sample/${sample}/reads", mode: 'copy'
  input: tuple val(sample), path(reads)
  output: tuple val(sample), path("${sample}.fastq.gz")
  script: "cat ${reads.join(' ')} > ${sample}.fastq.gz"
}

process QC_SEQKIT {
  tag "$sample"
  publishDir "${params.out_dir}/per_sample/${sample}/qc", mode: 'copy'
  input: tuple val(sample), path(fq)
  output: tuple val(sample), path("${sample}.seqkit.stats.tsv")
  script: "seqkit stats -a -T ${fq} > ${sample}.seqkit.stats.tsv"
}

process TRIM_CUTADAPT {
  tag "$sample"
  publishDir "${params.out_dir}/per_sample/${sample}/reads", mode: 'copy'
  input: tuple val(sample), path(fq)
  output: tuple val(sample), path("${sample}.trimmed.fastq.gz")
  script:
  def discard = params.require_primers ? "--discard-untrimmed" : ""
  """
  cutadapt --revcomp -e ${params.trim_error_rate} \
    -g ${params.primer_fwd} -a ${params.primer_rev} \
    ${discard} -o ${sample}.trimmed.fastq.gz ${fq}
  """
}

process FILTER_NANOFILT {
  tag "$sample"
  publishDir "${params.out_dir}/per_sample/${sample}/reads", mode: 'copy'
  input: tuple val(sample), path(fq)
  output: tuple val(sample), path("${sample}.filtered.fastq.gz")
  script:
  """
  NanoFilt -q ${params.min_q} -l ${params.min_len} --maxlength ${params.max_len} ${fq} | gzip > ${sample}.filtered.fastq.gz
  """
}

process FASTQ_TO_FASTA {
  tag "$sample"
  publishDir "${params.out_dir}/per_sample/${sample}/clusters", mode: 'copy'
  input: tuple val(sample), path(fq)
  output: tuple val(sample), path("${sample}.fasta")
  script: "seqkit fq2fa ${fq} > ${sample}.fasta"
}

process CLUSTER_VSEARCH {
  tag "$sample"
  publishDir "${params.out_dir}/per_sample/${sample}/clusters", mode: 'copy'
  input: tuple val(sample), path(fa)
  output: tuple val(sample), path("${sample}.centroids.fasta"), path("${sample}.clusters.uc")
  script:
  """
  vsearch --cluster_fast ${fa} --id ${params.cluster_id} --strand both \
    --threads ${task.cpus} --uc ${sample}.clusters.uc --centroids ${sample}.centroids.fasta
  """
}

process FILTER_CLUSTERS {
  tag "$sample"
  publishDir "${params.out_dir}/per_sample/${sample}/clusters", mode: 'copy'
  input: tuple val(sample), path(centroids), path(uc)
  output: tuple val(sample), path("${sample}.cluster_counts.tsv"), path("${sample}.centroids.kept.fasta")
  script:
  """
  awk -F'\\t' '\$1=="S"{cl=\$2; id=\$9; c[cl]=id; n[id]=1} \$1=="H"{n[c[\$2]]++} END{for(i in n) print i,n[i]}' ${uc} | sort -k2,2nr > ${sample}.cluster_counts.tsv
  awk -v m=${params.min_cluster_reads} '\$2>=m{print \$1}' ${sample}.cluster_counts.tsv > keep.txt
  if [ -s keep.txt ]; then seqkit grep -f keep.txt ${centroids} > ${sample}.centroids.kept.fasta; else touch ${sample}.centroids.kept.fasta; fi
  """
}

process MAKEBLASTDB {
  publishDir "${params.out_dir}/refdb", mode: 'copy'
  input: path(db_fasta)
  output: path("blastdb")
  script:
  """
  mkdir -p blastdb
  cp ${db_fasta} blastdb/db.fasta
  makeblastdb -in blastdb/db.fasta -dbtype nucl -out blastdb/fusarium_tef1
  """
}

process BLASTN {
  tag "$sample"
  publishDir "${params.out_dir}/per_sample/${sample}/taxonomy", mode: 'copy'
  input: tuple val(sample), path(counts), path(centroids_fa), path(dbdir)
  output: tuple val(sample), path(counts), path("${sample}.blast.tsv")
  script:
  """
  if [ ! -s ${centroids_fa} ]; then touch ${sample}.blast.tsv; exit 0; fi
  blastn -query ${centroids_fa} -db ${dbdir}/fusarium_tef1 -max_target_seqs ${params.blast_topn} \
    -num_threads ${task.cpus} -outfmt "6 qseqid sseqid pident length qlen bitscore evalue qcovs stitle" > ${sample}.blast.tsv
  """
}

process JOIN_COUNTS_BLAST {
  tag "$sample"
  publishDir "${params.out_dir}/per_sample/${sample}/taxonomy", mode: 'copy'
  input: tuple val(sample), path(counts), path(blast)
  output: tuple val(sample), path("${sample}.taxonomy.tsv")
  script:
  """
  python -c "import csv; b={r[0]:r for r in csv.reader(open('${blast}'),delimiter='\\t')} if '${blast}'!='None' and '${blast}'!='' and '${blast}'!=' ' else {}; print('cluster_id\\tcount\\tbest_hit\\tpident\\tstitle'); [print('\\t'.join([r[0],r[1]] + (b[r[0]][1:3]+[b[r[0]][8]] if r[0] in b else ['NA','NA','NA']))) for r in csv.reader(open('${counts}'),delimiter='\\t')]" > ${sample}.taxonomy.tsv
  """
}

process AGGREGATE_RESULTS {
  publishDir "${params.out_dir}/summary", mode: 'copy'
  input: path(tables)
  output: tuple path("all_samples.long.tsv"), path("abundance_matrix.tsv")
  script:
  """
  echo "sample\\tcluster_id\\tcount\\tbest_hit\\tpident\\tstitle" > all_samples.long.tsv
  for f in ${tables}; do s=\${f%%.taxonomy.tsv}; awk -v s=\$s 'NR>1{print s"\\t"\$0}' \$f >> all_samples.long.tsv; done
  
  python -c "import pandas as pd; df=pd.read_csv('all_samples.long.tsv', sep='\\t'); m=df.pivot_table(index='sample', columns='stitle', values='count', fill_value=0); m.to_csv('abundance_matrix.tsv', sep='\\t')"
  """
}

process REPORT_HTML {
  publishDir "${params.out_dir}", mode: 'copy'
  input: tuple path(l), path(m)
  output: path("report.html")
  script: "echo '<html><body><h1>Report fertig</h1></body></html>' > report.html"
}
