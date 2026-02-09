nextflow.enable.dsl=2

/*
  Haupt-Workflow: BYPASS MODUS
  Wir überspringen Cutadapt und NanoFilt komplett.
  Die Reads gehen direkt in das Clustering.
*/

workflow {
  // Check
  if( !params.reads || !params.db_fasta ) {
    log.info "Parameter fehlen."
    return
  }

  // 1. Input Kanäle
  Channel
    .fromPath("${params.reads}/*/*.fastq.gz", checkIfExists: true)
    .map { f -> tuple(f.parent.name, f) }
    .groupTuple()
    .map { sample, files -> tuple(sample, files.sort()) }
    .set { ch_samples }

  // 2. Workflow Logik (MIT BYPASS)
  merged   = MERGE_FASTQ(ch_samples)
  
  // QC läuft mit, aber wir nutzen das Ergebnis nicht für den Fluss
  qc_raw   = QC_SEQKIT(merged) 

  // --- BYPASS START ---
  // Wir ignorieren TRIM_CUTADAPT und FILTER_NANOFILT
  // Wir leiten 'merged' direkt weiter zu 'FASTQ_TO_FASTA'
  // --------------------
  
  fasta     = FASTQ_TO_FASTA(merged) // HIER IST DIE ABKÜRZUNG
  
  clustered = CLUSTER_VSEARCH(fasta)
  kept      = FILTER_CLUSTERS(clustered)
  
  ref_fa   = Channel.value( file(params.db_fasta) )
  dbdir_ch = MAKEBLASTDB(ref_fa)
  
  kept_with_db = kept.combine(dbdir_ch)
  blasted      = BLASTN(kept_with_db)
  
  tax = JOIN_COUNTS_BLAST(blasted)
  
  tax_tables_list = tax.map { sample, taxfile -> taxfile }.collect() 
  summary = AGGREGATE_RESULTS(tax_tables_list)
  
  tsv_ch = summary.map { long_tsv, mat_tsv -> [ long_tsv, mat_tsv ] }.flatten()
  csv_ch = TSV_TO_CSV(tsv_ch)
  
  csv_list = csv_ch.collect()
  csv_tuple = csv_list.map { files ->
      def m = files.collectEntries { [(it.name): it] }
      tuple(m["all_samples.long.csv"], m["abundance_matrix.csv"])
  }
  
  REPORT_HTML(csv_tuple)
}

// ------------------------------------------------------
// PROZESSE
// ------------------------------------------------------

process MERGE_FASTQ {
  tag "$sample"
  publishDir "output/per_sample/${sample}/reads", mode: 'copy'
  input: tuple val(sample), path(reads)
  output: tuple val(sample), path("${sample}.fastq.gz")
  script: "cat ${reads.join(' ')} > ${sample}.fastq.gz"
}

process QC_SEQKIT {
  tag "$sample"
  publishDir "output/per_sample/${sample}/qc", mode: 'copy'
  input: tuple val(sample), path(fq)
  output: tuple val(sample), path("${sample}.seqkit.stats.tsv")
  script: "seqkit stats -a -T ${fq} > ${sample}.seqkit.stats.tsv"
}

// TRIM und FILTER Prozesse sind hier gelöscht/ignoriert, 
// weil wir sie oben im Workflow übersprungen haben.

process FASTQ_TO_FASTA {
  tag "$sample"
  publishDir "output/per_sample/${sample}/clusters", mode: 'copy'
  input: tuple val(sample), path(fq)
  output: tuple val(sample), path("${sample}.fasta")
  script: "seqkit fq2fa ${fq} > ${sample}.fasta"
}

process CLUSTER_VSEARCH {
  tag "$sample"
  publishDir "output/per_sample/${sample}/clusters", mode: 'copy'
  input: tuple val(sample), path(fa)
  output: tuple val(sample), path("${sample}.centroids.fasta"), path("${sample}.clusters.uc")
  // WICHTIG: --minseqlength 10 erlaubt auch kurze Sequenzen
  script: "vsearch --cluster_fast ${fa} --id 0.85 --minseqlength 10 --strand both --threads ${task.cpus} --uc ${sample}.clusters.uc --centroids ${sample}.centroids.fasta"
}

process FILTER_CLUSTERS {
  tag "$sample"
  publishDir "output/per_sample/${sample}/clusters", mode: 'copy'
  input: tuple val(sample), path(centroids), path(uc)
  output: tuple val(sample), path("${sample}.cluster_counts.tsv"), path("${sample}.centroids.kept.fasta")
  script:
  """
  awk -F'\\t' 'BEGIN{OFS="\\t"} \$1=="S"{cl=\$2; id=\$9; c[cl]=id; n[id]=1} \$1=="H"{n[c[\$2]]++} END{for(i in n) print i,n[i]}' ${uc} | sort -k2,2nr > ${sample}.cluster_counts.tsv
  # WICHTIG: Filter auf >= 1 setzen (Hardcoded)
  awk '\$2>=1{print \$1}' ${sample}.cluster_counts.tsv > keep.txt
  if [ -s keep.txt ]; then seqkit grep -f keep.txt ${centroids} > ${sample}.centroids.kept.fasta; else touch ${sample}.centroids.kept.fasta; fi
  """
}

process MAKEBLASTDB {
  publishDir "output/refdb", mode: 'copy'
  container "ncbi/blast:2.16.0"
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
  publishDir "output/per_sample/${sample}/taxonomy", mode: 'copy'
  container "ncbi/blast:2.16.0"
  input: tuple val(sample), path(counts), path(centroids_fa), path(dbdir)
  output: tuple val(sample), path(counts), path("${sample}.blast.tsv")
  script:
  """
  if [ ! -s ${centroids_fa} ]; then touch ${sample}.blast.tsv; exit 0; fi
  blastn -query ${centroids_fa} -db ${dbdir}/fusarium_tef1 -max_target_seqs 5 -num_threads ${task.cpus} -outfmt "6 qseqid sseqid pident length qlen bitscore evalue qcovs stitle" > ${sample}.blast.tsv
  """
}

process JOIN_COUNTS_BLAST {
  tag "$sample"
  publishDir "output/per_sample/${sample}/taxonomy", mode: 'copy'
  input: tuple val(sample), path(counts_tsv), path(blast_tsv)
  output: tuple val(sample), path("${sample}.taxonomy.tsv")
  script:
  """
  python - << 'PY'
import csv
counts_file = "${counts_tsv}"
blast_file = "${blast_tsv}"
out_file = "${sample}.taxonomy.tsv"
hits = {}
try:
    with open(blast_file) as f:
        reader = csv.reader(f, delimiter='\\t')
        for row in reader:
            if not row: continue
            qid = row[0]
            if qid not in hits:
                stitle = row[8] if len(row) > 8 else 'NA'
                hits[qid] = {'sseqid': row[1], 'pident': row[2], 'qcovs': row[7], 'stitle': stitle}
except Exception: pass
with open(out_file, 'w') as out:
    out.write("cluster_id\\tread_count\\tbest_sseqid\\tpident\\tqcovs\\tbest_hit_title\\n")
    with open(counts_file) as f:
        reader = csv.reader(f, delimiter='\\t')
        for row in reader:
            if not row: continue
            cid = row[0]
            count = row[1]
            h = hits.get(cid, {'sseqid':'NA', 'pident':'NA', 'qcovs':'NA', 'stitle':'NA'})
            out.write(f"{cid}\\t{count}\\t{h['sseqid']}\\t{h['pident']}\\t{h['qcovs']}\\t{h['stitle']}\\n")
PY
  """
}

process AGGREGATE_RESULTS {
  publishDir "output/summary", mode: 'copy'
  input: path(tables)
  output: tuple path("all_samples.long.tsv"), path("abundance_matrix.tsv")
  script:
  """
  printf "sample\\tcluster_id\\tread_count\\tbest_sseqid\\tpident\\tqcovs\\tbest_hit_title\\n" > all_samples.long.tsv
  for f in ${tables}; do 
    fname=\$(basename \$f)
    s=\${fname%%.taxonomy.tsv}
    awk -v s=\$s 'NR>1{print s"\\t"\$0}' \$f >> all_samples.long.tsv
  done
  python - << 'PY'
import csv
samples = set()
taxa = set()
data = {}
with open('all_samples.long.tsv', 'r') as f:
    reader = csv.DictReader(f, delimiter='\\t')
    for row in reader:
        s = row['sample']
        t = row['best_hit_title']
        c = int(row['read_count'])
        samples.add(s)
        taxa.add(t)
        data[(s,t)] = data.get((s,t), 0) + c
sorted_samples = sorted(list(samples))
sorted_taxa = sorted(list(taxa))
with open('abundance_matrix.tsv', 'w') as out:
    out.write('sample\\t' + '\\t'.join(sorted_taxa) + '\\n')
    for s in sorted_samples:
        row_vals = [s]
        for t in sorted_taxa:
            row_vals.append(str(data.get((s,t), 0)))
        out.write('\\t'.join(row_vals) + '\\n')
PY
  """
}

process TSV_TO_CSV {
  publishDir "output/summary", mode: 'copy'
  input: path tsv
  output: path "${tsv.baseName}.csv"
  script: "sed 's/\\t/,/g' ${tsv} > ${tsv.baseName}.csv"
}

process REPORT_HTML {
  publishDir "output", mode: 'copy'
  input: tuple path(long_csv), path(mat_csv)
  output: path("wf-fusarium-tef1-report.html")
  script:
  """
  python - << 'PY'
from pathlib import Path
long_name = Path("${long_csv}").name
mat_name  = Path("${mat_csv}").name
html = f'''<!doctype html><html><head><meta charset="utf-8"/><title>Fusarium TEF1 Report</title>
<style>body {{ font-family: sans-serif; margin: 30px; }} .box {{ border: 1px solid #ddd; padding: 20px; background: #fafafa; }} a {{ color: #0366d6; }}</style>
</head><body><h1>Fusarium TEF1 Consensus Report</h1>
<div class="box"><h3>Ergebnisse</h3><ul>
<li><a href="summary/{long_name}">Detaillierte Tabelle (CSV)</a></li>
<li><a href="summary/{mat_name}">Abundanz Matrix (CSV)</a></li>
</ul></div></body></html>'''
with open("wf-fusarium-tef1-report.html", "w", encoding="utf-8") as f: f.write(html)
PY
  """
}
