nextflow.enable.dsl=2

/*
  Haupt-Workflow: FUSARIUM TEF1 (Final Version)
  - Auto-Datums-Stempel (Archivierung)
  - EPI2ME Sync (Anzeige im Browser)
  - Scharfstell-Filter (Identity/Coverage)
*/

// --- 1. DATUMS-GENERATOR ---
def timestamp = new java.util.Date().format( 'yyyy-MM-dd_HH-mm' )

// Logik: Zeitstempel-Ordner, außer der User tippt einen eigenen Namen
params.final_outdir = params.out_dir == 'output' ? "Fusarium_Run_${timestamp}" : params.out_dir

workflow {
  // Check Inputs
  if( !params.reads || !params.db_fasta ) {
    log.info "Parameter fehlen. Bitte Reads und DB angeben."
    return
  }
  
  log.info "------------------------------------------------"
  log.info " FUSARIUM TEF1 PIPELINE"
  log.info " Ergebnisse landen in: ${params.final_outdir}"
  log.info "------------------------------------------------"

  // 1. Input Kanäle
  Channel
    .fromPath("${params.reads}/*/*.fastq.gz", checkIfExists: true)
    .map { f -> tuple(f.parent.name, f) }
    .groupTuple()
    .map { sample, files -> tuple(sample, files.sort()) }
    .set { ch_samples }

  // 2. Workflow Logik
  merged   = MERGE_FASTQ(ch_samples)
  qc_raw   = QC_SEQKIT(merged) 

  // Bypass: Direkt Fasta
  fasta     = FASTQ_TO_FASTA(merged)
  
  // Clustering
  clustered = CLUSTER_VSEARCH(fasta)
  kept      = FILTER_CLUSTERS(clustered)
  
  // DB
  ref_fa   = Channel.value( file(params.db_fasta) )
  dbdir_ch = MAKEBLASTDB(ref_fa)
  
  // BLAST
  kept_with_db = kept.combine(dbdir_ch)
  blasted      = BLASTN(kept_with_db)
  
  // FILTERUNG (Scharfstellen)
  tax          = JOIN_COUNTS_BLAST(blasted)
  
  // Aggregation
  tax_tables_list = tax.map { sample, taxfile -> taxfile }.collect() 
  summary = AGGREGATE_RESULTS(tax_tables_list)
  
  // Bericht
  REPORT_HTML(summary)
}

// ------------------------------------------------------
// PROZESSE
// ------------------------------------------------------

process MERGE_FASTQ {
  tag "$sample"
  publishDir "${params.final_outdir}/per_sample/${sample}/reads", mode: 'copy'
  
  input: tuple val(sample), path(reads)
  output: tuple val(sample), path("${sample}.fastq.gz")
  script: "cat ${reads.join(' ')} > ${sample}.fastq.gz"
}

process QC_SEQKIT {
  tag "$sample"
  publishDir "${params.final_outdir}/per_sample/${sample}/qc", mode: 'copy'
  
  input: tuple val(sample), path(fq)
  output: tuple val(sample), path("${sample}.seqkit.stats.tsv")
  script: "seqkit stats -a -T ${fq} > ${sample}.seqkit.stats.tsv"
}

process FASTQ_TO_FASTA {
  tag "$sample"
  publishDir "${params.final_outdir}/per_sample/${sample}/clusters", mode: 'copy'
  
  input: tuple val(sample), path(fq)
  output: tuple val(sample), path("${sample}.fasta")
  script: "seqkit fq2fa ${fq} > ${sample}.fasta"
}

process CLUSTER_VSEARCH {
  tag "$sample"
  publishDir "${params.final_outdir}/per_sample/${sample}/clusters", mode: 'copy'
  
  input: tuple val(sample), path(fa)
  output: tuple val(sample), path("${sample}.centroids.fasta"), path("${sample}.clusters.uc")
  script: "vsearch --cluster_fast ${fa} --id 0.85 --minseqlength 10 --strand both --threads ${task.cpus} --uc ${sample}.clusters.uc --centroids ${sample}.centroids.fasta"
}

process FILTER_CLUSTERS {
  tag "$sample"
  publishDir "${params.final_outdir}/per_sample/${sample}/clusters", mode: 'copy'
  
  input: tuple val(sample), path(centroids), path(uc)
  output: tuple val(sample), path("${sample}.cluster_counts.tsv"), path("${sample}.centroids.kept.fasta")
  script:
  """
  awk -F'\\t' 'BEGIN{OFS="\\t"} \$1=="S"{cl=\$2; id=\$9; c[cl]=id; n[id]=1} \$1=="H"{n[c[\$2]]++} END{for(i in n) print i,n[i]}' ${uc} | sort -k2,2nr > ${sample}.cluster_counts.tsv
  awk '\$2>=1{print \$1}' ${sample}.cluster_counts.tsv > keep.txt
  if [ -s keep.txt ]; then seqkit grep -f keep.txt ${centroids} > ${sample}.centroids.kept.fasta; else touch ${sample}.centroids.kept.fasta; fi
  """
}

process MAKEBLASTDB {
  publishDir "${params.final_outdir}/refdb", mode: 'copy'
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
  publishDir "${params.final_outdir}/per_sample/${sample}/taxonomy", mode: 'copy'
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
  publishDir "${params.final_outdir}/per_sample/${sample}/taxonomy", mode: 'copy'
  
  input: tuple val(sample), path(counts_tsv), path(blast_tsv)
  output: tuple val(sample), path("${sample}.taxonomy.tsv")
  
  script:
  """
  python - << 'PY'
import csv
counts_file = "${counts_tsv}"
blast_file = "${blast_tsv}"
out_file = "${sample}.taxonomy.tsv"

# SCHWELLENWERTE aus UI
min_id = float("${params.min_identity}")
min_cov = float("${params.min_coverage}")

hits = {}
try:
    with open(blast_file) as f:
        reader = csv.reader(f, delimiter='\\t')
        for row in reader:
            if not row: continue
            qid = row[0]
            if qid in hits: continue
            
            pident = float(row[2])
            qcovs = float(row[7])
            
            # FILTER
            if pident < min_id or qcovs < min_cov:
                continue 
            
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
            h = hits.get(cid, {'sseqid':'NA', 'pident':'0', 'qcovs':'0', 'stitle':'Unclassified'})
            out.write(f"{cid}\\t{count}\\t{h['sseqid']}\\t{h['pident']}\\t{h['qcovs']}\\t{h['stitle']}\\n")
PY
  """
}

process AGGREGATE_RESULTS {
  publishDir "${params.final_outdir}/summary", mode: 'copy'
  publishDir "output/summary", mode: 'copy', overwrite: true

  input: path(tables)
  output: tuple path("wf-metagenomics-counts-species.csv"), path("abundance_matrix.csv")

  script:
  """
  printf "sample\\tcluster_id\\tread_count\\tbest_sseqid\\tpident\\tqcovs\\tbest_hit_title\\n" > raw_combined.tsv
  for f in ${tables}; do
    fname=\$(basename \$f)
    s=\${fname%%.taxonomy.tsv}
    awk -v s=\$s 'NR>1{print s"\\t"\$0}' \$f >> raw_combined.tsv
  done

  python - << 'PY'
import csv
import re
input_file = 'raw_combined.tsv'
output_target = 'wf-metagenomics-counts-species.csv'
output_matrix = 'abundance_matrix.csv'

data = {} 
all_samples = set()
all_species = set()

def clean_name(raw_title):
    if "Unclassified" in raw_title: return "Unclassified"
    
    match_fsp = re.search(r'(Fusarium\\s+[\\w\\.]+\\s+f\\.\\s*sp\\.\\s+[\\w\\.]+)', raw_title, re.IGNORECASE)
    if match_fsp:
        val = match_fsp.group(1)
        return val[0].upper() + val[1:]

    match_normal = re.search(r'(Fusarium\\s+[\\w\\.]+)', raw_title, re.IGNORECASE)
    if match_normal:
        return match_normal.group(1).capitalize()
    
    parts = raw_title.split()
    if len(parts) >= 2: return f"{parts[0]} {parts[1]}"
    return raw_title

with open(input_file, 'r') as f:
    reader = csv.DictReader(f, delimiter='\\t')
    for row in reader:
        sample = row['sample']
        raw_name = row['best_hit_title']
        try: count = int(row['read_count'])
        except: count = 0
        species = clean_name(raw_name)
        all_samples.add(sample)
        all_species.add(species)
        if species not in data: data[species] = {}
        data[species][sample] = data[species].get(sample, 0) + count

sorted_samples = sorted(list(all_samples))
sorted_species = sorted(list(all_species))
if "Unclassified" in sorted_species:
    sorted_species.remove("Unclassified")
    sorted_species.append("Unclassified")

with open(output_target, 'w') as out:
    header = ['species'] + sorted_samples + ['total', 'superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'tax']
    out.write(','.join(header) + '\\n')
    for sp in sorted_species:
        total = sum(data[sp].values())
        row = [sp]
        for sa in sorted_samples: row.append(str(data[sp].get(sa, 0.0)))
        row.append(str(float(total)))
        
        genus = sp.split()[0] if sp else "Unknown"
        sk, k, p, c, o, fa = "Eukaryota", "Fungi", "Ascomycota", "Sordariomycetes", "Hypocreales", "Nectriaceae"
        if genus != "Fusarium":
             p = c = o = fa = "NA"
             if sp == "Unclassified": k = "NA"

        tax_s = f"{sk};{k};{p};{c};{o};{fa};{genus};{sp}"
        row.extend([sk, k, p, c, o, fa, genus, tax_s])
        out.write(','.join(row) + '\\n')

with open(output_matrix, 'w') as out:
    header = ['species'] + sorted_samples
    out.write(','.join(header) + '\\n')
    for sp in sorted_species:
        row = [sp]
        for sa in sorted_samples: row.append(str(data[sp].get(sa, 0)))
        out.write(','.join(row) + '\\n')
PY
  """
}

process REPORT_HTML {
  publishDir "${params.final_outdir}", mode: 'copy'
  publishDir "output", mode: 'copy', overwrite: true

  input: tuple path(metagenomics_csv), path(matrix_csv)
  output: path("wf-fusarium-tef1-report.html")

  script:
  """
  python - << 'PY'
from pathlib import Path
csv_name = Path("${metagenomics_csv}").name
mat_name = Path("${matrix_csv}").name
html = f'''<!doctype html><html><head><meta charset="utf-8"/><title>Fusarium Report</title>
<style>body{{font-family:sans-serif;margin:30px;}}.box{{border:1px solid #ddd;padding:20px;background:#fafafa;}}a{{color:#0366d6;text-decoration:none;font-weight:bold;}}</style>
</head><body><h1>Fusarium TEF1 Report</h1><div class="box"><h3>Ergebnisse</h3><ul>
<li><a href="summary/{csv_name}" target="_blank">Metagenomics Counts (Final)</a></li>
<li><a href="summary/{mat_name}" target="_blank">Abundance Matrix</a></li>
</ul></div></body></html>'''
with open("wf-fusarium-tef1-report.html", "w", encoding="utf-8") as f: f.write(html)
PY
  """
}