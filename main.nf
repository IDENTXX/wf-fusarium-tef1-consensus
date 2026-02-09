nextflow.enable.dsl=2

workflow {
    // 1. Check ob Parameter da sind (für EPI2ME Refresh)
    if( !params.reads || !params.db_fasta ) {
        log.info "Warte auf Parameter..."
        return
    }

    // 2. Input Kanäle erstellen
    Channel
        .fromPath("${params.reads}/*/*.fastq.gz", checkIfExists: true)
        .map { f -> tuple(f.parent.name, f) }
        .groupTuple()
        .map { sample, files -> tuple(sample, files.sort()) }
        .set { ch_samples }

    // 3. Der Workflow-Ablauf
    merged   = MERGE_FASTQ(ch_samples)
    qc_raw   = QC_SEQKIT(merged)

    // Optionales Primer Trimming
    trimmed  = (params.primer_fwd?.trim() && params.primer_rev?.trim()) 
             ? TRIM_CUTADAPT(merged) 
             : merged
    
    // Qualitäts-Filter
    filtered = FILTER_NANOFILT(trimmed)
    
    // Konvertierung für VSEARCH
    fasta    = FASTQ_TO_FASTA(filtered)

    // WICHTIG: Consensus Generierung (Poliert die Sequenzen)
    consensus = CLUSTER_CONSENSUS(fasta)

    // Datenbank vorbereiten
    ref_fa   = Channel.value( file(params.db_fasta) )
    dbdir_ch = MAKEBLASTDB(ref_fa)

    // BLAST Suche
    blast_in = consensus.combine(dbdir_ch)
    blasted  = BLASTN(blast_in)
    tax      = JOIN_COUNTS_BLAST(blasted)

    // Zusammenfassung & Report
    tax_tables_list = tax.map { _, taxfile -> taxfile }.collect()
    summary_ch      = AGGREGATE_RESULTS(tax_tables_list)
    
    REPORT_HTML(summary_ch)
}

// ------------------------------------------------------
// PROZESS DEFINITIONEN
// ------------------------------------------------------

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
    -g ${params.primer_fwd} -a ${params.primer_rev} ${discard} \
    -o ${sample}.trimmed.fastq.gz ${fq} > ${sample}.cutadapt.log
    """
}

process FILTER_NANOFILT {
    tag "$sample"
    publishDir "${params.out_dir}/per_sample/${sample}/reads", mode: 'copy'

    input: tuple val(sample), path(fq)
    output: tuple val(sample), path("${sample}.filtered.fastq.gz")

    script:
    """
    python - << 'PY'
    import gzip
    infile, outfile = "${fq}", "${sample}.filtered.fastq.gz"
    min_q, min_len, max_len = int(${params.min_q}), int(${params.min_len}), int(${params.max_len})
    
    def mean(q): return sum(ord(c)-33 for c in q)/len(q) if q else 0

    with gzip.open(infile,"rt") as fi, gzip.open(outfile,"wt") as fo:
        while True:
            h = fi.readline()
            if not h: break
            seq, _, qual = fi.readline().strip(), fi.readline(), fi.readline().strip()
            if min_len <= len(seq) <= max_len and mean(qual) >= min_q:
                fo.write(h + seq + "\\n+\\n" + qual + "\\n")
    PY
    """
}

process FASTQ_TO_FASTA {
    tag "$sample"
    input: tuple val(sample), path(fq)
    output: tuple val(sample), path("${sample}.fasta")
    script: "seqkit fq2fa ${fq} > ${sample}.fasta"
}

process CLUSTER_CONSENSUS {
    tag "$sample"
    publishDir "${params.out_dir}/per_sample/${sample}/clusters", mode: 'copy'

    input: tuple val(sample), path(fa)
    output: tuple val(sample), path("${sample}.cons.fasta"), path("${sample}.cluster_counts.tsv")

    script:
    """
    # 1. Clustering & Consensus (MSA)
    vsearch --cluster_fast ${fa} --id ${params.cluster_id} --strand both \
            --minseqlength ${params.min_len} --consout temp_consensus.fasta \
            --sizeout --threads ${task.cpus}

    # 2. Filtern & Header bereinigen
    if [ -s temp_consensus.fasta ]; then
        seqkit seq -n temp_consensus.fasta | \
        awk -F'[=;]' -v min=${params.min_cluster_reads} '\$4 >= min {print \$0}' > keep_ids.txt
        
        if [ -s keep_ids.txt ]; then
            seqkit grep -f keep_ids.txt temp_consensus.fasta > ${sample}.cons.fasta
            
            # Count Table erstellen
            grep ">" ${sample}.cons.fasta | sed 's/>//' | \
            awk -F'[=;]' '{print \$2"\\t"\$4}' > ${sample}.cluster_counts.tsv
        else
            touch ${sample}.cons.fasta; touch ${sample}.cluster_counts.tsv
        fi
    else
        touch ${sample}.cons.fasta; touch ${sample}.cluster_counts.tsv
    fi
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
    input: tuple val(sample), path(cons_fa), path(counts), path(dbdir)
    output: tuple val(sample), path(counts), path("${sample}.blast.tsv")

    script:
    """
    if [ ! -s ${cons_fa} ]; then
        touch ${sample}.blast.tsv
    else
        blastn -query ${cons_fa} -db ${dbdir}/fusarium_tef1 \
            -max_target_seqs ${params.blast_topn} -num_threads ${task.cpus} \
            -outfmt "6 qseqid sseqid pident length qlen bitscore evalue qcovs stitle" \
            > ${sample}.blast.tsv
    fi
    """
}

process JOIN_COUNTS_BLAST {
    tag "$sample"
    publishDir "${params.out_dir}/per_sample/${sample}/taxonomy", mode: 'copy'

    input: tuple val(sample), path(counts_tsv), path(blast_tsv)
    output: tuple val(sample), path("${sample}.taxonomy.tsv")

    script:
    """
    python - << 'PY'
    import pandas as pd
    import sys

    try: counts = pd.read_csv("${counts_tsv}", sep="\\t", names=["id", "count"])
    except: counts = pd.DataFrame(columns=["id","count"])

    try: 
        blast = pd.read_csv("${blast_tsv}", sep="\\t", names=["qseqid", "sseqid", "pident", "length", "qlen", "bitscore", "evalue", "qcovs", "stitle"])
        blast = blast.sort_values("bitscore", ascending=False).drop_duplicates("qseqid")
    except: blast = pd.DataFrame(columns=["qseqid", "sseqid", "pident", "qcovs", "stitle"])

    merged = pd.merge(counts, blast, left_on="id", right_on="qseqid", how="left")
    merged["best_hit_title"] = merged["stitle"].fillna("No_Hit")
    
    out_cols = ["id", "count", "sseqid", "pident", "qcovs", "best_hit_title"]
    merged[out_cols].to_csv("${sample}.taxonomy.tsv", sep="\\t", index=False)
    PY
    """
}

process AGGREGATE_RESULTS {
    publishDir "${params.out_dir}/summary", mode: 'copy'

    input: path(tables)
    output: tuple path("all_samples.long.csv"), path("abundance_matrix.csv")

    script:
    """
    python - << 'PY'
    import pandas as pd
    import glob

    files = glob.glob("*.taxonomy.tsv")
    dfs = []
    for f in files:
        tmp = pd.read_csv(f, sep="\\t")
        tmp["sample"] = f.replace(".taxonomy.tsv", "")
        dfs.append(tmp)

    if not dfs:
        pd.DataFrame().to_csv("all_samples.long.csv", index=False)
        pd.DataFrame().to_csv("abundance_matrix.csv", index=False)
    else:
        full = pd.concat(dfs, ignore_index=True)
        full.to_csv("all_samples.long.csv", index=False)
        
        matrix = full.pivot_table(index="sample", columns="best_hit_title", values="count", aggfunc="sum", fill_value=0)
        matrix.to_csv("abundance_matrix.csv")
    PY
    """
}

process REPORT_HTML {
    publishDir "${params.out_dir}", mode: 'copy'

    input: tuple path(long_csv), path(mat_csv)
    output: path("wf-fusarium-tef1-report.html")

    script:
    """
    python - << 'PY'
    import pandas as pd
    from datetime import datetime

    try: mat = pd.read_csv("${mat_csv}")
    except: mat = pd.DataFrame()

    html = f'''<!doctype html>
    <html lang="en">
    <head><meta charset="utf-8"><title>Fusarium Report</title>
    <style>body{{font-family:sans-serif;padding:20px}} table{{border-collapse:collapse;width:100%}} th,td{{border:1px solid #ddd;padding:8px}}</style>
    </head><body>
    <h1>Fusarium <i>tef1</i> Consensus Report</h1>
    <p>Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}</p>
    <ul>
        <li><a href="summary/${long_csv}">Download Full Table</a></li>
        <li><a href="summary/${mat_csv}">Download Matrix</a></li>
    </ul>
    <h3>Abundance Preview</h3>
    {mat.to_html(classes="table", index=False, max_rows=10) if not mat.empty else "No data."}
    </body></html>'''
    
    with open("wf-fusarium-tef1-report.html", "w") as f: f.write(html)
    PY
    """
}