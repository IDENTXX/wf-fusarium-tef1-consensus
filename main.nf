process AGGREGATE_RESULTS {
  publishDir "${params.out_dir}/summary", mode: 'copy'
  // Wir nutzen hier pandas, das ist im python-slim container enthalten
  container "identxx/fusarium-toolbox:v1.0" 

  input:
  path(tables)

  output:
  tuple path("all_samples.long.tsv"), path("abundance_matrix.tsv")

  script:
  """
  # 1. Header sicher mit printf schreiben (damit \t als Tab erkannt wird)
  printf "sample\\tcluster_id\\tcount\\tbest_hit\\tpident\\tstitle\\n" > all_samples.long.tsv

  # 2. Alle Einzeldateien anfügen (Header überspringen)
  for f in ${tables}; do 
    s=\${f%%.taxonomy.tsv}
    awk -v s=\$s 'NR>1{print s"\\t"\$0}' \$f >> all_samples.long.tsv
  done
  
  # 3. Matrix erstellen mit Python (Fehlerabfangung eingebaut)
  python -c "import pandas as pd; 
try:
    df=pd.read_csv('all_samples.long.tsv', sep='\\t'); 
    if not df.empty and 'count' in df.columns:
        m=df.pivot_table(index='sample', columns='stitle', values='count', fill_value=0); 
        m.to_csv('abundance_matrix.tsv', sep='\\t')
    else:
        print('Warnung: Keine Daten vorhanden, erstelle leere Matrix.')
        pd.DataFrame(columns=['sample']).to_csv('abundance_matrix.tsv', sep='\\t')
except Exception as e:
    print(f'Fehler in Python: {e}')
    # Fallback leere Datei, damit Nextflow nicht abbricht
    pd.DataFrame(columns=['sample']).to_csv('abundance_matrix.tsv', sep='\\t')
"
  """
}
