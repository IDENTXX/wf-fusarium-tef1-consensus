# Basis-Image: Schlankes Python 3.11
FROM python:3.11-slim

# Metadaten
LABEL maintainer="IDENTXX"
LABEL description="Fusarium TEF1 Toolbox: VSEARCH, SeqKit, Cutadapt, Pandas"

# 1. System-Updates und Tools installieren
# procps ist wichtig für Nextflow-Monitoring, build-essential für manche Python-Pakete
RUN apt-get update && apt-get install -y \
    wget \
    curl \
    procps \
    build-essential \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# 2. VSEARCH installieren (Statische Binary v2.28.1)
# Wir laden es herunter, entpacken es und verschieben es in den Pfad
RUN wget https://github.com/torognes/vsearch/releases/download/v2.28.1/vsearch-2.28.1-linux-x86_64.tar.gz && \
    tar xzf vsearch-2.28.1-linux-x86_64.tar.gz && \
    cp vsearch-2.28.1-linux-x86_64/bin/vsearch /usr/local/bin/ && \
    rm -rf vsearch*

# 3. SeqKit installieren (v2.6.1)
# Das "Schweizer Taschenmesser" für Fasta/Fastq
RUN wget https://github.com/shenwei356/seqkit/releases/download/v2.6.1/seqkit_linux_amd64.tar.gz && \
    tar -zxvf seqkit_linux_amd64.tar.gz && \
    mv seqkit /usr/local/bin/ && \
    rm seqkit_linux_amd64.tar.gz

# 4. Python-Bibliotheken installieren
# Cutadapt für Primer-Trimming, Pandas für die Tabellen-Erstellung
RUN pip install --no-cache-dir cutadapt==4.9 pandas

# Arbeitsverzeichnis setzen (Standard für Nextflow)
WORKDIR /data