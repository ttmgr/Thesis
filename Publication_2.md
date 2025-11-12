# Appendix B: Bioinformatic Pipeline for Wetland Ecosystems

## Overview

This repository documents the computational framework used in *Publication II: Real-time Pathogen, Resistance, and Host Range Characterization*. The pipeline integrates multi-omic data from passive water samplers, utilizing both standard long-read assemblers and specialized De Bruijn graph methods for fragmented environmental DNA.

## 1\. Basecalling & Demultiplexing

### 1.1. Basecalling (Dorado)

*Model:* Super Accuracy (SUP) v5.0.0.

```bash
dorado basecaller dna_r10.4.1_e8.2_400bps_sup@v5.0.0 [input_pod5_dir] -r --kit-name SQK-RBK114-24 --no-trim --emit-fastq > [basecalled.fastq]
```

### 1.2. Demultiplexing (Dorado)

```bash
dorado demux --output-dir [output_dir] --kit-name SQK-RBK114-24 [basecalled.fastq] --emit-fastq
```

## 2\. Read Pre-processing

### 2.1. Adapter Trimming (Porechop)

*Purpose:* Removal of sequencing adapters and barcodes.

```bash
porechop -i [input_barcode.fastq] -o [trimmed.fastq] -t 20
```

### 2.2. Quality & Length Filtering (NanoFilt)

*Threshold:* Minimum length 100 bp.

```bash
cat [trimmed.fastq] | NanoFilt -l 100 > [filtered.fastq]
```

## 3\. Metagenomic Assembly Strategy I: OLC (Overlap-Layout-Consensus)

### 3.1. De Novo Assembly (MetaFlye)

*Strategy:* Long-read assembly using Nano-HQ mode optimized for metagenomes.

```bash
flye --meta --nano-hq [filtered.fastq] --threads [threads] -o [flye_output_dir]
```

### 3.2. Alignment for Polishing (Minimap2)

*Purpose:* Mapping reads to the draft assembly for consensus correction.

```bash
minimap2 -ax map-ont -t [threads] [assembly.fasta] [filtered.fastq] > [alignment.sam]
```

### 3.3. BAM Sorting (Samtools)

```bash
samtools view -b -@ [threads] [alignment.sam] | samtools sort -@ [threads] -o [sorted.bam]
```

### 3.4. Polishing (Racon)

*Purpose:* Consensus correction of the MetaFlye draft assembly.

```bash
racon -t [threads] [filtered.fastq] [alignment.sam] [assembly.fasta] > [polished_assembly.fasta]
```

## 4\. Metagenomic Assembly Strategy II: DBG (De Bruijn Graph)

### 4.1. De Novo Assembly (metaMDBG)

*Strategy:* Specialized k-mer based assembly for long reads.

```bash
metaMDBG asm --out-dir [mdbg_output_dir] --in-ont [filtered.fastq] --threads 12
```

## 5\. Taxonomic Classification (Kraken2)

### 5.1. Read-Level Classification

*Database:* Custom NCBI database.

```bash
kraken2 --db [kraken_db] --use-names --report [report_reads.txt] --output [output_reads.txt] [filtered.fastq] --threads 28
```

### 5.2. Assembly-Level Classification (Racon/Flye)

*Purpose:* Taxonomically classifying polished contigs from Strategy I.

```bash
kraken2 --db [kraken_db] --use-names --report [report_flye.txt] --output [output_flye.txt] [polished_assembly.fasta] --threads 28
```

### 5.3. Assembly-Level Classification (metaMDBG)

*Purpose:* Taxonomically classifying contigs from Strategy II (handling compressed input).

```bash
# Note: Process assumes unzipping if output is .gz
kraken2 --db [kraken_db] --use-names --report [report_mdbg.txt] --output [output_mdbg.txt] [mdbg_contigs.fasta] --threads 28
```

## 6\. Functional Annotation (AMR)

### 6.1. Format Conversion (Seqtk)

*Purpose:* Converting FASTQ reads to FASTA format for AMRFinderPlus.

```bash
seqtk seq -A [filtered.fastq] > [filtered.fasta]
```

### 6.2. Resistance Detection (AMRFinderPlus)

*Mode:* Plus (includes stress response and virulence genes).

```bash
amrfinder --plus --threads 10 -n [filtered.fasta] -d [database_path] > [output_amrfinder.txt]
```
