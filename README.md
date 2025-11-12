# Appendix A: Bioinformatic Pipeline for Air Monitoring

## Overview

This repository documents the computational framework used in *Publication I: Air Monitoring by Nanopore Sequencing*. The pipeline is designed for ultra-low-biomass metagenomics, characterizing microbial diversity, function, and antimicrobial resistance (AMR) from active air samples.

## 1\. Basecalling & Demultiplexing

### 1.1. Controlled/Natural Environments (Guppy)

*Model:* High Accuracy (HAC) for R10.4.1 flow cells.

```bash
guppy_basecaller -i [input_raw_data_dir] -r -s [output_dir] --detect_barcodes -c dna_r10.4.1_e8.2_400bps_hac.cfg -x "cuda:0"
```

### 1.2. Urban Environment (Dorado)

*Model:* Super Accuracy (SUP) v5.0.0 with modified base retention (4mC/5mC).

```bash
# Basecalling to BAM
dorado basecaller dna_r10.4.1_e8.2_400bps_sup@v5.0.0 [input_pod5_dir] -r --kit-name SQK-RBK114-24 --no-trim --emit-fastq > [basecalled.fastq]

# Demultiplexing
dorado demux --output-dir [output_demux_dir] --kit-name SQK-RBK114-24 [basecalled.fastsq]
```

## 2\. Read Pre-processing

### 2.1. Adapter Trimming (Porechop)

*Purpose:* Removal of sequencing adapters and barcodes.

```bash
porechop -i [input_barcode.fastq] -o [trimmed.fastq] -t 10
```

### 2.2. Quality & Length Filtering (NanoFilt)

*Threshold:* Minimum length 100 bp.

```bash
cat [trimmed.fastq] | NanoFilt -q 9 -l 100 > [filtered.fastq]
```

### 2.3. Format Conversion (Seqtk)

*Purpose:* Generating FASTA files for downstream functional screening.

```bash
seqtk seq -A [filtered.fastq] > [filtered.fasta]
```

## 3\. Metagenomic Assembly & Polishing

### 3.1. De Novo Assembly (MetaFlye)

*Strategy:* Long-read metagenomic assembly using Nano-HQ mode.

```bash
flye --meta --nano-hq [filtered.fastq] --threads [threads] -o [assembly_dir]
```

### 3.2. Read Mapping for Polishing (Minimap2)

*Purpose:* Mapping reads to the draft assembly for consensus correction.

```bash
minimap2 -ax map-ont -t [threads] [assembly.fasta] [filtered.fastq] > [alignment.sam]
```

### 3.3. Alignment Processing (Samtools)

```bash
samtools view -b -@ [threads] [alignment.sam] | samtools sort -@ [threads] -o [sorted.bam]
```

### 3.4. Polishing (Racon)

*Purpose:* Consensus correction of the draft assembly.

```bash
racon -t [threads] [filtered.fastq] [alignment.sam] [assembly.fasta] > [polished_assembly.fasta]
```

## 4\. Taxonomic Classification

### 4.1. Read-Level Classification (Kraken2)

```bash
kraken2 --db [kraken_db] --use-names --report [report.txt] --output [output.txt] [filtered.fastq] --memory-mapping --threads 28
```

### 4.2. Contig-Level Classification (Kraken2)

```bash
kraken2 --db [kraken_db] --use-names --report [report_contig.txt] --output [output_contig.txt] [polished_assembly.fasta] --memory-mapping --threads 28
```

## 5\. Metagenomic Binning

### 5.1. Alignment for Coverage (Minimap2)

*Note:* The `-L` flag is used for long-read mapping context.

```bash
minimap2 -ax map-ont -L -t 20 [polished_assembly.fasta] [filtered.fastq] > [binning_alignment.sam]
```

### 5.2. BAM Sorting (Samtools)

```bash
samtools sort -T tmp-samtools -@ 20 -O BAM -o [binning_sorted.bam] [binning_alignment.sam]
samtools index -@ 20 -b [binning_sorted.bam]
```

### 5.3. Initial Binning (MetaWRAP)

*Strategy:* Ensemble binning using MetaBAT2, MaxBin2, and CONCOCT.

```bash
metawrap binning --metabat2 --maxbin2 --concoct -t 20 -m 64 --single-end --universal --run-checkm -l 10000 -a [polished_assembly.fasta] -o [output_dir] [filtered.fastq]
```

### 5.4. Bin Refinement (MetaWRAP)

*Thresholds:* \>50% Completeness (`-c 50`), \<10% Contamination (`-x 10`).

```bash
metawrap bin_refinement -o [refinement_dir] -t 20 -A [metabat2_bins] -B [maxbin2_bins] -C [concoct_bins] -c 50 -x 10
```

## 6\. Functional Annotation

### 6.1. Gene Prediction (Prodigal)

*Purpose:* Predicting open reading frames (ORFs) on polished contigs using metagenomic mode (`-p meta`).

```bash
prodigal -i [polished_assembly.fasta] -o [genes.gff] -a [proteins.faa] -p meta
```

### 6.2. Orthology Assignment (eggNOG-mapper)

*Database:* bactNOG; *Search Mode:* DIAMOND.

```bash
emapper.py -m diamond --data_dir [data_dir] -d bactNOG -i [proteins.faa] --output [output_prefix] --override --target_orthologs all --query-cover 20 --subject-cover 20 --tax_scope auto --cpu [threads]
```

### 6.3. Protein Alignment (DIAMOND BLASTx)

*Purpose:* Aligning reads or contigs against large protein databases (e.g., NR or VFDB).

**For Reads (FASTQ):**

```bash
diamond blastx -d [database.dmnd] -q [filtered.fastq] -o [output.dmnd_out] --max-target-seqs 1 --threads 20 -f 6 qseqid sseqid pident length qstart qend sstart send evalue bitscore staxids salltitles sscinames
```

**For Contigs (FASTA):**

```bash
diamond blastx -d [database.dmnd] -q [polished_assembly.fasta] -o [output.dmnd_out] --max-target-seqs 1 --threads 20 -f 6 qseqid sseqid pident length qstart qend sstart send evalue bitscore staxids salltitles sscinames
```

### 6.4. Antimicrobial Resistance Detection

**AMRFinderPlus:**

```bash
amrfinder --threads 10 -n [filtered.fasta] -d [database_path] > [output_amrfinder.txt]
```

**ABRicate:**
*Note:* Run for specific databases (e.g., `ncbi`, `card`, `resfinder`).

```bash
abricate --db [database_name] [filtered.fasta] --threads 10 > [output_abricate.txt]
```

**ABRicate Summary:**

```bash
abricate --summary [list_of_output_files] > [summary.tab]
```
