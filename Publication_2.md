# Appendix B: Bioinformatic Pipeline for Wetland Ecosystems

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

*Strategy:* Long-read assembly using Nano-HQ mode.

```bash
flye --meta --nano-hq [filtered.fastq] --threads [threads] -o [flye_output_dir]
```

### 3.2. Consensus Correction (Racon)

*Purpose:* Initial polishing of the draft assembly using mapped reads.

```bash
minimap2 -ax map-ont -t [threads] [assembly.fasta] [filtered.fastq] > [alignment.sam]
samtools view -b -@ [threads] [alignment.sam] | samtools sort -@ [threads] -o [sorted.bam]
racon -t [threads] [filtered.fastq] [alignment.sam] [assembly.fasta] > [racon_polished.fasta]
```

### 3.3. Final Polishing (Medaka)

*Purpose:* Neural network-based polishing applied to the Racon-corrected contigs.

```bash
mini_align -t [threads] -i [filtered.fastq] -r [racon_polished.fasta] -m -p [output_prefix] -f
medaka inference [output_prefix].bam [output.hdf] --threads [threads] --model r1041_e82_400bps_sup_v5.0.0 --batch 200
medaka sequence [output.hdf] [racon_polished.fasta] [final_assembly.fasta]
```

## 4\. Metagenomic Assembly Strategy II: DBG (De Bruijn Graph)

### 4.1. De Novo Assembly (metaMDBG)

*Strategy:* Specialized k-mer based assembly for fragmented long reads.

```bash
metaMDBG asm --out-dir [mdbg_output_dir] --in-ont [filtered.fastq] --threads 12
```

### 4.2. Polishing (Medaka)

*Purpose:* Polishing the raw DBG assembly. Note: Input contigs are decompressed prior to alignment.

```bash
gunzip -c [mdbg_contigs.fasta.gz] > [temp_draft.fasta]
mini_align -t [threads] -i [filtered.fastq] -r [temp_draft.fasta] -m -p [output_prefix] -f
medaka inference [output_prefix].bam [output.hdf] --threads [threads] --model r1041_e82_400bps_sup_v5.0.0 --batch 200
medaka sequence [output.hdf] [temp_draft.fasta] [mdbg_polished.fasta]
```

## 5\. Taxonomic Classification (Kraken2)

### 5.1. Read-Level Classification

*Database:* Custom NCBI database.

```bash
kraken2 --db [kraken_db] --use-names --report [report_reads.txt] --output [output_reads.txt] [filtered.fastq] --threads 28
```

### 5.2. Assembly-Level Classification

*Applied to both Strategy I (Racon/Medaka) and Strategy II (MDBG) assemblies.*

```bash
kraken2 --db [kraken_db] --use-names --report [report_assembly.txt] --output [output_assembly.txt] [polished_assembly.fasta] --threads 28
```

## 6\. Pathogen Detection & MEGAN Post-Processing

### 6.1. FASTQ to FASTA Conversion (Seqtk)

*Purpose:* Converting sequencing reads for sorting and mapping.

```bash
seqtk seq -a [input.fastq] > [raw.fasta]
```

### 6.2. Read Sorting (SeqKit)

*Purpose:* Sorting FASTA files by read ID to optimize downstream alignment.

```bash
seqkit sort -n -w 0 --quiet [raw.fasta] -o [sorted.fasta]
```

### 6.3. Read-Level Alignment (Minimap2)

*Purpose:* Mapping reads to the NCBI-NT MMI index.
*Key Flags:* `-N 100` (retain many secondary alignments for LCA), `--lj-min-ratio 0.5` (stringency).

```bash
minimap2 -ax map-ont -k 19 -w 10 -I 10G -g 5000 -r 2000 -N 100 --lj-min-ratio 0.5 -A 2 -B 5 -O 5,56 -E 4,1 -z 400,50 --sam-hit-only -t [threads] --split-prefix [temp_idx] [minimap2_db_mmi] [sorted.fasta] > [aligned.sam]
```

### 6.4. SAM to RMA Conversion (MEGAN6)

*Purpose:* Converting SAM alignments to RMA format for MEGAN analysis.
*Filtered (Stringent):*

```bash
sam2rma -i [aligned.sam] -r [sorted.fasta] -o [filtered.rma] -lg -alg longReads -t [threads] -mdb [megan_db_nucl] -ram readCount --minSupportPercent 0.01
```

*Unfiltered:*

```bash
sam2rma -i [aligned.sam] -r [sorted.fasta] -o [unfiltered.rma] -lg -alg longReads -t [threads] -mdb [megan_db_nucl] -ram readCount --minSupportPercent 0
```

### 6.5. Taxonomic Information Extraction (rma2info)

*Purpose:* Extracting NCBI Taxonomy assignments from RMA files.
*Read-to-Class (r2c):*

```bash
rma2info -i [filtered.rma] -o [taxonomy.r2c.txt] -r2c Taxonomy -n
```

*Class-Count (c2c):*

```bash
rma2info -i [filtered.rma] -c2c Taxonomy -n -r -o [taxonomy.c2c.txt]
```

### 6.6. Report Conversion

*Purpose:* Converting MEGAN class counts (c2c) into standard Kraken and MetaPhlAn report formats.

```bash
python3 Convert_MEGAN_RMA_NCBI_c2c-snake.py --input [taxonomy.c2c.txt] --outname1 [names.temp.txt] --outname2 [codes.temp.txt] --mpa [report.mpa.txt] --kreport [report.kreport] --readsfile [read_counts.txt]
```

### 6.7. Assembly-Level Pathogen Detection (metaMDBG)

*Purpose:* Stringent pathogen identification applied directly to De Bruijn Graph assemblies to resolve gene-host context.

**Step 1: Decompression**

```bash
gunzip -c [mdbg_contigs.fasta.gz] > [mdbg_contigs.fasta]
```

**Step 2: Contig Alignment (Minimap2)**

```bash
minimap2 -ax map-ont -k 19 -w 10 -I 10G -g 5000 -r 2000 -N 100 --lj-min-ratio 0.5 -A 2 -B 5 -O 5,56 -E 4,1 -z 400,50 --sam-hit-only -t [threads] [minimap2_db_mmi] [mdbg_contigs.fasta] > [mdbg_aligned.sam]
```

**Step 3: RMA Conversion & Classification (MEGAN6)**

```bash
sam2rma -i [mdbg_aligned.sam] -r [mdbg_contigs.fasta] -o [mdbg.rma] -lg -alg longReads -t [threads] -mdb [megan_db_nucl] -ram readCount --minSupportPercent 0.01
```

## 7\. Genome Annotation & Gene Prediction

### 7.1. Gene Prediction (Prodigal)

*Purpose:* Predicting open reading frames (ORFs) on metaMDBG assemblies.

```bash
# Note: Uses the decompressed fasta from Step 6.7
prodigal -i [mdbg_contigs.fasta] -o [output.gff] -a [output.faa] -p meta -f gff
```

### 7.2. Functional Annotation (Prokka)

*Purpose:* Rapid annotation of prokaryotic genomes. Applied to both Racon-polished and MDBG assemblies.

```bash
prokka --force --quiet --outdir [output_dir] --prefix [sample_id] --cpus [threads] [input_assembly.fasta]
```

## 8\. Functional Annotation (Antimicrobial Resistance)

### 8.1. Format Conversion

```bash
seqtk seq -A [filtered.fastq] > [filtered.fasta]
```

### 8.2. Resistance Detection (AMRFinderPlus)

*Mode:* Plus (includes stress response and virulence genes).

```bash
amrfinder --plus --threads 10 -n [filtered.fasta] -d [database_path] > [output_amrfinder.txt]
```

## 9\. Targeted Analysis (eDNA & AIV)

### 9.1. Vertebrate Metabarcoding (OBITools/VSEARCH)

*Pipeline:* Demultiplexing, primer trimming, and OTU clustering.

```bash
obimultiplex -t [tag_file] -u [unidentified.fastq] [input.fastq] > [demultiplexed.fastq]
cutadapt -g [F_primer] -a [R_primer] -o [trimmed.fastq] [demultiplexed.fastq]
vsearch --cluster_size [trimmed.fasta] --id 0.97 --centroids [otus.fasta] --uchime_denovo [otus.fasta] --nonchimeras [otus_clean.fasta]
```

### 9.2. AIV Consensus Generation

*Pipeline:* Reference-based alignment and consensus calling.

```bash
minimap2 -ax map-ont [reference_segment.fasta] [filtered_aiv.fastq] | samtools sort > [aligned.bam]
bcftools mpileup -f [reference_segment.fasta] [aligned.bam] | bcftools call -c | vcfutils.pl vcf2fq > [consensus.fastq]
```
