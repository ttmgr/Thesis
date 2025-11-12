# Bioinformatic Appendices

This document details the computational pipelines utilized for the analysis of environmental nanopore sequencing data. It is organized by publication and analytical stage. All commands are presented as executable one-liners to facilitate reproducibility and adherence to FAIR data principles.

## Appendix A: Air Microbiome Surveillance (Publication I)

**Overview:** A tailored pipeline for ultra-low-biomass bioaerosol samples, transitioning from high-accuracy basecalling to de novo assembly and taxonomic classification.

### A.1. Basecalling

**1.1. Controlled and Natural Environments (Guppy)**
*Model:* High Accuracy (HAC) for R10.4.1 flow cells.

```bash
guppy_basecaller -i [input_raw_data_dir] -r -s [output_dir] --detect_barcodes -c dna_r10.4.1_e8.2_400bps_hac.cfg -x "cuda:0"
```

**1.2. Urban Environment (Dorado)**
*Model:* Super Accuracy (SUP) with modified base retention (4mC/5mC).

```bash
dorado basecaller dna_r10.4.1_e8.2_400bps_sup@v5.0.0 [input_pod5_dir] -r --kit-name SQK-RBK114-24 --no-trim --emit-fastq > [basecalled.fastq]
```

**1.3. Demultiplexing (Dorado)**

```bash
dorado demux --output-dir [output_demux_dir] --kit-name SQK-RBK114-24 [basecalled.fastq]
```

### A.2. Read Pre-processing

**2.1. Adapter Trimming (Porechop)**
*Purpose:* Removal of sequencing adapters and barcodes.

```bash
porechop -i [input_barcode.fastq] -o [output_trimmed.fastq] -t 10
```

**2.2. Quality and Length Filtering (NanoFilt)**
*Thresholds:* Minimum length 100 bp.

```bash
cat [input_trimmed.fastq] | NanoFilt -l 100 -q 9 > [output_filtered.fastq]
```

### A.3. Metagenomic Assembly and Polishing

**3.1. De Novo Assembly (MetaFlye)**
*Strategy:* Long-read metagenomic assembly using the Nano-HQ mode.

```bash
flye --meta --nano-hq [input_filtered.fastq] --threads [threads] -o [output_assembly_dir]
```

**3.2. Read Mapping for Polishing (Minimap2)**
*Purpose:* Aligning filtered reads to the draft assembly.

```bash
minimap2 -ax map-ont -t [threads] [assembly.fasta] [input_filtered.fastq] > [alignment.sam]
```

**3.3. Alignment Sorting (Samtools)**

```bash
samtools view -b -@ [threads] [alignment.sam] | samtools sort -@ [threads] -o [sorted_alignment.bam]
```

**3.4. Assembly Polishing (Racon)**
*Purpose:* Consensus correction of the draft assembly.

```bash
racon -t [threads] [input_filtered.fastq] [alignment.sam] [assembly.fasta] > [polished_assembly.fasta]
```

### A.4. Taxonomic Classification

**4.1. Read-Level Classification (Kraken2)**
*Database:* NCBI nt database with memory mapping enabled.

```bash
kraken2 --db [kraken_db_path] --use-names --report [report_read.txt] --output [output_read.txt] [input_filtered.fastq] --memory-mapping --threads 28
```

**4.2. Contig-Level Classification (Kraken2)**

```bash
kraken2 --db [kraken_db_path] --use-names --report [report_contig.txt] --output [output_contig.txt] [polished_assembly.fasta] --memory-mapping --threads 28
```

## Appendix B: Wetland Ecosystem Surveillance (Publication II)

**Overview:** An integrated multi-omic framework processing shotgun metagenomics, RNA viromics, and targeted amplicons (eDNA/AIV) from passive water samplers.

### B.1. Basecalling and Pre-processing

**1.1. Basecalling (Dorado)**
*Model:* Super Accuracy (SUP) v5.0.0.

```bash
dorado basecaller dna_r10.4.1_e8.2_400bps_sup@v5.0.0 [input_pod5_dir] -r --kit-name SQK-RBK114-24 --no-trim --emit-fastq > [basecalled.fastq]
```

**1.2. Quality Filtering (NanoFilt)**
*Standard Metagenomics/Virome:* Length \> 100 bp, Q-score \> 9.

```bash
cat [input.fastq] | NanoFilt -l 100 -q 9 > [filtered_metagenomics.fastq]
```

*Targeted AIV Sequencing:* Relaxed length (\> 150 bp) and Q-score (\> 8).

```bash
cat [input_aiv.fastq] | NanoFilt -l 150 -q 8 > [filtered_aiv.fastq]
```

### B.2. Taxonomic Classification

**2.1. Metagenomic Profiling (Kraken2)**

```bash
kraken2 --db [nt_core_db] --threads [threads] --output [output.kraken] --report [report.txt] [filtered_metagenomics.fastq]
```

**2.2. Normalization (SeqKit)**
*Threshold:* Downsampling to 87,000 reads.

```bash
seqkit sample -n 87000 -s 100 [input.fastq] > [normalized.fastq]
```

### B.3. Metagenomic Assembly and Polishing

**3.1. Assembly Strategy 1 (metaFlye)**

```bash
flye --nano-hq [filtered_metagenomics.fastq] --out-dir [output_dir] --threads [threads] --meta
```

**3.2. Assembly Strategy 2 (nanoMDBG)**
*Note:* Optimized De Bruijn graph assembly for fragmented environmental DNA.

```bash
nanoMDBG [filtered_metagenomics.fastq] [k-mer_size] [output_prefix]
```

**3.3. Polishing (Medaka)**
*Purpose:* Final consensus polishing using the specific R10.4.1 model.

```bash
medaka_consensus -i [input_reads.fastq] -d [assembly_draft.fasta] -o [output_dir] -m r1041_e82_400bps_sup_v5.0.0
```

### B.4. Functional Annotation (Pathogen & AMR)

**4.1. Antimicrobial Resistance Detection (AMRFinderPlus)**
*Mode:* "Plus" enabled for stress response and virulence genes.

```bash
amrfinder -n [input.fasta] --plus --threads [threads] > [amr_report.tsv]
```

**4.2. Virulence Factor Detection (DIAMOND)**
*Target:* Virulence Factor Database (VFDB) core proteins (e.g., *ctxA/B*).

```bash
diamond blastx -d [vfdb_core.dmnd] -q [input_reads.fasta] -o [matches.tsv] -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
```

### B.5. Targeted Analysis (eDNA & AIV)

**5.1. Vertebrate Metabarcoding (OBITools/VSEARCH)**
*Pipeline:* Demultiplexing, primer trimming, and OTU clustering.

```bash
obimultiplex -t [tag_file] -u [unidentified.fastq] [input.fastq] > [demultiplexed.fastq]
cutadapt -g [F_primer] -a [R_primer] -o [trimmed.fastq] [demultiplexed.fastq]
vsearch --cluster_size [trimmed.fasta] --id 0.97 --centroids [otus.fasta] --uchime_denovo [otus.fasta] --nonchimeras [otus_clean.fasta]
```

**5.2. AIV Consensus Generation**

```bash
minimap2 -ax map-ont [reference_segment.fasta] [filtered_aiv.fastq] | samtools sort > [aligned.bam]
bcftools mpileup -f [reference_segment.fasta] [aligned.bam] | bcftools call -c | vcfutils.pl vcf2fq > [consensus.fastq]
```

### B.6. Pathogen Detection & MEGAN Post-Processing

**6.1. FASTQ to FASTA Conversion (Seqtk)**

```bash
seqtk seq -a [input.fastq] > [raw.fasta]
```

**6.2. Read Sorting (SeqKit)**

```bash
seqkit sort -n -w 0 --quiet [raw.fasta] -o [sorted.fasta]
```

**6.3. Read-Level Alignment (Minimap2)**
*Purpose:* Mapping reads to the NCBI-NT MMI index with high stringency.

```bash
minimap2 -ax map-ont -k 19 -w 10 -I 10G -g 5000 -r 2000 -N 100 --lj-min-ratio 0.5 -A 2 -B 5 -O 5,56 -E 4,1 -z 400,50 --sam-hit-only -t [threads] --split-prefix [temp_idx] [minimap2_db_mmi] [sorted.fasta] > [aligned.sam]
```

**6.4. SAM to RMA Conversion (MEGAN6)**
*Filtered (Stringent):*

```bash
sam2rma -i [aligned.sam] -r [sorted.fasta] -o [filtered.rma] -lg -alg longReads -t [threads] -mdb [megan_db_nucl] -ram readCount --minSupportPercent 0.01
```

**6.5. Taxonomic Information Extraction (rma2info)**
*Purpose:* Extracting assignments from Read-Level RMA files.

```bash
rma2info -i [filtered.rma] -o [taxonomy.r2c.txt] -r2c Taxonomy -n
rma2info -i [filtered.rma] -c2c Taxonomy -n -r -o [taxonomy.c2c.txt]
```

**6.6. Report Conversion**

```bash
python3 Convert_MEGAN_RMA_NCBI_c2c-snake.py --input [taxonomy.c2c.txt] --outname1 [names.temp.txt] --outname2 [codes.temp.txt] --mpa [report.mpa.txt] --kreport [report.kreport] --readsfile [read_counts.txt]
```

**6.7. Assembly-Level Pathogen Detection (metaMDBG)**
*Purpose:* Stringent pathogen identification applied directly to De Bruijn Graph assemblies.

**Step 1: Decompression**

```bash
gunzip -c [mdbg_contigs.fasta.gz] > [mdbg_contigs.fasta]
```

**Step 2: Contig Alignment (Minimap2)**

```bash
minimap2 -ax map-ont -k 19 -w 10 -I 10G -g 5000 -r 2000 -N 100 --lj-min-ratio 0.5 -A 2 -B 5 -O 5,56 -E 4,1 -z 400,50 --sam-hit-only -t [threads] [minimap2_db_mmi] [mdbg_contigs.fasta] > [mdbg_aligned.sam]
```

**Step 3: RMA Conversion (MEGAN6)**

```bash
sam2rma -i [mdbg_aligned.sam] -r [mdbg_contigs.fasta] -o [mdbg.rma] -lg -alg longReads -t [threads] -mdb [megan_db_nucl] -ram readCount --minSupportPercent 0.01
```

**Step 4: Taxonomic Information Extraction (rma2info)**

```bash
rma2info -i [mdbg.rma] -o [mdbg.r2c.txt] -r2c Taxonomy -n
rma2info -i [mdbg.rma] -c2c Taxonomy -n -r -o [mdbg.c2c.txt]
```

**Step 5: Report Conversion**

```bash
python3 Convert_MEGAN_RMA_NCBI_c2c-snake.py --input [mdbg.c2c.txt] --outname1 [names.temp.txt] --outname2 [codes.temp.txt] --mpa [report.mpa.txt] --kreport [report.kreport] --readsfile [contig_counts.txt]
```

## B.7. Genome Annotation & Gene Prediction

**7.1. Gene Prediction (Prodigal)**
*Purpose:* Predicting open reading frames (ORFs) on metaMDBG assemblies.

```bash
# Note: Uses the decompressed fasta from Step 6.7
prodigal -i [mdbg_contigs.fasta] -o [output.gff] -a [output.faa] -p meta -f gff
```

**7.2. Functional Annotation (Prokka)**
*Purpose:* Rapid annotation of prokaryotic genomes. Applied to both Racon-polished and MDBG assemblies.

```bash
prokka --force --quiet --outdir [output_dir] --prefix [sample_id] --cpus [threads] [input_assembly.fasta]
```
