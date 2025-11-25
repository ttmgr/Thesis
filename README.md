# Bioinformatic Appendices

This document provides a comprehensive specification of the computational pipelines established for the real-time genomic surveillance of aerial and aquatic microbiomes. It consolidates the workflows validated in Reska et al. (2024) and Perlas et al. (2025). All processing steps are presented as executable command-line instructions to facilitate the replication of these analyses in diverse environmental monitoring contexts.

## Appendix A: Air Microbiome Surveillance (Publication I)

**Overview:** A specialized pipeline optimized for the analysis of ultra-low-biomass bioaerosol samples[cite: 16], addressing the specific challenges of high DNA fragmentation and low input yields.

### A.1. Basecalling and Demultiplexing

**1.1. Controlled and Natural Environments (Guppy)**
*Model:* High Accuracy (HAC) model for R10.4.1 flow cells[cite: 83].

```bash
guppy_basecaller -i [input_raw_data_dir] -r -s [output_dir] --detect_barcodes -c dna_r10.4.1_e8.2_400bps_hac.cfg -x "cuda:0"
````

**1.2. Urban Environment (Dorado)**
[cite\_start]*Model:* High Accuracy (HAC) model (`dna_r10.4.1_e8.2_400bps_hac`) was used for validation[cite: 83].

```bash
dorado basecaller dna_r10.4.1_e8.2_400bps_hac@v4.3.0 [input_pod5_dir] -r --kit-name SQK-RBK114-24 --no-trim --emit-fastq > [basecalled.fastq]
```

**1.3. Demultiplexing (Dorado)**

```bash
dorado demux --output-dir [output_demux_dir] --kit-name SQK-RBK114-24 [basecalled.fastq]
```

### A.2. Read Pre-processing

**2.1. Adapter Trimming (Porechop)**
[cite\_start]*Purpose:* Removal of sequencing adapters and barcodes[cite: 85].

```bash
porechop -i [input_barcode.fastq] -o [output_trimmed.fastq] -t 10
```

**2.2. Quality and Length Filtering (NanoFilt)**
[cite\_start]*Thresholds:* Minimum length 100 bp; minimum average Q-score 8[cite: 85].

```bash
cat [input_trimmed.fastq] | NanoFilt -l 100 -q 8 > [output_filtered.fastq]
```

**2.3. Normalization (SeqKit)**
[cite\_start]*Purpose:* Downsampling reads for comparable taxonomic assessments (e.g., 30k reads for urban samples)[cite: 86].

```bash
seqkit sample -n 30000 -s 100 [output_filtered.fastq] > [output_normalized.fastq]
```

### A.3. Metagenomic Assembly and Polishing

**3.1. De Novo Assembly (MetaFlye)**
[cite\_start]*Strategy:* Long-read metagenomic assembly utilizing the Nano-HQ mode[cite: 89].

```bash
flye --meta --nano-hq [input_filtered.fastq] --threads [threads] -o [output_assembly_dir]
```

**3.2. Read Mapping for Polishing (Minimap2)**
[cite\_start]*Purpose:* Aligning filtered reads to the draft assembly for consensus correction[cite: 89].

```bash
minimap2 -ax map-ont -t [threads] [assembly.fasta] [input_filtered.fastq] > [alignment.sam]
```

**3.3. Assembly Polishing (Racon)**
[cite\_start]*Purpose:* Iterative consensus correction (3 rounds)[cite: 89].

```bash
racon -t [threads] [input_filtered.fastq] [alignment.sam] [assembly.fasta] > [polished_assembly.fasta]
```

### A.4. Taxonomic Classification

**4.1. Read-Level Classification (Kraken2)**
[cite\_start]*Database:* NCBI nt database with memory mapping enabled[cite: 86].

```bash
kraken2 --db [kraken_db_path] --use-names --report [report_read.txt] --output [output_read.txt] [output_normalized.fastq] --memory-mapping --threads 28
```

**4.2. Contig-Level Classification (Kraken2)**
[cite\_start]*Note:* Applied to bins or assembled contigs[cite: 515].

```bash
kraken2 --db [kraken_db_path] --use-names --report [report_contig.txt] --output [output_contig.txt] [polished_assembly.fasta] --memory-mapping --threads 28
```

### A.5. Binning & MAG Quality Control

**5.1. Metagenomic Binning (MetaWRAP)**
[cite\_start]*Purpose:* Integrating results from MetaBAT2, MaxBin2, and CONCOCT[cite: 90].

```bash
metawrap binning -o [output_dir] -t [threads] -a [assembly.fasta] --metabat2 --maxbin2 --concoct [clean_reads.fastq]
```

**5.2. Quality Assessment (CheckM)**
[cite\_start]*Thresholds:* Minimum completeness 30%, maximum contamination 10%[cite: 91, 92].

```bash
checkm lineage_wf -t [threads] -x fa [bin_directory] [output_directory]
```

### A.6. Functional Screening

**6.1. Mass Screening for Resistance/Virulence Genes (ABRicate/AMRFinderPlus)**
[cite\_start]*Purpose:* Screening of reads, contigs, and bins[cite: 95].

```bash
abricate --db [ncbi/card/vfdb] [input.fasta] > [output.tab]
amrfinder -n [input.fasta] --plus --threads [threads] > [output.amr]
```

-----

## Appendix B: Wetland Ecosystem Surveillance (Publication II)

[cite\_start]**Overview:** An integrated multi-omic framework capable of processing shotgun metagenomics, RNA viromics, and targeted amplicons (eDNA/AIV) from passive water samplers[cite: 743].

### B.1. Basecalling and Pre-processing

**1.1. Basecalling (Dorado)**
[cite\_start]*Model:* Super Accuracy (SUP) v5.0.0[cite: 1039].

```bash
dorado basecaller dna_r10.4.1_e8.2_400bps_sup@v5.0.0 [input_pod5_dir] -r --kit-name SQK-RBK114-24 --no-trim --emit-fastq > [basecalled.fastq]
```

**1.2. Adapter Trimming (Porechop)**
[cite\_start]*Purpose:* Demultiplexing and adapter/barcode trimming for shotgun and virome reads[cite: 1040].

```bash
porechop -i [basecalled.fastq] -o [trimmed.fastq] --threads [threads]
```

**1.3. Quality Filtering (NanoFilt)**
[cite\_start]*Standard Metagenomics/Virome:* Length \> 100 bp, Q-score \> 9[cite: 1041].

```bash
cat [trimmed.fastq] | NanoFilt -l 100 -q 9 > [filtered_metagenomics.fastq]
```

[cite\_start]*Targeted AIV Sequencing:* Relaxed length (\> 150 bp) and Q-score (\> 8)[cite: 1041].

```bash
cat [trimmed_aiv.fastq] | NanoFilt -l 150 -q 8 > [filtered_aiv.fastq]
```

### B.2. Taxonomic Classification

**2.1. Metagenomic Profiling (Kraken2)**
[cite\_start]*Database:* NCBI nt\_core[cite: 1050].

```bash
kraken2 --db [nt_core_db] --threads [threads] --output [output.kraken] --report [report.txt] [filtered_metagenomics.fastq]
```

**2.2. Normalization (SeqKit)**
[cite\_start]*Threshold:* Downsampling to 87,000 reads for comparative analysis[cite: 1051].

```bash
seqkit sample -n 87000 -s 100 [filtered_metagenomics.fastq] > [normalized.fastq]
```

### B.3. Metagenomic Assembly and Polishing

**3.1. Workflow A: metaFlye with Hybrid Polishing**
[cite\_start]*Pipeline:* metaFlye -\> Minimap2 -\> Racon (3 rounds) -\> Medaka [cite: 1059, 1061-1064].

```bash
# Assembly
flye --nano-hq [filtered_metagenomics.fastq] --out-dir [flye_out] --threads [threads] --meta

# Polishing Round 1-3 (Racon)
minimap2 -ax map-ont -t [threads] [flye_out/assembly.fasta] [filtered_metagenomics.fastq] > [aln.sam]
racon -t [threads] [filtered_metagenomics.fastq] [aln.sam] [flye_out/assembly.fasta] > [racon_1.fasta]
# (Repeat block 2 more times for 3 rounds total)

# Polishing Round 4 (Medaka)
medaka_consensus -i [filtered_metagenomics.fastq] -d [racon_3.fasta] -o [final_flye_dir] -m r1041_e82_400bps_sup_v5.0.0
```

**3.2. Workflow B: nanoMDBG with Medaka Polishing**
[cite\_start]*Pipeline:* nanoMDBG -\> Medaka only[cite: 1059, 1064].

```bash
# Assembly
nanoMDBG [filtered_metagenomics.fastq] [k-mer_size] [output_prefix]

# Polishing
medaka_consensus -i [filtered_metagenomics.fastq] -d [mdbg_contigs.fasta] -o [final_mdbg_dir] -m r1041_e82_400bps_sup_v5.0.0
```

### B.4. RNA Virome Analysis

**4.1. Viral Assembly (nanoMDBG)**
[cite\_start]*Note:* nanoMDBG was used for viral de novo assembly followed by Medaka polishing[cite: 1108].

**(See B.3.2 for commands)**

**4.2. Viral Taxonomy Assignment (DIAMOND BLASTx)**
[cite\_start]*Database:* NCBI non-redundant protein database (NR)[cite: 1109].
[cite\_start]*Threshold:* Contigs \>80% identity to kingdom "Viruses" (taxid: 10239)[cite: 1110].

```bash
diamond blastx -d [nr_db.dmnd] -q [viral_contigs.fasta] -o [viral_matches.tsv] -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore --sensitive
```

### B.5. Functional Annotation (Pathogen & AMR)

**5.1. Antimicrobial Resistance Detection (AMRFinderPlus)**
[cite\_start]*Mode:* "Plus" enabled for stress response and virulence genes; nucleotide and protein analysis [cite: 1077-1080].

```bash
amrfinder -n [input.fasta] --plus --threads [threads] > [amr_report.tsv]
```

**5.2. Virulence Factor Detection (DIAMOND)**
[cite\_start]*Target:* Virulence Factor Database (VFDB) core proteins (e.g., *ctxA/B*)[cite: 1074].

```bash
diamond blastx -d [vfdb_core.dmnd] -q [input_reads.fasta] -o [matches.tsv] -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
```

**5.3. General Functional Annotation (Prokka)**
[cite\_start]*Purpose:* Rapid annotation of prokaryotic genomes (on contigs)[cite: 1066].

```bash
prokka --force --quiet --outdir [output_dir] --prefix [sample_id] --cpus [threads] [polished_assembly.fasta]
```

### B.6. Targeted Analysis (eDNA & AIV)

**6.1. Vertebrate Metabarcoding (OBITools/VSEARCH)**
[cite\_start]*Pipeline:* Demultiplexing, primer trimming, and OTU clustering [cite: 1085-1089].

```bash
obimultiplex -t [tag_file] -u [unidentified.fastq] [input.fastq] > [demultiplexed.fastq]
cutadapt -g [F_primer] -a [R_primer] -o [trimmed.fastq] [demultiplexed.fastq]
vsearch --cluster_size [trimmed.fasta] --id 0.97 --centroids [otus.fasta] --uchime_denovo [otus.fasta] --nonchimeras [otus_clean.fasta]
```

**6.2. AIV Consensus Generation**
[cite\_start]*Pipeline:* Alignment to NCBI Influenza Virus Database segments [cite: 1123-1128].

```bash
minimap2 -ax map-ont [reference_segment.fasta] [filtered_aiv.fastq] | samtools sort > [aligned.bam]
bcftools mpileup -f [reference_segment.fasta] [aligned.bam] | bcftools call -c | vcfutils.pl vcf2fq > [consensus.fastq]
```

### B.7. Pathogen Detection & MEGAN Post-Processing

**7.1. FASTQ to FASTA Conversion (Seqtk)**
[cite\_start]*Purpose:* Preparation for alignment[cite: 96].

```bash
seqtk seq -a [input.fastq] > [raw.fasta]
```

**7.2. Read-Level Alignment (Minimap2)**
[cite\_start]*Purpose:* Mapping reads to the NCBI-NT MMI index with high stringency[cite: 1067].

```bash
minimap2 -ax map-ont -k 19 -w 10 -I 10G -g 5000 -r 2000 -N 100 --lj-min-ratio 0.5 -A 2 -B 5 -O 5,56 -E 4,1 -z 400,50 --sam-hit-only -t [threads] --split-prefix [temp_idx] [minimap2_db_mmi] [sorted.fasta] > [aligned.sam]
```

**7.3. SAM to RMA Conversion (MEGAN6)**
[cite\_start]*Purpose:* Lowest Common Ancestor (LCA) assignment[cite: 1067].
[cite\_start]*Criteria:* Taxonomic assignment accepted only if \>50% of near-best alignments match the same genus[cite: 1068].

```bash
sam2rma -i [aligned.sam] -r [sorted.fasta] -o [filtered.rma] -lg -alg longReads -t [threads] -mdb [megan_db_nucl] -ram readCount --minSupportPercent 0.01
```

**7.4. Taxonomic Information Extraction (rma2info)**

```bash
rma2info -i [filtered.rma] -o [taxonomy.r2c.txt] -r2c Taxonomy -n
rma2info -i [filtered.rma] -c2c Taxonomy -n -r -o [taxonomy.c2c.txt]
```

### B.8. Phylogenetic Analysis

**8.1. Multiple Sequence Alignment (MAFFT)**
[cite\_start]*Purpose:* Aligning the consensus H4 HA sequence with GISAID reference sequences[cite: 1143].

```bash
mafft --auto --thread [threads] [combined_sequences.fasta] > [alignment.aln]
```

**8.2. Phylogeny Inference (IQ-TREE2)**
[cite\_start]*Purpose:* Constructing Maximum-Likelihood trees with ultrafast bootstrap support[cite: 1143].

```bash
iqtree2 -s [alignment.aln] -m MFP -bb 1000 -alrt 1000 -nt [threads]
```

### B.9. Mobile Genetic Element Analysis

**9.1. Plasmid Detection (PlasmidFinder)**
[cite\_start]*Purpose:* Identifying plasmid replicons in assembled contigs[cite: 1083].

```bash
python3 plasmidfinder.py -i [input_assembly.fasta] -o [output_dir] -p [database_path]
```

## Appendix C: Web-Based Analytical Interfaces

In addition to command-line processing, the following web-based platforms were integral to the methodology:

  * [cite\_start]**CZID (Chan Zuckerberg ID):** Utilized for hybrid taxonomic classification benchmarking in Publication I [cite: 88] [cite\_start]and stringent pathogen species cross-referencing in Publication II[cite: 1073].
  * [cite\_start]**GISAID BLAST & FluSurver:** Employed for the subtyping and mutation analysis of the assembled Avian Influenza Virus (AIV) sequences [cite: 1129-1130].
  * [cite\_start]**ITOL (Interactive Tree of Life):** Used for the visualization and annotation of the AIV phylogenetic trees[cite: 1147].

<!-- end list -->

```
```
