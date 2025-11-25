# Bioinformatic Appendices

This document provides a comprehensive specification of the computational pipelines established for the cumulative dissertation titled **"Genomic Surveillance of Aerial and Aquatic Microbiomes by Nanopore Sequencing."** It consolidates the bioinformatic workflows validated in the constituent studies **Reska et al. (2024)** and **Perlas et al. (2025)**. All processing steps are presented as executable command-line instructions to ensure the reproducibility of these real-time genomic surveillance frameworks in diverse environmental monitoring contexts.

## Software Versions & Tools

The following table details the specific software versions employed in each study.

| Tool | Publication I (Air) | Publication II (Wetlands) | Purpose |
| :--- | :--- | :--- | :--- |
| **Basecalling** | | | |
| MinKNOW | v23.04.3 / v23.04.5 | v24.11.10 | Data acquisition & device control |
| Guppy | v6.3.2 | — | High-accuracy basecalling (R10.4.1) |
| Dorado | v4.3.0 | v5.0.0 | Super-accuracy basecalling (R10.4.1) |
| **Pre-processing** | | | |
| Porechop | v0.2.3 | v0.2.4 | Adapter and barcode trimming |
| NanoFilt | v2.8.0 | v2.8.0 | Read quality and length filtering |
| SeqKit | v2.8.2 | v2.3.0 | Read sampling, sorting, and formatting |
| **Taxonomy** | | | |
| Kraken2 | v2.0.7 | v2.1.2 | Metagenomic taxonomic classification |
| MEGAN-CE | — | v6.21.1 | Lowest Common Ancestor (LCA) analysis |
| OBITools4 | — | v1.3.1 | Metabarcoding demultiplexing |
| VSEARCH | — | v2.21 | OTU clustering and chimera removal |
| Cutadapt | — | v4.2 | Primer trimming (amplicon data) |
| **Assembly** | | | |
| MetaFlye | v2.9.1 | v2.9.6 | Long-read de novo assembly |
| nanoMDBG | — | v1.1 | De Bruijn graph assembly (low biomass) |
| Minimap2 | v2.17 | v2.28 | Read mapping / polishing alignment |
| Racon | v1.5 | v1.5.0 | Assembly polishing (Consensus) |
| Medaka | — | v1.7.2 | Assembly polishing (Neural network) |
| **Downstream** | | | |
| MetaWRAP | v1.3 | — | Metagenomic binning wrapper |
| CheckM | v1.2.2 | — | MAG quality assessment |
| AMRFinderPlus | v3.12.8 | v4.0.23 | Antimicrobial resistance gene detection |
| ABRicate | v1.0.1 | — | Mass screening of contigs |
| DIAMOND | — | v2.1.13 | Protein alignment (Virulence/Viral) |
| Prokka | — | v1.14.5 | Prokaryotic genome annotation |
| Prodigal | — | v2.6.3 | Gene prediction |
| PlasmidFinder | — | v2.1.6 | Plasmid detection |
| MAFFT | — | v7.526 | Multiple sequence alignment |
| IQ-TREE2 | — | v2.3.4 | Phylogeny inference |
| SAMtools | — | v1.17 | Alignment file processing |
| BCFtools | — | v1.17 | Variant calling / Consensus generation |

---

## Appendix A: Air Microbiome Surveillance (Publication I)

**Overview:** A specialized pipeline optimized for the analysis of ultra-low-biomass bioaerosol samples, addressing the specific challenges of high DNA fragmentation and low input yields through sensitive basecalling, rigorous assembly, and genome binning.

### A.1. Basecalling and Demultiplexing

**1.1. Controlled and Natural Environments (Guppy)**
*Model:* High Accuracy (HAC) model for R10.4.1 flow cells.

```bash
guppy_basecaller -i [input_raw_data_dir] -r -s [output_dir] --detect_barcodes -c dna_r10.4.1_e8.2_400bps_hac.cfg -x "cuda:0"
````

**1.2. Urban Environment (Dorado)**
*Model:* High Accuracy (HAC) model (`dna_r10.4.1_e8.2_400bps_hac`) was used for validation.

```bash
dorado basecaller dna_r10.4.1_e8.2_400bps_hac@v4.3.0 [input_pod5_dir] -r --kit-name SQK-RBK114-24 --no-trim --emit-fastq > [basecalled.fastq]
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
*Thresholds:* Minimum length 100 bp; minimum average Q-score 8.

```bash
cat [input_trimmed.fastq] | NanoFilt -l 100 -q 8 > [output_filtered.fastq]
```

**2.3. Normalization (SeqKit)**
*Purpose:* Downsampling reads for comparable taxonomic assessments (e.g., 30k reads for urban samples).

```bash
seqkit sample -n 30000 -s 100 [output_filtered.fastq] > [output_normalized.fastq]
```

### A.3. Metagenomic Assembly and Polishing

**3.1. De Novo Assembly (MetaFlye)**
*Strategy:* Long-read metagenomic assembly utilizing the Nano-HQ mode.

```bash
flye --meta --nano-hq [input_filtered.fastq] --threads [threads] -o [output_assembly_dir]
```

**3.2. Read Mapping for Polishing (Minimap2)**
*Purpose:* Aligning filtered reads to the draft assembly for consensus correction.

```bash
minimap2 -ax map-ont -t [threads] [assembly.fasta] [input_filtered.fastq] > [alignment.sam]
```

**3.3. Assembly Polishing (Racon)**
*Purpose:* Iterative consensus correction (3 rounds).

```bash
racon -t [threads] [input_filtered.fastq] [alignment.sam] [assembly.fasta] > [polished_assembly.fasta]
```

### A.4. Taxonomic Classification

**4.1. Read-Level Classification (Kraken2)**
*Database:* NCBI nt database with memory mapping enabled.

```bash
kraken2 --db [kraken_db_path] --use-names --report [report_read.txt] --output [output_read.txt] [output_normalized.fastq] --memory-mapping --threads 28
```

**4.2. Contig-Level Classification (Kraken2)**
*Note:* Applied to bins or assembled contigs.

```bash
kraken2 --db [kraken_db_path] --use-names --report [report_contig.txt] --output [output_contig.txt] [polished_assembly.fasta] --memory-mapping --threads 28
```

### A.5. Binning & MAG Quality Control

**5.1. Metagenomic Binning (MetaWRAP)**
*Purpose:* Integrating results from MetaBAT2, MaxBin2, and CONCOCT.

```bash
metawrap binning -o [output_dir] -t [threads] -a [assembly.fasta] --metabat2 --maxbin2 --concoct [clean_reads.fastq]
```

**5.2. Quality Assessment (CheckM)**
*Thresholds:* Minimum completeness 30%, maximum contamination 10%.

```bash
checkm lineage_wf -t [threads] -x fa [bin_directory] [output_directory]
```

### A.6. Functional Screening

**6.1. Mass Screening for Resistance/Virulence Genes (ABRicate/AMRFinderPlus)**
*Purpose:* Screening of reads, contigs, and bins.

```bash
abricate --db [ncbi/card/vfdb] [input.fasta] > [output.tab]
amrfinder -n [input.fasta] --plus --threads [threads] > [output.amr]
```

-----

## Appendix B: Wetland Ecosystem Surveillance (Publication II)

**Overview:** An integrated multi-omic framework capable of processing shotgun metagenomics, RNA viromics, and targeted amplicons (eDNA/AIV) from passive water samplers.

### B.1. Basecalling and Pre-processing

**1.1. Basecalling (Dorado)**
*Model:* Super Accuracy (SUP) v5.0.0.

```bash
dorado basecaller dna_r10.4.1_e8.2_400bps_sup@v5.0.0 [input_pod5_dir] -r --kit-name SQK-RBK114-24 --no-trim --emit-fastq > [basecalled.fastq]
```

**1.2. Adapter Trimming (Porechop)**
*Purpose:* Demultiplexing and adapter/barcode trimming for shotgun and virome reads.

```bash
porechop -i [basecalled.fastq] -o [trimmed.fastq] --threads [threads]
```

**1.3. Quality Filtering (NanoFilt)**
*Standard Metagenomics/Virome:* Length \> 100 bp, Q-score \> 9.

```bash
cat [trimmed.fastq] | NanoFilt -l 100 -q 9 > [filtered_metagenomics.fastq]
```

*Targeted AIV Sequencing:* Relaxed length (\> 150 bp) and Q-score (\> 8).

```bash
cat [trimmed_aiv.fastq] | NanoFilt -l 150 -q 8 > [filtered_aiv.fastq]
```

### B.2. Taxonomic Classification

**2.1. Metagenomic Profiling (Kraken2)**
*Database:* NCBI nt\_core.

```bash
kraken2 --db [nt_core_db] --threads [threads] --output [output.kraken] --report [report.txt] [filtered_metagenomics.fastq]
```

**2.2. Normalization (SeqKit)**
*Threshold:* Downsampling to 87,000 reads for comparative analysis.

```bash
seqkit sample -n 87000 -s 100 [filtered_metagenomics.fastq] > [normalized.fastq]
```

### B.3. Metagenomic Assembly and Polishing

**3.1. Workflow A: metaFlye with Hybrid Polishing**
*Pipeline:* metaFlye -\> Minimap2 -\> Racon (3 rounds) -\> Medaka.

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
*Pipeline:* nanoMDBG -\> Medaka only.

```bash
# Assembly
nanoMDBG [filtered_metagenomics.fastq] [k-mer_size] [output_prefix]

# Polishing
medaka_consensus -i [filtered_metagenomics.fastq] -d [mdbg_contigs.fasta] -o [final_mdbg_dir] -m r1041_e82_400bps_sup_v5.0.0
```

### B.4. RNA Virome Analysis

**4.1. Viral Assembly (nanoMDBG)**
*Note:* nanoMDBG was used for viral de novo assembly followed by Medaka polishing.

**(See B.3.2 for commands)**

**4.2. Viral Taxonomy Assignment (DIAMOND BLASTx)**
*Database:* NCBI non-redundant protein database (NR).
*Threshold:* Contigs \>80% identity to kingdom "Viruses" (taxid: 10239).

```bash
diamond blastx -d [nr_db.dmnd] -q [viral_contigs.fasta] -o [viral_matches.tsv] -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore --sensitive
```

### B.5. Functional Annotation (Pathogen & AMR)

**5.1. Antimicrobial Resistance Detection (AMRFinderPlus)**
*Mode:* "Plus" enabled for stress response and virulence genes; nucleotide and protein analysis.

```bash
amrfinder -n [input.fasta] --plus --threads [threads] > [amr_report.tsv]
```

**5.2. Virulence Factor Detection (DIAMOND)**
*Target:* Virulence Factor Database (VFDB) core proteins (e.g., *ctxA/B*).

```bash
diamond blastx -d [vfdb_core.dmnd] -q [input_reads.fasta] -o [matches.tsv] -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
```

**5.3. General Functional Annotation (Prokka)**
*Purpose:* Rapid annotation of prokaryotic genomes (on contigs).

```bash
prokka --force --quiet --outdir [output_dir] --prefix [sample_id] --cpus [threads] [polished_assembly.fasta]
```

### B.6. Targeted Analysis (eDNA & AIV)

**6.1. Vertebrate Metabarcoding (OBITools/VSEARCH)**
*Pipeline:* Demultiplexing, primer trimming, and OTU clustering.

```bash
obimultiplex -t [tag_file] -u [unidentified.fastq] [input.fastq] > [demultiplexed.fastq]
cutadapt -g [F_primer] -a [R_primer] -o [trimmed.fastq] [demultiplexed.fastq]
vsearch --cluster_size [trimmed.fasta] --id 0.97 --centroids [otus.fasta] --uchime_denovo [otus.fasta] --nonchimeras [otus_clean.fasta]
```

**6.2. AIV Consensus Generation**
*Pipeline:* Alignment to NCBI Influenza Virus Database segments.

```bash
minimap2 -ax map-ont [reference_segment.fasta] [filtered_aiv.fastq] | samtools sort > [aligned.bam]
bcftools mpileup -f [reference_segment.fasta] [aligned.bam] | bcftools call -c | vcfutils.pl vcf2fq > [consensus.fastq]
```

### B.7. Pathogen Detection & MEGAN Post-Processing

**7.1. FASTQ to FASTA Conversion (Seqtk)**
*Purpose:* Preparation for alignment.

```bash
seqtk seq -a [input.fastq] > [raw.fasta]
```

**7.2. Read-Level Alignment (Minimap2)**
*Purpose:* Mapping reads to the NCBI-NT MMI index with high stringency.

```bash
minimap2 -ax map-ont -k 19 -w 10 -I 10G -g 5000 -r 2000 -N 100 --lj-min-ratio 0.5 -A 2 -B 5 -O 5,56 -E 4,1 -z 400,50 --sam-hit-only -t [threads] --split-prefix [temp_idx] [minimap2_db_mmi] [sorted.fasta] > [aligned.sam]
```

**7.3. SAM to RMA Conversion (MEGAN6)**
*Purpose:* Lowest Common Ancestor (LCA) assignment.
*Criteria:* Taxonomic assignment accepted only if \>50% of near-best alignments match the same genus.

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
*Purpose:* Aligning the consensus H4 HA sequence with GISAID reference sequences.

```bash
mafft --auto --thread [threads] [combined_sequences.fasta] > [alignment.aln]
```

**8.2. Phylogeny Inference (IQ-TREE2)**
*Purpose:* Constructing Maximum-Likelihood trees with ultrafast bootstrap support.

```bash
iqtree2 -s [alignment.aln] -m MFP -bb 1000 -alrt 1000 -nt [threads]
```

### B.9. Mobile Genetic Element Analysis

**9.1. Plasmid Detection (PlasmidFinder)**
*Purpose:* Identifying plasmid replicons in assembled contigs.

```bash
python3 plasmidfinder.py -i [input_assembly.fasta] -o [output_dir] -p [database_path]
```

## Appendix C: Web-Based Analytical Interfaces

In addition to command-line processing, the following web-based platforms were integral to the methodology:

  * **CZID (Chan Zuckerberg ID):** Utilized for hybrid taxonomic classification benchmarking in Publication I and stringent pathogen species cross-referencing in Publication II.
  * **GISAID BLAST & FluSurver:** Employed for the subtyping and mutation analysis of the assembled Avian Influenza Virus (AIV) sequences.
  * **ITOL (Interactive Tree of Life):** Used for the visualization and annotation of the AIV phylogenetic trees.

<!-- end list -->

```
```
