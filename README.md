# Denovo Transcriptome QC and Quantification Pipeline

A Nextflow pipeline for quality control and quantification of de novo transcriptome assemblies.

## Features

- **Flexible Input:** Supports short reads, long reads, hybrid assemblies, or assembly-only analysis.
- **Comprehensive QC:** Includes Transrate, BUSCO, and SQANTI3 for thorough assembly evaluation.
- **Quantification:** Performs quantification using Salmon, with options for both merged and individual sample analysis.
- **Clustering:** Uses Corset to cluster transcripts into gene-level groups.
- **Customizable:** Easily configurable through parameters and a `nextflow.config` file.

## Quick Start

### 1. Assembly-only analysis
```bash
nextflow run main.nf \
    --assembly /path/to/transcripts.fasta \
    --reference_fasta /path/to/reference.fasta \
    --reference_gtf /path/to/reference.gtf \
    --reference_genome /path/to/genome.fasta \
    --busco_lineage /path/to/busco_lineage \
    --input_type assembly_only \
    --outdir results
```

### 2. Hybrid analysis with individual sample quantification
```bash
nextflow run main.nf \
    --assembly /path/to/transcripts.fasta \
    --reference_fasta /path/to/reference.fasta \
    --reference_gtf /path/to/reference.gtf \
    --reference_genome /path/to/genome.fasta \
    --busco_lineage /path/to/busco_lineage \
    --input_type hybrid \
    --short_list /path/to/short_read_samples.txt \
    --long_list /path/to/long_read_samples.txt \
    --outdir results
```

## Parameters

### Required Parameters
| Parameter | Description |
|---|---|
| `--assembly` | Path to the transcriptome assembly file (FASTA). |
| `--reference_fasta` | Path to the reference transcriptome file (FASTA). |
| `--reference_gtf` | Path to the reference annotation file (GTF). |
| `--reference_genome` | Path to the reference genome file (FASTA). |
| `--busco_lineage` | Path to the BUSCO lineage database. |

### Input Options
| Parameter | Description | Default |
|---|---|---|
| `--input_type` | Type of input data. Options: `short`, `long`, `hybrid`, `assembly_only`. | `hybrid` |
| `--single_end` | Whether the short reads are single-end. | `true` |
| `--stranded` | Whether the data is strand-specific. | `false` |
| `--long_read_tech` | Long read technology. Options: `ont`, `hifi`. | `ont` |
| `--short_reads` | Path to the merged short reads file (R1 for paired-end). | `null` |
| `--short_reads2` | Path to the merged short reads file (R2 for paired-end). | `null` |
| `--long_reads` | Path to the merged long reads file. | `null` |
| `--short_list` | Path to a text file containing a list of short read files for individual quantification. | `null` |
| `--long_list` | Path to a text file containing a list of long read files for individual quantification. | `null` |
| `--short_suffix1` | Suffix for R1 short read files. | `_hybrid_sub_1.fastq.gz` |
| `--short_suffix2` | Suffix for R2 short read files. | `_hybrid_sub_2.fastq.gz` |
| `--long_suffix` | Suffix for long read files. | `_hybrid_LR.fq.gz` |

### General Options
| Parameter | Description | Default |
|---|---|---|
| `--outdir` | The output directory. | `results` |
| `--cpus` | Number of CPUs to use for each process. | `48` |

## Pipeline Steps

1.  **REFORMAT_SEQUENCES:** Converts 'U' to 'T' in the assembly.
2.  **TRANSRATE:** Assesses the quality of the assembly.
3.  **BUSCO:** Assesses the completeness of the assembly.
4.  **SQANTI3:** Structurally annotates the assembly against a reference.
5.  **BUSCO_CORRECTED:** Assesses the completeness of the SQANTI3-corrected assembly.
6.  **SHORT_READ_QUANT:** Quantifies the merged short reads against the assembly.
7.  **LONG_READ_QUANT:** Quantifies the merged long reads against the assembly.
8.  **CORSET:** Clusters transcripts based on sequence similarity.
9.  **INDIVIDUAL_SHORT_QUANT:** Quantifies individual short read samples.
10. **INDIVIDUAL_LONG_QUANT:** Quantifies individual long read samples.

## Output Structure

```
<outdir>/
├── 01_reformat/
├── 02_transrate/
├── 03_busco/
├── 04_sqanti3/
├── 05_busco_corrected/
├── 06_salmon_short/
├── 07_salmon_long/
├── 08_corset/
├── 09_dge_short/
├── 10_dge_long/
└── pipeline_info/
```

## Sample List File Format

The sample list files (`--short_list` and `--long_list`) should contain one file path per line.

**Example `short_read_samples.txt`:**
```
/path/to/sample1_R1.fastq.gz
/path/to/sample2_R1.fastq.gz
```

**Example `long_read_samples.txt`:**
```
/path/to/sample1.fastq.gz
/path/to/sample2.fastq.gz
```

