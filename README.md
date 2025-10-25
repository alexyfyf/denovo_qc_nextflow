# Simple Transcriptome Analysis Pipeline

A simplified Nextflow pipeline based on Alex Yan's original bash script. This pipeline uses your existing local conda environments and includes **individual sample quantification** from list files.

## Features

✅ **Uses your existing conda environments** - No need to install new tools  
✅ **Simple structure** - Easy to understand and modify  
✅ **Flexible input types** - Short reads, long reads, hybrid, or assembly-only  
✅ **Individual sample quantification** - Process multiple samples from list files  
✅ **Proper paired-end support** - Handles both single-end and paired-end reads correctly  
✅ **All original functionality** - Transrate, BUSCO, SQANTI3, quantification, clustering  

## New Features Added

### 🔥 **Individual Sample Quantification**
- Process multiple samples from text files (like your `fq_SR.txt` and `fq_LR.txt`)
- Proper handling of file naming patterns with suffixes
- Separate quantification for each sample

### 🔥 **Proper Paired-End Support**
- Correctly handles single-end vs paired-end reads in Salmon
- Automatic detection of read2 files based on suffix patterns
- Proper Salmon parameters for each read type

## Quick Start

### 1. Basic Usage (Assembly-only analysis)
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

### 2. Individual Sample Quantification (like your original script)
```bash
nextflow run main.nf \
    --assembly /path/to/transcripts.fasta \
    --reference_fasta /path/to/reference.fasta \
    --reference_gtf /path/to/reference.gtf \
    --reference_genome /path/to/genome.fasta \
    --busco_lineage /path/to/busco_lineage \
    --input_type hybrid \
    --short_list /path/to/fq_SR.txt \
    --long_list /path/to/fq_LR.txt \
    --single_end true \
    --stranded false \
    --outdir results
```

### 3. Merged Quantification with Paired-End Reads
```bash
nextflow run main.nf \
    --assembly /path/to/transcripts.fasta \
    --reference_fasta /path/to/reference.fasta \
    --reference_gtf /path/to/reference.gtf \
    --reference_genome /path/to/genome.fasta \
    --busco_lineage /path/to/busco_lineage \
    --input_type short \
    --short_reads /path/to/reads_R1.fastq.gz \
    --short_reads2 /path/to/reads_R2.fastq.gz \
    --single_end false \
    --outdir results
```

### 3. Using your pea configuration
```bash
nextflow run main.nf -c pea_example.config -profile slurm
```

## Parameters

### Required Parameters
| Parameter | Description | Example |
|-----------|-------------|---------|
| `--assembly` | Transcriptome assembly (FASTA) | `/path/to/transcripts.fasta` |
| `--reference_fasta` | Reference transcriptome | `/path/to/reference.fasta` |
| `--reference_gtf` | Reference annotation | `/path/to/reference.gtf` |
| `--reference_genome` | Reference genome | `/path/to/genome.fasta` |
| `--busco_lineage` | BUSCO lineage database | `/path/to/lineage_odb10` |

### Analysis Options
| Parameter | Description | Default | Options |
|-----------|-------------|---------|---------|
| `--input_type` | Analysis type | `hybrid` | `short`, `long`, `hybrid`, `assembly_only` |
| `--single_end` | Single-end reads | `true` | `true`, `false` |
| `--stranded` | Strand-specific data | `false` | `true`, `false` |

### Individual Sample Quantification (NEW!)
| Parameter | Description | Example |
|-----------|-------------|---------|
| `--short_list` | Text file with short read file paths | `/path/to/fq_SR.txt` |
| `--long_list` | Text file with long read file paths | `/path/to/fq_LR.txt` |
| `--short_suffix1` | Suffix for read1 files | `_hybrid_sub_1.fastq.gz` |
| `--short_suffix2` | Suffix for read2 files (paired-end) | `_hybrid_sub_2.fastq.gz` |
| `--long_suffix` | Suffix for long read files | `_hybrid_LR.fq.gz` |

### Optional Merged Quantification
| Parameter | Description | Required for |
|-----------|-------------|--------------|
| `--short_reads` | Short read file (R1 for paired-end) | Merged short read analysis |
| `--short_reads2` | Short read file R2 (paired-end only) | Paired-end merged analysis |
| `--long_reads` | Single long read file | Merged long read analysis |

## Pipeline Steps

1. **Sequence Reformatting** - Convert U→T, prepare sequences
2. **Transrate** - Quality assessment
3. **BUSCO** - Completeness assessment  
4. **SQANTI3** - Structural annotation (handles stranded/unstranded)
5. **BUSCO Corrected** - Assessment of corrected transcriptome
6. **Corset** - Transcript clustering
7. **Merged Quantification** - BWA + Salmon (if single files provided)
8. **Individual Sample Quantification** - Process each sample from list files

## Output Structure

```
results/
├── 01_reformat/           # Reformatted sequences
├── 02_transrate/          # Transrate results
├── 03_busco/              # BUSCO results
├── 04_sqanti3/            # SQANTI3 annotation
├── 05_busco_corrected/    # BUSCO on corrected transcriptome
├── 06_salmon_short/       # Merged short read quantification
├── 07_salmon_long/        # Merged long read quantification
├── 08_corset/             # Transcript clustering
├── 09_dge_short/          # Individual short read samples
├── 10_dge_long/           # Individual long read samples
└── pipeline_info/         # Execution reports
```

## Sample List File Format

Your list files (`fq_SR.txt`, `fq_LR.txt`) should contain one file path per line:

**fq_SR.txt:**
```
/path/to/sample1_hybrid_sub_1.fastq.gz
/path/to/sample2_hybrid_sub_1.fastq.gz
/path/to/sample3_hybrid_sub_1.fastq.gz
```

**fq_LR.txt:**
```
/path/to/sample1_hybrid_LR.fq.gz
/path/to/sample2_hybrid_LR.fq.gz
/path/to/sample3_hybrid_LR.fq.gz
```

## Paired-End vs Single-End Handling

### Single-End Reads (`--single_end true`)
```bash
salmon quant -i index -l A -r reads.fastq -o output
```

### Paired-End Reads (`--single_end false`)
```bash
salmon quant -i index -l A -1 reads_1.fastq -2 reads_2.fastq -o output
```

The pipeline automatically:
- Detects read2 files by replacing `suffix1` with `suffix2`
- Uses appropriate Salmon parameters for each read type
- Handles BWA alignment for both single and paired reads

## Example for Your Pea Data

Based on your original script with individual sample processing:

```bash
nextflow run main.nf \
    --assembly /vast/projects/lab_davidson/yan.a/pea_fastq_merged/hybrid_merged/transcripts.fasta \
    --reference_fasta /vast/projects/lab_davidson/yan.a/ref/ncbi/pea_zw6/GCF_024323335.1_CAAS_Psat_ZW6_1.0_rna.fna \
    --reference_gtf /home/users/allstaff/yan.a/lab_davidson/yan.a/ref/ncbi/pea_zw6/GCF_024323335.1_CAAS_Psat_ZW6_1.0_genomic.gtf \
    --reference_genome /home/users/allstaff/yan.a/lab_davidson/yan.a/ref/ncbi/pea_zw6/GCF_024323335.1_CAAS_Psat_ZW6_1.0_genomic.fna \
    --busco_lineage /vast/projects/lab_davidson/yan.a/pea_fastq_merged/rnabloom2/busco_downloads/lineages/fabales_odb10 \
    --input_type hybrid \
    --short_list /vast/projects/lab_davidson/yan.a/pea_fastq_merged/hybrid_reads/fq_SR.txt \
    --long_list /vast/projects/lab_davidson/yan.a/pea_fastq_merged/hybrid_reads/fq_LR.txt \
    --short_suffix1 "_hybrid_sub_1.fastq.gz" \
    --long_suffix "_hybrid_LR.fq.gz" \
    --single_end true \
    --stranded false \
    --outdir results_pea \
    -profile slurm
```

## What's Fixed

### ✅ **Individual Sample Processing**
- Added `INDIVIDUAL_SHORT_QUANT` and `INDIVIDUAL_LONG_QUANT` processes
- Reads sample lists from text files (like your `fq_SR.txt` and `fq_LR.txt`)
- Processes each sample individually with proper naming

### ✅ **Proper Paired-End Support**
- Correctly detects read1 and read2 files based on suffix patterns
- Uses appropriate Salmon commands for single-end vs paired-end
- Handles BWA alignment for both read types

### ✅ **File Naming Logic**
- Extracts sample names by removing suffixes
- Automatically finds paired files for paired-end reads
- Maintains consistent naming throughout the pipeline

## Troubleshooting

1. **Update paths** in `nextflow.config` to match your system
2. **Check sample list files** contain correct file paths
3. **Verify suffix patterns** match your actual file names
4. **Check conda environments** are accessible
5. **Verify file paths** exist and are readable

This pipeline now fully replicates your original script functionality with the added benefits of Nextflow's workflow management!