# Denovo Transcriptome QC and Quantification Pipeline

A Nextflow pipeline for quality control and quantification of de novo transcriptome assemblies.

## Features

- **Flexible Input:** Supports short reads, long reads, hybrid assemblies, or assembly-only analysis.
- **Comprehensive QC:** Includes Transrate, BUSCO, and SQANTI3 for thorough assembly evaluation.
- **Quantification:** Performs quantification using Salmon (with BWA or minimap2 alignment) and Oarfish. Supports both merged and individual sample analysis.
- **Clustering:** Uses Corset to cluster transcripts into gene-level groups.
- **Customizable:** Easily configurable through parameters and a `nextflow.config` file.

## Quick Start

### 1. Assembly-only analysis
This mode runs only the QC steps on the provided assembly (no read quantification inputs provided).
```bash
nextflow run main.nf \
    --assembly /path/to/transcripts.fasta \
    --reference_fasta /path/to/reference.fasta \
    --reference_gtf /path/to/reference.gtf \
    --reference_genome /path/to/genome.fasta \
    --busco_lineage /path/to/busco_lineage \
    --input_type hybrid \
    --outdir results
```

### 2. Hybrid analysis with individual sample quantification
This example runs the full pipeline with both short and long reads, performing quantification for each sample provided in the list files.
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

## Download Example FASTA Files (Zenodo)

Example FASTA files are available from Zenodo record [17538009](https://zenodo.org/records/17538009).

### Option 1: Download manually
1. Open [https://zenodo.org/records/17538009](https://zenodo.org/records/17538009)
2. Download the `.tar.gz` archive from the **Files** section.

### Option 2: Download with command line
The example FASTA files are packaged in a single `.tar.gz` archive on Zenodo.

```bash
# Optional: inspect files available in the record
curl -s https://zenodo.org/api/records/17538009 | jq -r '.files[] | .key'

# Download the tar.gz archive (replace <archive_name.tar.gz> with the actual filename)
wget "https://zenodo.org/records/17538009/files/<archive_name.tar.gz>?download=1" -O <archive_name.tar.gz>

# Extract FASTA files
tar -xzf <archive_name.tar.gz>
```

After extraction, use the local FASTA path for the `--assembly` parameter.

## Download Reference Files

For human runs, use a single GENCODE release for `--reference_fasta`, `--reference_gtf`, and `--reference_genome` so sequence names stay consistent.

### 1. Reference transcript FASTA (`--reference_fasta`)
- Source: [GENCODE Human releases](https://www.gencodegenes.org/human/)
- Typical filename: `gencode.v<RELEASE>.transcripts.fa.gz`

### 2. Reference annotation GTF (`--reference_gtf`)
- Source: [GENCODE Human releases](https://www.gencodegenes.org/human/)
- Typical filename: `gencode.v<RELEASE>.annotation.gtf.gz`

### 3. Reference genome FASTA (`--reference_genome`, GRCh38)
- Source: [GENCODE Human releases](https://www.gencodegenes.org/human/)
- Typical filename: `GRCh38.primary_assembly.genome.fa.gz` (or `GRCh38.p14.genome.fa.gz` for full assembly files in newer releases)

### 4. BUSCO lineage dataset (`--busco_lineage`)
- Source: [BUSCO user guide](https://busco.ezlab.org/busco_userguide)
- List datasets:
```bash
busco --list-datasets
```
- Download a lineage dataset:
```bash
busco --download primates_odb10 --download_path busco_downloads
```
- Use it in this pipeline as:
  - `--busco_lineage busco_downloads/lineages/primates_odb10`

### Example download commands (GENCODE FTP pattern)
Replace `<RELEASE>` with your chosen release number (for example, `47`):

```bash
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_<RELEASE>/gencode.v<RELEASE>.transcripts.fa.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_<RELEASE>/gencode.v<RELEASE>.annotation.gtf.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_<RELEASE>/GRCh38.primary_assembly.genome.fa.gz
```

## Parameters

### Required Parameters
| Parameter | Description | Example |
|---|---|---|
| `--assembly` | Path to the transcriptome assembly file (FASTA). | `/path/to/zenodo_extracted/example_assembly.fasta` |
| `--reference_fasta` | Path to the reference transcriptome file (FASTA). | `gencode.v47.transcripts.fa.gz` |
| `--reference_gtf` | Path to the reference annotation file (GTF). | `gencode.v47.annotation.gtf.gz` |
| `--reference_genome` | Path to the reference genome file (FASTA). | `GRCh38.primary_assembly.genome.fa.gz` |
| `--busco_lineage` | Path to the BUSCO lineage database. | `busco_downloads/lineages/primates_odb10` |

### Input Options
| Parameter | Description | Default |
|---|---|---|
| `--input_type` | Type of input data. Options: `short`, `long`, `hybrid`. | `hybrid` |
| `--single_end` | Whether the short reads are single-end. | `true` |
| `--stranded` | Whether the data is strand-specific. This affects the SQANTI3 analysis. | `false` |
| `--long_read_tech` | Long read technology. Options: `ont`, `hifi`. | `ont` |
| `--short_reads` | Path to the merged short reads file (R1 for paired-end) for merged quantification. | `null` |
| `--short_reads2` | Path to the merged short reads file (R2 for paired-end) for merged quantification. | `null` |
| `--long_reads` | Path to the merged long reads file for merged quantification. | `null` |
| `--short_list` | Path to a text file containing a list of short read files for individual quantification. | `null` |
| `--long_list` | Path to a text file containing a list of long read files for individual quantification. | `null` |
| `--short_suffix1` | Suffix for R1 short read files to be removed for sample naming. | `_hybrid_sub_1.fastq.gz` |
| `--short_suffix2` | Suffix for R2 short read files. | `_hybrid_sub_2.fastq.gz` |
| `--long_suffix` | Suffix for long read files to be removed for sample naming. | `_hybrid_LR.fq.gz` |
| `--short_rep` | Number of short read replicates. | `12` |
| `--long_rep` | Number of long read replicates. | `12` |

### General Options
| Parameter | Description | Default |
|---|---|---|
| `--outdir` | The output directory. | `results` |
| `--cpus` | Number of CPUs to use for each process. | `48` |

### Configuration Parameters
These parameters define the paths to conda environments and tool installations. They are set in `nextflow.config` and may need to be updated to match your system.
| Parameter | Description |
|---|---|
| `transrate_env` | Path to the Transrate conda environment. |
| `busco_env` | Path to the BUSCO conda environment. |
| `sqanti_env` | Path to the SQANTI3 conda environment. |
| `isoncorrect_env` | Path to the isoncorrect conda environment. |
| `sqanti_path` | Path to the SQANTI3 installation directory. |
| `cupcake_path` | Path to the cDNA_Cupcake installation directory. |
| `seqtk_path` | Path to the seqtk executable. |
| `corset_path` | Path to the Corset executable. |
| `oarfish_path` | Path to the Oarfish executable. |


## Pipeline Steps

1.  **REFORMAT_SEQUENCES:** Converts 'U' to 'T' in the assembly using `seqtk`.
2.  **TRANSRATE:** Assesses the quality of the assembly against a reference transcriptome.
3.  **BUSCO:** Assesses the completeness of the assembly against a specified lineage.
4.  **SQANTI3:** Structurally annotates the assembly against a reference genome and annotation. For unstranded data, it performs a second pass to correct antisense transcripts.
5.  **BUSCO_CORRECTED:** Assesses the completeness of the SQANTI3-corrected assembly.
6.  **CORSET:** Clusters transcripts into gene-level groups based on self-alignment with `minimap2`.
7.  **SHORT_READ_QUANT (Optional):** Quantifies merged short reads against the assembly using `bwa` and `salmon`. Triggered by `--short_reads`.
8.  **LONG_READ_QUANT (Optional):** Quantifies merged long reads against the assembly using `minimap2`, `salmon`, and `oarfish`. Triggered by `--long_reads`.
9.  **INDIVIDUAL_SHORT_QUANT (Optional):** Quantifies individual short read samples from a list file. Triggered by `--short_list`.
10. **INDIVIDUAL_LONG_QUANT (Optional):** Quantifies individual long read samples from a list file. Triggered by `--long_list`.

## Output Structure

```
<outdir>/
├── 01_reformat/
├── 02_transrate/
├── 03_busco/
├── 04_sqanti3/
├── 05_busco_corrected/
├── 06_salmon_short/ (if --short_reads is provided)
├── 07_salmon_long/ (if --long_reads is provided)
├── 08_corset/
├── 09_dge_short/ (if --short_list is provided)
├── 10_dge_long/ (if --long_list is provided)
└── pipeline_info/
```

## Sample List File Format

The sample list files (`--short_list` and `--long_list`) should contain one file path per line. The pipeline will generate sample names by removing the specified suffixes (`--short_suffix1` or `--long_suffix`) and adding a prefix (`SR_` for short reads, `LR_` for long reads).

**Example `short_read_samples.txt`:**
```
/path/to/sample1_hybrid_sub_1.fastq.gz
/path/to/sample2_hybrid_sub_1.fastq.gz
```
This will generate samples named `SR_sample1` and `SR_sample2`.

**Example `long_read_samples.txt`:**
```
/path/to/sampleA_hybrid_LR.fq.gz
/path/to/sampleB_hybrid_LR.fq.gz
```
This will generate samples named `LR_sampleA` and `LR_sampleB`.

## Configuration

The pipeline relies on several conda environments and tool installations. The paths to these are defined as parameters in the `nextflow.config` file. Before running the pipeline, you should copy the `nextflow.config` file and modify the paths to match your system's configuration.

Example configuration for tool paths in `nextflow.config`:
```nextflow
params {
    // Your existing conda environments
    transrate_env = '/path/to/your/conda_env/transrate'
    busco_env = '/path/to/your/conda_env/busco'
    sqanti_env = '/path/to/your/conda_env/squanti3'
    isoncorrect_env = '/path/to/your/conda_env/isoncorrect'
    
    // Tool paths (update these to your actual paths)
    sqanti_path = '/path/to/your/SQANTI3-5.1.2/'
    cupcake_path = '/path/to/your/cDNA_Cupcake/'
    seqtk_path = '/path/to/your/seqtk/seqtk'
    corset_path = '/path/to/your/corset-1.09-linux64/corset'
    oarfish_path = '/path/to/your/oarfish'
}
```
