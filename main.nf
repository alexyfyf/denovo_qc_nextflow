#!/usr/bin/env nextflow

/*
 * Simple Transcriptome Analysis Pipeline
 * Based on Alex Yan's original script
 * Uses local conda environments
 */

nextflow.enable.dsl = 2

// Parameters
params.assembly = null
params.reference_fasta = null
params.reference_gtf = null
params.reference_genome = null
params.busco_lineage = null

// Input type: 'short', 'long', 'hybrid'
params.input_type = 'hybrid'
params.single_end = true
params.stranded = false

// Long read technology: 'ont' or 'hifi'
params.long_read_tech = 'ont'

// Read files (optional - for merged quantification)
params.short_reads = null
params.short_reads2 = null  // For paired-end merged quantification
params.long_reads = null

// Sample lists for individual quantification
params.short_list = null  // Text file with list of short read files
params.long_list = null   // Text file with list of long read files
params.short_suffix1 = "_hybrid_sub_1.fastq.gz"
params.short_suffix2 = "_hybrid_sub_2.fastq.gz"
params.long_suffix = "_hybrid_LR.fq.gz"
params.short_rep = 12
params.long_rep = 12

// Output directory
params.outdir = 'results'

// Conda environments (use your existing ones)
params.transrate_env = '/stornext/Bioinf/data/lab_davidson/yan.a/conda_env/transrate'
params.busco_env = '/stornext/Bioinf/data/lab_davidson/yan.a/conda_env/busco'
params.sqanti_env = '/stornext/Bioinf/data/lab_davidson/yan.a/conda_env/squanti3'
params.isoncorrect_env = '/stornext/Bioinf/data/lab_davidson/yan.a/conda_env/isoncorrect'

// Tool paths
params.sqanti_path = '/vast/projects/lab_davidson/yan.a/pea_fastq_merged/SQANTI3-5.1.2/'
params.cupcake_path = '$HOME/yan.a/software/cDNA_Cupcake/'
params.seqtk_path = '$HOME/yan.a/software/seqtk/seqtk'
params.corset_path = '$HOME/yan.a/software/corset-1.09-linux64/corset'
params.oarfish_path = '$HOME/yan.a/software/oarfish-x86_64-unknown-linux-gnu/oarfish'

// Resources
params.cpus = 48

/*
 * Process 1: Sequence reformatting
 */
process REFORMAT_SEQUENCES {
    publishDir "${params.outdir}/01_reformat", mode: 'copy'
    
    input:
    path assembly
    
    output:
    path "transcriptome.reformat.fasta", emit: reformat_fasta
    
    script:
    """
    # Convert U to T if needed (for dRNA assemblies)
    ${params.seqtk_path} seq -a ${assembly} | awk '/^[^>]/{ gsub(/U/,"T") }1' > transcriptome.reformat.fasta
    """
}

/*
 * Process 2: Transrate quality assessment
 */
process TRANSRATE {
    publishDir "${params.outdir}/02_transrate", mode: 'copy'
    
    input:
    path reformat_fasta
    path reference_fasta
    
    output:
    path "transrate/*", emit: results
    
    script:
    """
    source activate ${params.transrate_env}
    
    transrate \\
        --assembly ${reformat_fasta} \\
        --threads ${params.cpus} \\
        --output transrate \\
        --reference ${reference_fasta}
    """
}

/*
 * Process 3: BUSCO completeness assessment
 */
process BUSCO {
    publishDir "${params.outdir}/03_busco", mode: 'copy'
    
    input:
    path reformat_fasta
    path busco_lineage
    
    output:
    path "busco/*", emit: results
    
    script:
    def lineage_name = busco_lineage.name
    """
    source activate ${params.busco_env}
    
    busco \\
        -i ${reformat_fasta} \\
        -l ${params.busco_lineage} \\
        -o busco \\
        -m transcriptome \\
        -c ${params.cpus} \\
        -f \\
        --offline
    
    # Compress results
    tar -czf busco/run_${lineage_name}.tar.gz busco/run_${lineage_name}
    rm -rf busco/run_${lineage_name}
    """
}

/*
 * Process 4: SQANTI3 structural annotation
 */
process SQANTI3 {
    publishDir "${params.outdir}/04_sqanti3", mode: 'copy'
    
    input:
    path reformat_fasta
    path reference_gtf
    path reference_genome
    
    output:
    path "firstpass/", emit: results
    path "unstranded_classification.txt", emit: classification, optional: true
    path "firstpass/firstpass_corrected.fasta", emit: corrected_fasta
    
    script:
    
    """
    set +u
    source activate ${params.sqanti_env}
    
    export PYTHONPATH=\$PYTHONPATH:${params.cupcake_path}
    export PYTHONPATH=\$PYTHONPATH:${params.cupcake_path}/sequence
    set -u

    echo "Starting first pass"
    python ${params.sqanti_path}/sqanti3_qc.py \\
        ${reformat_fasta} \\
        ${reference_gtf} \\
        ${reference_genome} \\
        --CAGE_peak ${params.sqanti_path}/data/ref_TSS_annotation/human.refTSS_v3.1.hg38.bed \\
        --polyA_motif_list ${params.sqanti_path}/data/polyA_motifs/mouse_and_human.polyA_motif.txt \\
        -o firstpass \\
        -d firstpass \\
        --cpus ${params.cpus} \\
        --report skip \\
        --fasta \\
        --force_id_ignore \\
        --skipORF
    
    if [ "${params.stranded}" != "true" ]; then
        echo "Processing unstranded data - running second pass for antisense transcripts"
        
        # Get antisense transcripts
        awk '\$6=="antisense" {print \$1}' firstpass/firstpass_classification.txt | sed 's/^/"/;s/\$/"/' > asid.txt
        
        # Flip strand for antisense transcripts
        grep -F -f asid.txt firstpass/firstpass_corrected.gtf | \\
            awk 'BEGIN{FS=OFS="\\t"} {if(\$7 == "+") \$7="-"; else if(\$7 == "-") \$7="+"; print}' > antisense.gtf
        
        echo "Starting second pass for antisense transcripts"
        python ${params.sqanti_path}/sqanti3_qc.py \\
            antisense.gtf \\
            ${reference_gtf} \\
            ${reference_genome} \\
            --CAGE_peak ${params.sqanti_path}/data/ref_TSS_annotation/human.refTSS_v3.1.hg38.bed \\
            --polyA_motif_list ${params.sqanti_path}/data/polyA_motifs/mouse_and_human.polyA_motif.txt \\
            -o secondpass \\
            -d secondpass \\
            --cpus ${params.cpus} \\
            --report skip \\
            --force_id_ignore \\
            --skipORF
        
        # Merge results
        awk 'BEGIN{FS=OFS="\\t"} {if(\$7 ~ /novelGene/) \$7 = \$7 "_as"; print}' secondpass/secondpass_classification.txt > secondpass/secondpass_classification_tmp.txt
        awk 'BEGIN{FS=OFS="\\t"} NR==FNR{a[\$1]=\$0; next} \$1 in a{\$0=a[\$1]} {print \$0}' secondpass/secondpass_classification_tmp.txt firstpass/firstpass_classification.txt > unstranded_classification.txt
    fi
    
    """
}

/*
 * Process 5: BUSCO on corrected transcriptome
 */
process BUSCO_CORRECTED {
    publishDir "${params.outdir}/05_busco_corrected", mode: 'copy'
    
    input:
    path corrected_fasta
    path busco_lineage
    
    output:
    path "busco_corrected", emit: results
    
    script:
    def lineage_name = busco_lineage.name
    def base_name = corrected_fasta.baseName
    """
    source activate ${params.busco_env}
    
    # Deduplicate headers in corrected fasta
    awk '/^>/ {seen[\$0]++; if(seen[\$0] > 1) {print \$0 "_" seen[\$0]-1} else {print}; next} {print}' \\
        ${corrected_fasta} > ${base_name}_corrected_dedup.fasta
    
    busco \\
        -i ${base_name}_corrected_dedup.fasta \\
        -l ${params.busco_lineage} \\
        -o busco_corrected \\
        -m transcriptome \\
        -c ${params.cpus} \\
        -f
    
    # Compress results
    tar -czf busco_corrected/run_${lineage_name}.tar.gz busco_corrected/run_${lineage_name}
    rm -rf busco_corrected/run_${lineage_name}
    """
}

/*
 * Process 6: Short read quantification (BWA + Salmon) - Merged reads
 */
process SHORT_READ_QUANT {
    publishDir "${params.outdir}/06_salmon_short", mode: 'copy'
    
    input:
    path reformat_fasta
    tuple val(meta), path(reads)
    
    output:
    path "*_quant_*", emit: results
    
    when:
    params.input_type == 'short' || params.input_type == 'hybrid'
    
    script:
    def reads1 = params.single_end ? reads[0] : reads[0]
    def reads2 = params.single_end ? "" : reads[1]
    """
    source activate ${params.isoncorrect_env}
    module load salmon/1.10.2
    module load samtools/1.19.2
    module load bwa/0.7.17
    
    # BWA alignment and quantification
    bwa index ${reformat_fasta}
    
    if [ "${params.single_end}" == "true" ]; then
        bwa mem -t ${params.cpus} ${reformat_fasta} ${reads1} | samtools view -@ ${params.cpus} -Sb > all_unsorted.bam
    else
        bwa mem -t ${params.cpus} ${reformat_fasta} ${reads1} ${reads2} | samtools view -@ ${params.cpus} -Sb > all_unsorted.bam
    fi
    
    samtools view -@ ${params.cpus} -f 2 -F3840 -Sb all_unsorted.bam > all_filtered.bam
    
    salmon quant -p ${params.cpus} -t ${reformat_fasta} -l A -a all_filtered.bam -o all_quant_align
    
    # Salmon direct quantification
    salmon index -t ${reformat_fasta} -i salmon_index --keepDuplicates -p ${params.cpus}
    
    if [ "${params.single_end}" == "true" ]; then
        salmon quant -i salmon_index -l A -r ${reads1} -o all_quant_map -p ${params.cpus} --gcBias --seqBias --posBias --thinningFactor 64
    else
        salmon quant -i salmon_index -l A -1 ${reads1} -2 ${reads2} -o all_quant_map -p ${params.cpus} --gcBias --seqBias --posBias --thinningFactor 64
    fi
    
    rm -r salmon_index
    """
}

/*
 * Process 7: Long read quantification (Minimap2 + Salmon) - Merged reads
 */
process LONG_READ_QUANT {
    publishDir "${params.outdir}/07_salmon_long", mode: 'copy'
    
    input:
    path reformat_fasta
    path long_reads
    
    output:
    path "*_quant_*", emit: results
    
    when:
    params.input_type == 'long' || params.input_type == 'hybrid'
    
    script:
    def minimap_preset = params.long_read_tech == 'hifi' ? 'map-hifi' : 'map-ont'
    """
    source activate ${params.isoncorrect_env}
    module load salmon/1.10.2
    module load samtools/1.19.2
    
    # Minimap2 alignment
    minimap2 -ax ${minimap_preset} -Y -p 1.0 -N 100 -I 1000G -t ${params.cpus} ${reformat_fasta} ${long_reads} | samtools view -@ ${params.cpus} -Sb > all_unsorted.bam
    
    # Quantification with secondary alignments
    salmon quant --ont -p ${params.cpus} -t ${reformat_fasta} -l A --numBootstraps 100 -a all_unsorted.bam -o all_quant_onts
    
    # Quantification with primary alignments only
    samtools view -F 256 -b -@ ${params.cpus} all_unsorted.bam > all_primary.bam
    salmon quant --ont -p ${params.cpus} -t ${reformat_fasta} -l A --numBootstraps 100 -a all_primary.bam -o all_quant_ontp
    
    # Quant with oarfish
    ${params.oarfish_path} -j ${params.cpus} -a all_unsorted.bam -o all_quant_oarfish/all --filter-group no-filters --model-coverage

    """
}

/*
 * Process 8: Corset clustering
 */
process CORSET {
    publishDir "${params.outdir}/08_corset", mode: 'copy'
    
    input:
    path reformat_fasta
    
    output:
    path "corset*", emit: results
    
    script:
    """
    source activate ${params.isoncorrect_env}
    module load samtools/1.19.2
    
    # Create self-alignment BAM
    minimap2 -t ${params.cpus} -a -k15 -w5 -e0 -m100 -r2k -P --dual=yes --no-long-join -I 1000G ${reformat_fasta} ${reformat_fasta} | \\
        samtools sort -@ ${params.cpus} -O BAM -o tx_ovlp.bam
    
    samtools index tx_ovlp.bam
    samtools idxstats tx_ovlp.bam | cut -f1 | grep -v '^*' > refid.txt
    
    # Run Corset clustering
    ${params.corset_path} -f true -m 1 -r true -p corset tx_ovlp.bam
    
    # Handle unmapped transcripts
    comm -13 <(cut -f1 corset-clusters.txt | sort) <(sort refid.txt) | \\
        awk -v OFS='\\t' '{print \$1, "nomap_" NR, 0}' > corset_unmapped_table.txt
    
    # Combine mapped and unmapped results
    cat corset-clusters.txt <(cut -f1,2 corset_unmapped_table.txt) > corset-clusters_mod.txt
    
    """
}

/*
 * Process 9: Individual short read sample quantification - FIXED
 */
process INDIVIDUAL_SHORT_QUANT {
    publishDir "${params.outdir}/09_dge_short", mode: 'copy'
    tag "${sample_name}"
    
    input:
    each path(reformat_fasta)  // Use 'each' to broadcast the fasta to all samples
    tuple val(sample_name), path(reads1), path(reads2) 
    
    output:
    path "${sample_name}_quant_*", emit: results
    
    when:
    params.input_type == 'short' || params.input_type == 'hybrid'
    
    script:
    """
    source activate ${params.isoncorrect_env}
    module load salmon/1.10.2
    module load samtools/1.19.2
    module load bwa/0.7.17
    
    # BWA alignment and quantification
    bwa index ${reformat_fasta}
    
    if [ "${params.single_end}" == "true" ]; then
        bwa mem -t ${params.cpus} ${reformat_fasta} ${reads1} | samtools view -@ ${params.cpus} -Sb > ${sample_name}_unsorted.bam
    else
        bwa mem -t ${params.cpus} ${reformat_fasta} ${reads1} ${reads2} | samtools view -@ ${params.cpus} -Sb > ${sample_name}_unsorted.bam
    fi
    
    samtools view -@ ${params.cpus} -f 2 -F3840 -Sb ${sample_name}_unsorted.bam > ${sample_name}_filtered.bam
    
    salmon quant -p ${params.cpus} -t ${reformat_fasta} -l A --numBootstraps 100 -a ${sample_name}_filtered.bam -o ${sample_name}_quant_align
    rm ${sample_name}_unsorted.bam ${sample_name}_filtered.bam
    
    # Salmon direct quantification
    salmon index -t ${reformat_fasta} -i salmon_index --keepDuplicates -p ${params.cpus}
    
    if [ "${params.single_end}" == "true" ]; then
        salmon quant -i salmon_index -l A --numBootstraps 100 -r ${reads1} -o ${sample_name}_quant_map -p ${params.cpus} --gcBias --seqBias --posBias --thinningFactor 64
    else
        salmon quant -i salmon_index -l A --numBootstraps 100 -1 ${reads1} -2 ${reads2} -o ${sample_name}_quant_map -p ${params.cpus} --gcBias --seqBias --posBias --thinningFactor 64
    fi
    
    rm -r salmon_index
    """
}

/*
 * Process 10: Individual long read sample quantification - FIXED
 */
process INDIVIDUAL_LONG_QUANT {
    publishDir "${params.outdir}/10_dge_long", mode: 'copy'
    tag "${sample_name}"
    
    input:
    each path(reformat_fasta)  // Use 'each' to broadcast the fasta to all samples
    tuple val(sample_name), path(long_reads)
    
    output:
    path "${sample_name}_quant_*", emit: results
    
    when:
    params.input_type == 'long' || params.input_type == 'hybrid'
    
    script:
    def minimap_preset = params.long_read_tech == 'hifi' ? 'map-hifi' : 'map-ont'
    """
    source activate ${params.isoncorrect_env}
    module load salmon/1.10.2
    module load samtools/1.19.2
    
    # Minimap2 alignment
    minimap2 -ax ${minimap_preset} -Y -p 1.0 -N 100 -I 1000G -t ${params.cpus} ${reformat_fasta} ${long_reads} | samtools view -@ ${params.cpus} -Sb > ${sample_name}_unsorted.bam
    
    samtools sort -@ ${params.cpus} ${sample_name}_unsorted.bam > ${sample_name}_sorted.bam
    samtools index -@ ${params.cpus} ${sample_name}_sorted.bam
    samtools flagstat ${sample_name}_sorted.bam > ${sample_name}.flagstat
    
    # Quantification with secondary alignments
    salmon quant --ont -p ${params.cpus} -t ${reformat_fasta} -l A --numBootstraps 100 -a ${sample_name}_unsorted.bam -o ${sample_name}_quant_onts
    
    # Quantification with primary alignments only
    samtools view -F 256 -b -@ ${params.cpus} ${sample_name}_unsorted.bam > ${sample_name}.primary.bam
    salmon quant --ont -p ${params.cpus} -t ${reformat_fasta} -l A --numBootstraps 100 -a ${sample_name}.primary.bam -o ${sample_name}_quant_ontp
    
    # Quant with oarfish
    ${params.oarfish_path} -j ${params.cpus} -a ${sample_name}_unsorted.bam -o ${sample_name}_quant_oarfish/${sample_name} --filter-group no-filters --model-coverage

    """
}

/*
 * Main workflow
 */
workflow {
    // Input validation
    if (!params.assembly) {
        error "Please provide --assembly parameter"
    }
    if (!params.reference_fasta) {
        error "Please provide --reference_fasta parameter"
    }
    if (!params.reference_gtf) {
        error "Please provide --reference_gtf parameter"
    }
    if (!params.reference_genome) {
        error "Please provide --reference_genome parameter"
    }
    if (!params.busco_lineage) {
        error "Please provide --busco_lineage parameter"
    }
    
    // Create input channels
    assembly_ch = Channel.fromPath(params.assembly, checkIfExists: true)
    reference_fasta_ch = Channel.fromPath(params.reference_fasta, checkIfExists: true)
    reference_gtf_ch = Channel.fromPath(params.reference_gtf, checkIfExists: true)
    reference_genome_ch = Channel.fromPath(params.reference_genome, checkIfExists: true)
    busco_lineage_ch = Channel.fromPath(params.busco_lineage, checkIfExists: true)
    
    // Step 1: Reformat sequences
    REFORMAT_SEQUENCES(assembly_ch)

    // Step 2: Transrate quality assessment
    TRANSRATE(REFORMAT_SEQUENCES.out.reformat_fasta, reference_fasta_ch)
    
    // Step 3: BUSCO completeness assessment
    BUSCO(REFORMAT_SEQUENCES.out.reformat_fasta, busco_lineage_ch)
    
    // Step 4: SQANTI3 structural annotation
    SQANTI3(REFORMAT_SEQUENCES.out.reformat_fasta, reference_gtf_ch, reference_genome_ch)
    
    // Step 5: BUSCO on corrected transcriptome
    BUSCO_CORRECTED(SQANTI3.out.corrected_fasta, busco_lineage_ch)
    
    // Step 6: Corset clustering
    CORSET(REFORMAT_SEQUENCES.out.reformat_fasta)

    // Step 7: Merged quantification (if single files provided)
    if (params.short_reads && (params.input_type == 'short' || params.input_type == 'hybrid')) {
        if (params.single_end) {
            // Single-end reads
            short_reads_ch = Channel.fromPath(params.short_reads, checkIfExists: true)
                .map { file -> [ [id: 'merged'], [file] ] }
        } else {
            // Paired-end reads
            if (!params.short_reads2) {
                error "For paired-end reads, please provide both --short_reads and --short_reads2"
            }
            short_reads_ch = Channel.fromPath([params.short_reads, params.short_reads2], checkIfExists: true)
                .collect()
                .map { files -> [ [id: 'merged'], files ] }
        }
        SHORT_READ_QUANT(REFORMAT_SEQUENCES.out.reformat_fasta, short_reads_ch)
    }
    
    if (params.long_reads && (params.input_type == 'long' || params.input_type == 'hybrid')) {
        long_reads_ch = Channel.fromPath(params.long_reads, checkIfExists: true)
        LONG_READ_QUANT(REFORMAT_SEQUENCES.out.reformat_fasta, long_reads_ch)
    }
    
    // Step 8: Individual sample quantification (from list files) - FIXED VERSION
    if (params.short_list && (params.input_type == 'short' || params.input_type == 'hybrid')) {
        // Create channel from short read list file
        short_samples_ch = Channel
            .fromPath(params.short_list, checkIfExists: true)
            .splitText()
            .map { line ->
                def read1_path = line.trim()
                if (!read1_path) return null  // Skip empty lines
                
                // Remove suffix from the full path
                def path_without_suffix = read1_path
                if (read1_path.endsWith(params.short_suffix1)) {
                    path_without_suffix = read1_path.replace(params.short_suffix1, "")
                }
                
                // Get basename of the remaining path
                def sample_name = "SR_" + file(path_without_suffix).name

                if (params.single_end) {
                    return tuple(sample_name, file(read1_path), file("NO_FILE"))
                } else {
                    def read2_path = read1_path.replaceAll(params.short_suffix1, params.short_suffix2)
                    return tuple(sample_name, file(read1_path), file(read2_path))
                }
            }
            .filter { it != null }  // Remove null entries from empty lines

        INDIVIDUAL_SHORT_QUANT(REFORMAT_SEQUENCES.out.reformat_fasta, short_samples_ch)
    }

    if (params.long_list && (params.input_type == 'long' || params.input_type == 'hybrid')) {
        // Create channel from long read list file
        long_samples_ch = Channel
            .fromPath(params.long_list, checkIfExists: true)
            .splitText()
            .map { line ->
                def long_read_path = line.trim()
                if (!long_read_path) return null  // Skip empty lines
                
                // Remove suffix from the full path
                def path_without_suffix = long_read_path
                if (long_read_path.endsWith(params.long_suffix)) {
                    path_without_suffix = long_read_path.replace(params.long_suffix, "")
                }

                // Get basename of the remaining path
                def sample_name = "LR_" + file(path_without_suffix).name
                
                return tuple(sample_name, file(long_read_path))
            }
            .filter { it != null }  // Remove null entries from empty lines

        INDIVIDUAL_LONG_QUANT(REFORMAT_SEQUENCES.out.reformat_fasta, long_samples_ch)
    }
}