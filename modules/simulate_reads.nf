process badread_simulate {

    tag { assembly_id + ' / ' + fold_coverage + 'x' + ' / ' + 'replicate=' + replicate }

    publishDir "${params.outdir}/${output_subdir}", pattern: "${assembly_id}-${md5_fragment}*_read_simulation_parameters.csv", mode: 'copy'

    input:
    tuple val(assembly_id), path(assembly), val(fold_coverage), val(replicate)

    output:
    tuple val(assembly_id), val(md5_fragment), path("${assembly_id}-${md5_fragment}*_RL.fastq"), emit: reads
    tuple val(assembly_id), val(md5_fragment), path("${assembly_id}-${md5_fragment}_read_simulation_parameters.csv"), emit: metrics

    script:
    mean_read_length = params.mean_read_length
    stdev_read_length = params.stdev_read_length
    junk_reads = params.percent_junk_reads
    random_reads = params.percent_random_reads
    chimeras = params.percent_chimeras
    seed = Math.round(Math.random() * 1000000)
    md5_input = assembly_id + fold_coverage.toString() + mean_read_length.toString() + stdev_read_length.toString() + junk_reads.toString() + random_reads.toString() + chimeras.toString() + seed.toString()
    md5_fragment = md5_input.md5()[0..3]
    output_subdir = params.flat ? '' : assembly_id + '-' + md5_fragment
    """
    
    badread simulate \
	--seed ${seed} \
	--reference ${assembly} \
	--length ${mean_read_length},${stdev_read_length} \
	--quantity ${fold_coverage}x \
	--junk_reads ${junk_reads} \
	--random_reads ${random_reads} \
	--chimeras ${chimeras} \
	> ${assembly_id}-${md5_fragment}_RL.fastq

    echo 'sample_id,replicate,random_seed,fold_coverage,mean_read_length,stdev_read_length,junk_reads_percent,random_reads_percent,chimeras_percent' > ${assembly_id}-${md5_fragment}_read_simulation_parameters.csv
    echo '${assembly_id}-${md5_fragment},${replicate},${seed},${fold_coverage},${mean_read_length},${stdev_read_length},${junk_reads},${random_reads},${chimeras}' >> ${assembly_id}-${md5_fragment}_read_simulation_parameters.csv
    """
}

process simulate_contaminant_reads {

    tag { contaminant_id + ' / ' + proportion }

    input:
    tuple val(contaminant_id), path(assembly), val(proportion)

    output:
    tuple val(contaminant_id), path("${contaminant_id}*_RL.fastq"), val(proportion)

    script:
    mean_read_length = params.mean_read_length
    stdev_read_length = params.stdev_read_length
    junk_reads = params.percent_junk_reads
    random_reads = params.percent_random_reads
    chimeras = params.percent_chimeras
    seed = Math.round(Math.random() * 1000000)
    """
    badread simulate \
	--seed ${seed} \
	--reference ${assembly} \
	--length ${mean_read_length},${stdev_read_length} \
	--quantity 30x \
	--junk_reads ${junk_reads} \
	--random_reads ${random_reads} \
	--chimeras ${chimeras} \
	> ${assembly_id}-${md5_fragment}_RL.fastq
    """
}

process downsample_simulated_reads {

    tag { assembly_id + '-' + md5_fragment + ' / ' + proportion_uncontaminated }

    input:
    tuple val(assembly_id), val(md5_fragment), path(assembly_reads), val(proportion_contaminants)

    output:
    tuple val(assembly_id), val(md5_fragment), path("${assembly_id}-${md5_fragment}*_RL.fastq")

    script:
    proportion_uncontaminated = 1.0 - proportion_contaminants
    seed = Math.round(Math.random() * 1000000)
    """
    seqkit sample \
	-s ${seed} \
	-p ${proportion_uncontaminated} \
	${assembly_reads} \
	-o ${assembly_id}-${md5_fragment}_sample_RL.fastq
    """
}

process downsample_contaminant_reads {

    tag { assembly_id + '-' + md5_fragment + ' / ' + contaminant_id + ' / ' + contaminant_proportion }

    publishDir "${params.outdir}/${output_subdir}/contaminants", pattern: "${contaminant_id}_contaminant_RL*.fastq.gz", mode: 'copy'
    publishDir "${params.outdir}/${output_subdir}/contaminants", pattern: "${assembly_id}-${md5_fragment}-${contaminant_id}_num_contaminant_reads.csv", mode: 'copy'
    
    input:
    tuple val(contaminant_id), path(contaminant_reads), val(contaminant_proportion), val(assembly_id), val(md5_fragment), path(assembly_reads)

    output:
    tuple val(assembly_id), val(md5_fragment), val(contaminant_id), path("${contaminant_id}_contaminant_RL.fastq"), emit: uncompressed_reads
    tuple val(assembly_id), val(md5_fragment), val(contaminant_id), path("${contaminant_id}_contaminant_RL.fastq.gz"), emit: compressed_reads
    tuple val(assembly_id), val(md5_fragment), val(contaminant_id), path("${assembly_id}-${md5_fragment}-${contaminant_id}_num_contaminant_reads.csv"), emit: num_reads_csv

    script:
    output_subdir = params.flat ? '' : assembly_id + '-' + md5_fragment
    seed = Math.round(Math.random() * 1000000)
    """
    seqkit stats -T ${assembly_reads} | tail -n 1 | cut -f 4 | tr -d ',' > num_simulated_reads

    echo 'sample_id,contaminant_id,num_simulated_reads,num_contaminant_reads,target_contaminant_proportion' > ${assembly_id}-${md5_fragment}-${contaminant_id}_num_contaminant_read_pairs.csv

    python -c "import sys; print(int(round(int(sys.stdin.read().strip()) * ${contaminant_proportion})))" < num_simulated_read_pairs > num_contaminant_read_pairs

    paste -d ',' <(echo "${assembly_id}-${md5_fragment}") <(echo "${contaminant_id}") num_simulated_reads num_contaminant_reads <(echo "${contaminant_proportion}") >> ${assembly_id}-${md5_fragment}-${contaminant_id}_num_contaminant_reads.csv

    seqkit sample -s ${seed} -n \$(cat num_contaminant_reads) ${contaminant_reads} > ${contaminant_id}_contaminant_RL.fastq

    gzip --keep ${contaminant_id}_contaminant_RL*.fastq
    """
}

process introduce_contaminants {
    
    tag { assembly_id + '-' + md5_fragment }
    
    input:
    tuple val(assembly_id), val(md5_fragment), path(assembly_reads), val(contaminant_ids), path(contaminant_reads)
    
    output:
    tuple val(assembly_id), val(md5_fragment), path("${assembly_id}-${md5_fragment}_RL.fastq"), emit: reads
    
    script:
    seed = Math.round(Math.random() * 1000000)
    """
    mv ${assembly_reads} uncontaminated_RL.fastq
    cat uncontaminated_RL.fastq ${contaminant_reads} > ${assembly_id}-${md5_fragment}_unshuffled_RL.fastq
    cat ${assembly_id}-${md5_fragment}_unshuffled_RL.fastq \
	| paste - - - - | shuf | awk -F'\\t' '{OFS="\\n"; print \$1,\$2,\$3,\$4 > "${assembly_id}-${md5_fragment}_RL.fastq";}'
    """
}

process fastp {

    tag { assembly_id + '-' + md5_fragment }
    
    publishDir "${params.outdir}/${output_subdir}", pattern: "${assembly_id}-${md5_fragment}_RL.fastq.gz", mode: 'copy'
    publishDir "${params.outdir}/${output_subdir}", pattern: "${assembly_id}-${md5_fragment}_fastp.json", mode: 'copy'
    
    input:
    tuple val(assembly_id), val(md5_fragment), path(reads)
    
    output:
    tuple val(assembly_id), val(md5_fragment), path("${assembly_id}-${md5_fragment}_fastp.json"), emit: json
    tuple val(assembly_id), val(md5_fragment), path("${assembly_id}-${md5_fragment}_fastp.csv"), emit: fastp_csv
    tuple val(assembly_id), val(md5_fragment), path("${assembly_id}-${md5_fragment}_RL.fastq.gz"), emit: untrimmed_reads
    
    script:
    output_subdir = params.flat ? '' : assembly_id + '-' + md5_fragment
    """
    fastp -i ${reads} -o ${assembly_id}-${md5_fragment}_trimmed_RL.fastq.gz

    mv fastp.json ${assembly_id}-${md5_fragment}_fastp.json

    fastp_json_to_csv.py -s ${assembly_id}-${md5_fragment} ${assembly_id}-${md5_fragment}_fastp.json > ${assembly_id}-${md5_fragment}_fastp.csv

    cp ${reads} untrimmed_RL.fastq

    gzip -c untrimmed_RL.fastq > ${assembly_id}-${md5_fragment}_RL.fastq.gz
    """
}

process minimap2_align {

    tag { assembly_id + '-' + md5_fragment + ' / ' + assembly_id }
    
    publishDir "${params.outdir}/${output_subdir}", pattern: "${assembly_id}-${md5_fragment}.{bam,bam.bai}", mode: 'copy', enabled: params.keep_bams
    
    input:
    tuple val(assembly_id), val(md5_fragment), path(reads), path(ref)
    
    output:
    tuple val(assembly_id), val(md5_fragment), path("${assembly_id}-${md5_fragment}.bam"), path("${assembly_id}-${md5_fragment}.bam.bai")
    
    script:
    output_subdir = params.flat ? '' : assembly_id + '-' + md5_fragment
    """
    minimap2 -ax map-ont \
	-t ${task.cpus} \
	${ref} \
	${reads} | \
	samtools sort -o ${assembly_id}-${md5_fragment}.bam -
	samtools index ${assembly_id}-${md5_fragment}.bam
    """
}

process qualimap_bamqc {

    tag { assembly_id + '-' + md5_fragment }
    
    publishDir "${params.outdir}/${output_subdir}", mode: 'copy', pattern: "${assembly_id}-${md5_fragment}_bamqc"
    
    input:
    tuple val(assembly_id), val(md5_fragment), file(alignment), file(alignment_index)
    
    output:
    tuple val(assembly_id), val(md5_fragment), path("${assembly_id}-${md5_fragment}_bamqc/genome_results.txt"), emit: genome_results
    tuple val(assembly_id), val(md5_fragment), path("${assembly_id}-${md5_fragment}_bamqc"), emit: bamqc_dir
    tuple val(assembly_id), val(md5_fragment), path("${assembly_id}-${md5_fragment}_qualimap_alignment_qc.csv"), emit: alignment_qc

    script:
    output_subdir = params.flat ? '' : assembly_id + '-' + md5_fragment
    """
    qualimap bamqc \
	-bam ${alignment} \
	--outdir ${assembly_id}-${md5_fragment}_bamqc

    qualimap_bamqc_genome_results_to_csv.py \
	-s ${assembly_id}-${md5_fragment} \
	${assembly_id}-${md5_fragment}_bamqc/genome_results.txt \
	> ${assembly_id}-${md5_fragment}_qualimap_alignment_qc.csv
    """
}


process samtools_stats {

    tag { assembly_id + '-' + md5_fragment }

    publishDir "${params.outdir}/${output_subdir}", mode: 'copy', pattern: "${assembly_id}-${md5_fragment}_samtools_stats*"

    input:
    tuple val(assembly_id), val(md5_fragment), path(alignment), path(alignment_index)

    output:
    tuple val(assembly_id), val(md5_fragment), path("${assembly_id}-${md5_fragment}_samtools_stats.txt"), emit: stats
    tuple val(assembly_id), val(md5_fragment), path("${assembly_id}-${md5_fragment}_samtools_stats_summary.txt"), emit: stats_summary
    tuple val(assembly_id), val(md5_fragment), path("${assembly_id}-${md5_fragment}_samtools_stats_summary.csv"), emit: stats_summary_csv
    tuple val(assembly_id), val(md5_fragment), path("${assembly_id}-${md5_fragment}_samtools_stats_insert_sizes.tsv"), emit: insert_sizes
    tuple val(assembly_id), val(md5_fragment), path("${assembly_id}-${md5_fragment}_samtools_stats_coverage_distribution.tsv"), emit: coverage_distribution

    script:
    output_subdir = params.flat ? '' : assembly_id + '-' + md5_fragment
    """
    samtools stats \
	--threads ${task.cpus} \
	${alignment[0]} > ${assembly_id}-${md5_fragment}_samtools_stats.txt

    grep '^SN' ${assembly_id}-${md5_fragment}_samtools_stats.txt | cut -f 2-  > ${assembly_id}-${md5_fragment}_samtools_stats_summary.txt

    parse_samtools_stats_summary.py -i ${assembly_id}-${md5_fragment}_samtools_stats_summary.txt -s ${assembly_id} > ${assembly_id}-${md5_fragment}_samtools_stats_summary.csv

    echo "insert_size,pairs_total,inward_oriented_pairs,outward_oriented_pairs,other_pairs" | tr ',' '\t' > ${assembly_id}-${md5_fragment}_samtools_stats_insert_sizes.tsv
    grep '^IS' ${assembly_id}-${md5_fragment}_samtools_stats.txt | cut -f 2-  >> ${assembly_id}-${md5_fragment}_samtools_stats_insert_sizes.tsv

    echo "coverage,depth" | tr ',' '\t' > ${assembly_id}-${md5_fragment}_samtools_stats_coverage_distribution.tsv
    grep '^COV' ${assembly_id}-${md5_fragment}_samtools_stats.txt | cut -f 2- >> ${assembly_id}-${md5_fragment}_samtools_stats_coverage_distribution.tsv	
    """
}


process combine_alignment_qc {

    tag { assembly_id + '-' + md5_fragment }

    publishDir "${params.outdir}/${output_subdir}", mode: 'copy', pattern: "${assembly_id}-${md5_fragment}_combined_alignment_qc.csv"

    input:
    tuple val(assembly_id), val(md5_fragment), path(qualimap_genome_results_csv), path(samtools_stats_summary_csv)

    output:
    tuple val(assembly_id), val(md5_fragment), path("${assembly_id}-${md5_fragment}_combined_alignment_qc.csv")

    script:
    output_subdir = params.flat ? '' : assembly_id + '-' + md5_fragment
    """
    combine_alignment_qc.py \
	--sample-id ${assembly_id}-${md5_fragment} \
	--read-type "long" \
	--qualimap-bamqc-genome-results ${qualimap_genome_results_csv} \
	--samtools-stats-summary ${samtools_stats_summary_csv} \
	> ${assembly_id}-${md5_fragment}_combined_alignment_qc.csv
    """
}


process mosdepth {
    
    tag { assembly_id + '-' + md5_fragment }
    
    publishDir "${params.outdir}/${output_subdir}", mode: 'copy', pattern: "${assembly_id}-${md5_fragment}_samtools_stats_summary.txt"
    
    input:
    tuple val(assembly_id), val(md5_fragment), file(alignment), file(alignment_index), val(depth_by)

    output:
    tuple val(assembly_id), val(md5_fragment), path("${assembly_id}-${md5_fragment}_by_${depth_by}.mosdepth.summary.txt"), emit: summary
    tuple val(assembly_id), val(md5_fragment), path("${assembly_id}-${md5_fragment}_by_${depth_by}.mosdepth.summary.txt"), emit: regions
    
    script:
    output_subdir = params.flat ? '' : assembly_id + '-' + md5_fragment
    """
    mosdepth -t ${task.cpus} --fast-mode --by ${depth_by} --no-per-base ${assembly_id}-${md5_fragment}_by_${depth_by} ${alignment}
    gunzip ${assembly_id}-${md5_fragment}_by_${depth_by}.regions.bed.gz
    mv ${assembly_id}-${md5_fragment}_by_${depth_by}.regions.bed ${assembly_id}-${md5_fragment}_by_${depth_by}.mosdepth.regions.bed
    """
}
