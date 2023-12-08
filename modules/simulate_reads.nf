process badread_simulate {

    tag { assembly_id + ' / ' + fold_coverage + 'x' + ' / ' + 'replicate=' + replicate }

    publishDir "${params.outdir}/${output_subdir}", pattern: "${assembly_id}-${md5_fragment}*_read_simulation_parameters.csv", mode: 'copy'

    input:
    tuple val(assembly_id), path(assembly), val(fold_coverage), val(replicate)

    output:
    tuple val(assembly_id), val(md5_fragment), path("${assembly_id}-${md5_fragment}*_L.fastq"), emit: reads
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
	> ${assembly_id}-${md5_fragment}_L.fastq

    echo 'sample_id,replicate,random_seed,fold_coverage,mean_read_length,stdev_read_length,junk_reads_percent,random_reads_percent,chimeras_percent' > ${assembly_id}-${md5_fragment}_read_simulation_parameters.csv
    echo '${assembly_id}-${md5_fragment},${replicate},${seed},${fold_coverage},${mean_read_length},${stdev_read_length},${junk_reads},${random_reads},${chimeras}' >> ${assembly_id}-${md5_fragment}_read_simulation_parameters.csv
    """
}

process simulate_contaminant_reads {

    tag { contaminant_id + ' / ' + proportion }

    input:
    tuple val(contaminant_id), path(assembly), val(proportion)

    output:
    tuple val(contaminant_id), path("${contaminant_id}*_L.fastq"), val(proportion)

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
	> ${assembly_id}-${md5_fragment}_L.fastq
    """
}

process downsample_simulated_reads {

    tag { assembly_id + '-' + md5_fragment + ' / ' + proportion_uncontaminated }

    input:
    tuple val(assembly_id), val(md5_fragment), path(assembly_reads), val(proportion_contaminants)

    output:
    tuple val(assembly_id), val(md5_fragment), path("${assembly_id}-${md5_fragment}*_L.fastq")

    script:
    proportion_uncontaminated = 1.0 - proportion_contaminants
    seed = Math.round(Math.random() * 1000000)
    """
    seqkit sample \
	-s ${seed} \
	-p ${proportion_uncontaminated} \
	${assembly_reads} \
	-o ${assembly_id}-${md5_fragment}_sample_L.fastq
    """
}

process downsample_contaminant_reads {

    tag { assembly_id + '-' + md5_fragment + ' / ' + contaminant_id + ' / ' + contaminant_proportion }

    publishDir "${params.outdir}/${output_subdir}/contaminants", pattern: "${contaminant_id}_contaminant_L*.fastq.gz", mode: 'copy'
    publishDir "${params.outdir}/${output_subdir}/contaminants", pattern: "${assembly_id}-${md5_fragment}-${contaminant_id}_num_contaminant_reads.csv", mode: 'copy'
    
    input:
    tuple val(contaminant_id), path(contaminant_reads), val(contaminant_proportion), val(assembly_id), val(md5_fragment), path(assembly_reads)

    output:
    tuple val(assembly_id), val(md5_fragment), val(contaminant_id), path("${contaminant_id}_contaminant_L.fastq"), emit: uncompressed_reads
    tuple val(assembly_id), val(md5_fragment), val(contaminant_id), path("${contaminant_id}_contaminant_L.fastq.gz"), emit: compressed_reads
    tuple val(assembly_id), val(md5_fragment), val(contaminant_id), path("${assembly_id}-${md5_fragment}-${contaminant_id}_num_contaminant_reads.csv"), emit: num_reads_csv

    script:
    output_subdir = params.flat ? '' : assembly_id + '-' + md5_fragment
    seed = Math.round(Math.random() * 1000000)
    """
    seqkit stats -T ${assembly_reads} | tail -n 1 | cut -f 4 | tr -d ',' > num_simulated_reads
    echo 'sample_id,contaminant_id,num_simulated_reads,num_contaminant_reads,target_contaminant_proportion' > ${assembly_id}-${md5_fragment}-${contaminant_id}_num_contaminant_read_pairs.csv
    python -c "import sys; print(int(round(int(sys.stdin.read().strip()) * ${contaminant_proportion})))" < num_simulated_read_pairs > num_contaminant_read_pairs
    paste -d ',' <(echo "${assembly_id}-${md5_fragment}") <(echo "${contaminant_id}") num_simulated_reads num_contaminant_reads <(echo "${contaminant_proportion}") >> ${assembly_id}-${md5_fragment}-${contaminant_id}_num_contaminant_reads.csv
    seqkit sample -s ${seed} -n \$(cat num_contaminant_reads) ${contaminant_reads} > ${contaminant_id}_contaminant_L.fastq
    gzip --keep ${contaminant_id}_contaminant_L*.fastq
    """
}

process introduce_contaminants {
    
    tag { assembly_id + '-' + md5_fragment }
    
    input:
    tuple val(assembly_id), val(md5_fragment), path(assembly_reads), val(contaminant_ids), path(contaminant_reads)
    
    output:
    tuple val(assembly_id), val(md5_fragment), path("${assembly_id}-${md5_fragment}_L.fastq"), emit: reads
    
    script:
    seed = Math.round(Math.random() * 1000000)
    """
    mv ${assembly_reads} uncontaminated_L.fastq
    cat uncontaminated_L.fastq ${contaminant_reads} > ${assembly_id}-${md5_fragment}_unshuffled_L.fastq
    cat ${assembly_id}-${md5_fragment}_unshuffled_L.fastq \
	| paste - - - - | shuf | awk -F'\\t' '{OFS="\\n"; print \$1,\$2,\$3,\$4 > "${assembly_id}-${md5_fragment}_L.fastq";}'
    """
}

process fastp {

    tag { assembly_id + '-' + md5_fragment }
    
    publishDir "${params.outdir}/${output_subdir}", pattern: "${assembly_id}-${md5_fragment}_L.fastq.gz", mode: 'copy'
    publishDir "${params.outdir}/${output_subdir}", pattern: "${assembly_id}-${md5_fragment}_fastp.json", mode: 'copy'
    
    input:
    tuple val(assembly_id), val(md5_fragment), path(reads)
    
    output:
    tuple val(assembly_id), val(md5_fragment), path("${assembly_id}-${md5_fragment}_fastp.json"), emit: json
    tuple val(assembly_id), val(md5_fragment), path("${assembly_id}-${md5_fragment}_fastp.csv"), emit: csv
    tuple val(assembly_id), val(md5_fragment), path("${assembly_id}-${md5_fragment}_L.fastq.gz"), emit: untrimmed_reads
    
    script:
    output_subdir = params.flat ? '' : assembly_id + '-' + md5_fragment
    """
    fastp -i ${reads} -o ${assembly_id}-${md5_fragment}_trimmed_L.fastq.gz
    mv fastp.json ${assembly_id}-${md5_fragment}_fastp.json
    fastp_json_to_csv.py -s ${assembly_id}-${md5_fragment} ${assembly_id}-${md5_fragment}_fastp.json > ${assembly_id}-${md5_fragment}_fastp.csv
    cp ${reads} untrimmed_L.fastq
    gzip -c untrimmed_L.fastq > ${assembly_id}-${md5_fragment}_L.fastq.gz
    """
}

process minimap2_align {

    tag { assembly_id + '-' + md5_fragment + ' / ' + assembly_id }
    
    publishDir "${params.outdir}/${output_subdir}", pattern: "${assembly_id}-${md5_fragment}.{bam,bam.bai}", mode: 'copy'
    
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
    
    script:
    output_subdir = params.flat ? '' : assembly_id + '-' + md5_fragment
    """
    qualimap bamqc -bam ${alignment} --outdir ${assembly_id}-${md5_fragment}_bamqc
    """
}

process qualimap_bamqc_genome_results_to_csv {

    tag { assembly_id + '-' + md5_fragment }
    
    publishDir "${params.outdir}/${output_subdir}", mode: 'copy', pattern: "${assembly_id}-${md5_fragment}_qualimap_bamqc_genome_results.csv"
    
    executor 'local'
    
    input:
    tuple val(assembly_id), val(md5_fragment), path(qualimap_bamqc_genome_results)
    
    output:
    tuple val(assembly_id), val(md5_fragment), path("${assembly_id}-${md5_fragment}_qualimap_bamqc_genome_results.csv")
    
    script:
    output_subdir = params.flat ? '' : assembly_id + '-' + md5_fragment
    """
    qualimap_bamqc_genome_results_to_csv.py -s ${assembly_id}-${md5_fragment} ${qualimap_bamqc_genome_results} > ${assembly_id}-${md5_fragment}_qualimap_bamqc_genome_results.csv
    """
}

process samtools_stats {

    tag { assembly_id + '-' + md5_fragment }

    publishDir "${params.outdir}/${output_subdir}", mode: 'copy', pattern: "${assembly_id}-${md5_fragment}_samtools_stats_summary.txt"

    input:
    tuple val(assembly_id), val(md5_fragment), file(alignment), file(alignment_index)
    
    output:
    tuple val(assembly_id), val(md5_fragment), path("${assembly_id}-${md5_fragment}_samtools_stats_summary.txt"), emit: summary
  
    script:
    output_subdir = params.flat ? '' : assembly_id + '-' + md5_fragment
    """
    samtools stats ${alignment} > ${assembly_id}-${md5_fragment}_samtools_stats.txt
    grep ^SN ${assembly_id}-${md5_fragment}_samtools_stats.txt | cut -f 2- > ${assembly_id}-${md5_fragment}_samtools_stats_summary.txt
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
