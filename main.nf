#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { badread_simulate as simulate_reads }   from './modules/simulate_reads.nf'
include { simulate_contaminant_reads }           from './modules/simulate_reads.nf'
include { downsample_simulated_reads }           from './modules/simulate_reads.nf'
include { downsample_contaminant_reads }         from './modules/simulate_reads.nf'
include { introduce_contaminants }               from './modules/simulate_reads.nf'
include { fastp }                                from './modules/simulate_reads.nf'
include { minimap2_align }                       from './modules/simulate_reads.nf'
include { qualimap_bamqc }                       from './modules/simulate_reads.nf'
include { samtools_stats }                       from './modules/simulate_reads.nf'
include { combine_alignment_qc }                 from './modules/simulate_reads.nf'


workflow {
    ch_assemblies = Channel.fromPath( params.assembly_search_path ).map{ it -> [it.baseName, it] }.unique{ it -> it[0] }

    if (params.depths_file == 'NO_FILE'){
	ch_depths = Channel.fromList([params.depth])
    } else {
	ch_depths = Channel.fromPath(params.depths_file).splitCsv(header: false)
    }

    ch_replicates = Channel.fromList([1..params.replicates][0])

    if (params.contaminants != 'NO_FILE') {
	ch_contaminants = Channel.fromPath(params.contaminants).splitCsv(header: true).map{ it -> [it['ID'], it['ASSEMBLY'], it['PROPORTION']] }
	simulate_contaminant_reads(ch_contaminants)
	simulate_reads(ch_assemblies.combine(ch_depths).combine(ch_replicates))
	ch_proportion_uncontaminated = ch_contaminants.map{ it -> Float.parseFloat(it[2]) }.reduce{ x, y -> (x + y).round(6) }
	downsample_simulated_reads(simulate_reads.out.reads.combine(ch_proportion_uncontaminated))
	downsample_contaminant_reads(simulate_contaminant_reads.out.combine(simulate_reads.out.reads))
	introduce_contaminants(downsample_simulated_reads.out.join(downsample_contaminant_reads.out.uncompressed_reads.groupTuple(by: [0, 1]), by: [0, 1]))
	ch_reads = introduce_contaminants.out.reads
    } else {
	simulate_reads(ch_assemblies.combine(ch_depths).combine(ch_replicates))
	ch_reads = simulate_reads.out.reads
    }

    main:
    fastp(ch_reads)
    minimap2_align(ch_assemblies.cross(ch_reads).map{ it -> [it[1][0], it[1][1], it[1][2], it[0][1]] })
    samtools_stats(minimap2_align.out)
    qualimap_bamqc(minimap2_align.out)
    combine_alignment_qc(qualimap_bamqc.out.alignment_qc.join(samtools_stats.out.stats_summary_csv, by: [0, 1]))

    if (params.collect_outputs) {
	fastp.out.fastp_csv.map{ it -> it[2] }.collectFile(keepHeader: true, sort: { it.text }, name: "${params.collected_outputs_prefix}_fastp.csv", storeDir: "${params.outdir}")
	qualimap_bamqc.out.alignment_qc.map{ it -> it[2] }.collectFile(keepHeader: true, sort: { it.text }, name: "${params.collected_outputs_prefix}_qualimap_alignment_qc.csv", storeDir: "${params.outdir}")
	samtools_stats.out.stats_summary_csv.map{ it -> it[2] }.collectFile(keepHeader: true, sort: { it.text }, name: "${params.collected_outputs_prefix}_samtools_stats_summary.csv", storeDir: "${params.outdir}")
	combine_alignment_qc.out.map{ it -> it[2] }.collectFile(keepHeader: true, sort: { it.text }, name: "${params.collected_outputs_prefix}_combined_alignment_qc.csv", storeDir: "${params.outdir}")
    }

}
