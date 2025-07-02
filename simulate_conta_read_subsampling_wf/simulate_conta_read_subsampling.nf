nextflow.enable.dsl=2

// Mandatory parameters
if (!params.sequencing_data_dir)  error "Missing required parameter: --sequencing_data_dir"
if (!params.conta_desc_table)  error "Missing required parameter: --conta_desc_table"
if (!params.output_dir)  error "Missing required parameter: --output_dir"


// Optional/default parameters
params.conda_env = "conda_envs/simulate_conta_read_subsampling_env.yml"
params.fastq_extension   = "fastq.gz"            // e.g. "fastq.gz" or "fq.gz"
params.pigz_threads      = 4                     // number of threads for pigz
params.seed_base         = 42                    // base seed to ensure reproducibility with variation per sample

Channel
    .fromPath(params.conta_desc_table)
    .splitCsv(header:true, sep:'\t')
    .map { row ->
        def sequencing_data_dir = file(params.sequencing_data_dir)
        def hash = row.contaminated_sample_name.hashCode().abs() % 10000
        def seed = params.seed_base + hash
        tuple(
            row.contaminated_sample_name,
            sequencing_data_dir/row.source_sample_name,
            sequencing_data_dir/row.sink_sample_name,
            row.num_reads_from_source.toInteger(),
            row.num_reads_from_sink.toInteger(),
            row.library_layout.toLowerCase(),  // expects 'single' or 'paired'
            seed
        )
    }
    .set { contamination_tasks }

process simulate_contamination {
    conda "${params.conda_env}"

    tag "${contaminated_sample_name}"

    input:
    tuple val(contaminated_sample_name), path(source_dir), path(sink_dir), val(n_source), val(n_sink), val(format), val(seed)

    output:
    path("*.fastq.gz")
 
    publishDir "${params.output_dir}/${contaminated_sample_name}", mode: 'copy'

    script:
    def ext = params.fastq_extension
    def pigz_cmd = "pigz -p ${params.pigz_threads}"

    if (format == 'single') {
        return """
        # Subsample source reads (single-end)
        ls ${source_dir}/*.${ext} | sort | xargs zcat -f | seqtk sample -s${seed} - ${n_source} > source.fastq

        # Subsample sink reads (single-end)
        ls ${sink_dir}/*.${ext} | sort | xargs zcat -f | seqtk sample -s${seed} - ${n_sink} > sink.fastq

        # Merge source and sink
        cat source.fastq sink.fastq | ${pigz_cmd} > ${contaminated_sample_name}.fastq.gz

        # Cleanup
        rm -f source.fastq sink.fastq
        """
    } else {
        return """
        # Subsample source reads (paired-end)
        ls ${source_dir}/*_1.${ext} | sort | xargs zcat -f | seqtk sample -s${seed} - ${n_source} > source_1.fastq
        ls ${source_dir}/*_2.${ext} | sort | xargs zcat -f | seqtk sample -s${seed} - ${n_source} > source_2.fastq

        # Subsample sink reads (paired-end)
        ls ${sink_dir}/*_1.${ext} | sort | xargs zcat -f | seqtk sample -s${seed} - ${n_sink} > sink_1.fastq
        ls ${sink_dir}/*_2.${ext} | sort | xargs zcat -f | seqtk sample -s${seed} - ${n_sink} > sink_2.fastq

        # Merge source and sink
        cat source_1.fastq sink_1.fastq | ${pigz_cmd} > ${contaminated_sample_name}_1.fastq.gz
        cat source_2.fastq sink_2.fastq | ${pigz_cmd} > ${contaminated_sample_name}_2.fastq.gz

        # Cleanup
        rm -f source_1.fastq source_2.fastq sink_1.fastq sink_2.fastq
        """
    }
}

workflow {
	simulate_contamination(contamination_tasks)
}

