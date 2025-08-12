nextflow.enable.dsl=2

// Mandatory parameters
if (!params.conta_desc_table)  error "Missing required parameter: --conta_desc_table"
if (!params.output_dir)  error "Missing required parameter: --output_dir"
if (!params.gene_count_table)  error "Missing required parameter: --gene_count_table"

// Optional/default parameters
params.R_env             = "${projectDir}/conda_envs/simulate_conta_rarefaction_env.yml"
params.simulation_script = "${projectDir}/Rscript/simulationContamination.R"

process simulate_contamination {
    conda "${params.R_env}"

    tag "${contaminated_sample_name}"

    input:
    tuple path(script), path(gene_count), val(contaminated_sample_name), val(source_name), val(sink_name), val(n_source), val(n_contaminated_sample), val(contamination_rate)

    output:
    path("*.tsv")
 
    script:
    """
    Rscript -e '
    if (!requireNamespace("momr", quietly=TRUE)) {
      devtools::install_git(url="https://forge.inrae.fr/metagenopolis/momr.git")
    }'

    Rscript $script $gene_count $contaminated_sample_name $source_name $sink_name $contamination_rate $n_contaminated_sample $n_source 
    """

}


process merge_files {
    conda "${params.R_env}"
    
    input:
    path(tsv_files)

    output:
    path("*tsv")

    publishDir "${params.output_dir}", mode: 'copy'

    script:
    """
    Rscript -e '
    library(dplyr)
    library(readr)
    
    args <- commandArgs(trailingOnly = TRUE)
      tables <- lapply(args, read_tsv)
      merged <- Reduce(function(x, y) full_join(x, y, by = colnames(x)[1]), tables)
      write_tsv(merged, "${params.project_name}.contaminated_samples.meteor2_species_ab_profiles.tsv")
    ' ${tsv_files.join(' ')}
    """

    

}

workflow {

    Channel
    .fromPath(params.gene_count_table)
    .set { gene_count_table }

    Channel
    .fromPath(params.conta_desc_table)
    .splitCsv(header:true, sep:'\t')
    .combine(gene_count_table)
    .map { row, gene_count_file ->
        tuple(
            file(params.simulation_script), 
            file(gene_count_file),
            row.contaminated_sample_name,
            row.source_sample_name,
            row.sink_sample_name,
            row.num_reads_source.toInteger(),
            row.num_reads_contaminated_sample.toInteger(),
            row.contamination_rate
        )
    }
    .set { contamination_tasks }

	contamination_res = simulate_contamination(contamination_tasks)

    merge_files(contamination_res.collect())
}