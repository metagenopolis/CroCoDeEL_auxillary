// Mandatory parameters
if (!params.sequencing_data_dir)  error "Missing required parameter: --sequencing_data_dir"
if (!params.output_dir)           error "Missing required parameter: --output_dir"
if (!params.project_name)         error "Missing required parameter: --project_name"

// Optional/default parameters
params.meteor_env = "conda_envs/meteor_env.yml"
params.gene_catalogue = "hs_10_4_gut"
params.fastq_extension = '.fastq.gz'
params.cpus = 8

def fastq_file_glob = "*${params.fastq_extension}"
def fastq_file_regex = /^(.+?)(?:_([12]))?\Q${params.fastq_extension}\E$/



process meteor_download {
    conda "${params.meteor_env}"

    tag { gene_catalogue }

    input:
    val(gene_catalogue)

    output:
    path("${gene_catalogue}")

    script:
    """
    meteor download \\
        -i ${gene_catalogue} \\
        -o ./
    """
}

process index_fastq_files {
    tag "${sample_name}/${fastq_file}"

    input:
    tuple val(sample_name), val(run_name), val(tag), path(fastq_file)

    output:
    tuple val(sample_name), path(fastq_file), path("*_census_stage_0.json")

    script:
    """
    cat <<EOF > "${run_name}_${tag}_census_stage_0.json"
    {
      "sample_info": {
        "sample_name": "${sample_name}",
        "tag": "${tag}",
        "full_sample_name": "${run_name}"
      },
      "sample_file": {
        "fastq_file": "${fastq_file}"
      }
    }
    EOF
    """
}

process meteor_mapping {
    cpus params.cpus
    memory "100G"
    conda "${params.meteor_env}"

    tag { sample_name }

    input:
    tuple val(sample_name), path(fastq_files), path(json_files), path(gene_catalogue)

    output:
    tuple val(sample_name), path("mapping"), path(gene_catalogue)

    script:
    """
    meteor mapping \\
        -i ./ \\
        -r ${gene_catalogue} \\
        -t ${params.cpus} \\
        -o mapping
    """
}

process meteor_profile {
    cpus params.cpus
    memory "10G"
    conda "${params.meteor_env}"

    tag { sample_name }

    input:
    tuple val(sample_name), path(mapping_dir), path(gene_catalogue)

    output:
    path("profile/*")

    script:
    """
    meteor profile \\
        -i ${mapping_dir}/${sample_name} \\
        -r ${gene_catalogue} \\
        -o profile
    """
}

process meteor_merge {
    cpus params.cpus
    memory "10G"
    conda "${params.meteor_env}"

    publishDir "${params.output_dir}", mode: 'copy'

    input:
    path(profile_dirs)
    path(gene_catalogue)
    val(project_name)

    output:
    path { "*.tsv" }

    script:
    """
        PREFIX="${project_name}_vs_${gene_catalogue}.meteor2"
	meteor merge \\
        -i ./ \\
        -r ${gene_catalogue} \\
        -p ${project_name} \\
        -g \\
        -o merge

        mv merge/${project_name}_raw.tsv \\
        \${PREFIX}_raw_gene_profiles.tsv 

        mv merge/${project_name}_msp.tsv \\
        \${PREFIX}_species_ab_profiles.tsv 

        mv merge/${project_name}_msp_taxonomy.tsv \\
        \${PREFIX}_species_tax.tsv 
    """
}


workflow {
    Channel.value(params.gene_catalogue)
        .set { gene_catalogue_ch }
    gene_catalogue_path = meteor_download(gene_catalogue_ch) 

    Channel
        .fromPath("${params.sequencing_data_dir}/**/${fastq_file_glob}")
        .map { fastq_file ->
            def file_name = fastq_file.getName()
            def sample_dir = fastq_file.getParent()
            def sample_name = sample_dir.getName()
            def match = file_name =~ fastq_file_regex
            if (!match) error "Filename doesn't match expected pattern: ${file_name}"
    
            def run_name = match[0][1]
            def tag = match[0][2] ?: "single"
    
            return [sample_name, run_name, tag, fastq_file]
        }
        .set { fastq_files_infos }
    indexed_fastq_files = index_fastq_files(fastq_files_infos)

    indexed_fastq_files
        .groupTuple()
	.combine(gene_catalogue_path)
        .set { grouped_indexed_fastq_files}
    mapping_res = meteor_mapping(grouped_indexed_fastq_files)

    profile_res = meteor_profile(mapping_res)

    Channel.value(params.project_name)
        .set { project_name_ch }
    meteor_merge(profile_res.collect(), gene_catalogue_path, project_name_ch)
}