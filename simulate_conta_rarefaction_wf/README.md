
# Nextflow Workflow for Creating Semi-Simulated Contaminated Samples by rarefying and combining gene counts
This Nextflow workflow generates semi-simulated contaminated samples by rarefying and combining gene counts from an existing gene count table.
For all contaminated samples, it generates species abundance tables compatible with [CroCoDeEL](https://github.com/metagenopolis/CroCoDeEL).

## Requirements

* [Nextflow](https://www.nextflow.io/)
* [Conda](https://www.anaconda.com/download/) (Anaconda or Miniconda)

Nextflow will automatically create the appropriate Conda environment on the first run.

### Tested with

* Nextflow 24.04.3 
* Conda 24.5.0

---

## Input Files

### Gene count table

The workflow requires a gene count table (rows: gene IDs, columns: sample names), such as one generated with the [`meteor2`](../meteor2_wf) workflow.

### Contamination Description Table

Prepare a tab-separated values (TSV) file listing the contaminated samples to generate, with the following columns:

| Column                          | Description                                                        |
| ------------------------------- | ------------------------------------------------------------------ |
| `contaminated_sample_name`      | Name of the contaminated sample to generate                        |
| `source_sample_name`            | Name of the contamination source sample (must exist in input data) |
| `sink_sample_name`              | Name of the sink sample (must exist in input data)                 |
| `num_reads_source`              | Number of reads in the source sample                               |
| `num_reads_contaminated_sample` | Number of reads in the contaminated sample                         |
| `contamination_rate`            | Contamination rate (float between 0 and 1)                         |

📄 Example table: [conta\_desc\_example.tsv](conta_desc_example.tsv)

---

## Running the Workflow

```bash
nextflow run simulate.nf \
  --project_name <output_files_prefix> \
  --gene_count_table <path_to_gene_count> \
  --conta_desc_table <path_to_description_table.tsv> \
  --output_dir <path_to_output_dir>
 ```

---

## Output Files

The output directory will contain the following file:
* `<project_name>.contaminated_samples.meteor2_species_ab_profiles.tsv`  
  This is a species abundance table, where rows represent species (as MSPs) and columns represent samples.  
  This file is compatible with CroCoDeEL.

---

