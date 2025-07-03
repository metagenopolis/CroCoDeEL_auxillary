# Nextflow Workflow for Creating Semi-Simulated Contaminated Samples from Real Metagenomic Data

This Nextflow pipeline simulates contamination between metagenomic samples by subsampling and mixing reads from real datasets.

---

## Requirements

* [Nextflow](https://www.nextflow.io/)
* [Conda](https://www.anaconda.com/download/) (Anaconda or Miniconda)

Nextflow will automatically create the required Conda environment and install dependencies (`seqtk`, `pigz`) on first run.

### Tested With

* Nextflow 25.04.4
* Conda 23.3.1

---

## Input Data Structure

Organize your sequencing data as follows:

* Place all sequencing data in a single directory.
* Each sample must have its own subdirectory.
* Each subdirectory must contain one or more FASTQ files.

> ⚠️ For paired-end data, file names must include `_1` and `_2` to indicate forward and reverse reads, respectively.  
> ⚠️ Each sample must contain **either** single-end **or** paired-end data, not both.

### Example

```
sequencing_data/  
├── sample1_pe/  
│   ├── runA_1.fastq.gz  
│   ├── runA_2.fastq.gz  
│   ├── runB_1.fastq.gz  
│   ├── runB_2.fastq.gz  
├── sample2_pe/  
│   ├── runC_1.fastq.gz  
│   ├── runC_2.fastq.gz  
│   ├── runD_1.fastq.gz  
│   ├── runD_2.fastq.gz  
├── sample3_se/  
│   ├── runE.fastq.gz  
│   ├── runF.fastq.gz  
```

---

## Contamination Description Table

Prepare a tab-separated values (TSV) file listing the contaminated samples to generate, with the following columns:

| Column                     | Description                                                                |
| -------------------------- | -------------------------------------------------------------------------- |
| `contaminated_sample_name` | Name of the contaminated sample to be created                              |
| `source_sample_name`       | Name of the sample used as contamination source (must exist in input data) |
| `sink_sample_name`         | Name of the sample acting as the sink (must exist in input data)           |
| `num_reads_from_source`    | Number of reads to draw from the source                                    |
| `num_reads_from_sink`      | Number of reads to draw from the sink                                      |
| `library_layout`           | `"paired"` or `"single"` — both source and sink must use the same layout   |

📄 Example table: [conta\_desc\_example.tsv](conta_desc_example.tsv)

---

## Running the Workflow

The main script is `simulate_conta_read_subsampling.nf`. For each row in the contamination table, it will:

1. Randomly subsample reads from the source
2. Randomly subsample reads from the sink
3. Merge reads from both
4. Output the contaminated sample in gzip-compressed FASTQ format

### Basic Command

```bash
nextflow run simulate_conta_read_subsampling.nf \
  --sequencing_data_dir <path_to_sequencing_data> \
  --conta_desc_table <path_to_description_table.tsv> \
  --output_dir <path_to_output_directory>
```

### Optional Parameters

You can customize the expected FASTQ extension if your files differ from the default.

```bash
--fastq_extension <extension>   # Default: ".fastq.gz"
```
---

## Output Structure

The specified output directory will contain one subdirectory per contaminated sample. Each subdirectory contains:

* One FASTQ file for single-end samples
* Two FASTQ files (`_1` and `_2`) for paired-end samples

