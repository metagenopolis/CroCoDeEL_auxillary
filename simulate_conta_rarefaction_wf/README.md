
# Nextflow Workflow for Creating Semi-Simulated Contaminated Samples by rarefaction

## Requirements

* Nextflow ([https://www.nextflow.io/](https://www.nextflow.io/))
* Conda (Anaconda or Miniconda) ([https://www.anaconda.com/download/](https://www.anaconda.com/download/))

Nextflow will automatically create the appropriate Conda environment and install Meteor2 on the first run.

### Tested with

* Nextflow 24.04.3 
* Conda 24.5.0

---

## Input Data Structure

Organize your metagenomic data as follows:

* Store all sequencing data in a single directory.
* Each sample must be placed in its own subdirectory.
* Each subdirectory should contain one or more FASTQ files.

> ⚠️ File naming must follow the convention with `_1` and `_2` indicating forward and reverse reads, respectively.
> ⚠️ If these suffixes are not present, files will be treated as single-end reads.

Example:

sequencing\_data/  
├── sample1\_pe/  
│   ├── runA\_1.fastq.gz  
│   ├── runA\_2.fastq.gz  
│   ├── runB\_1.fastq.gz  
│   ├── runB\_2.fastq.gz  
├── sample2\_pe/  
│   ├── runC\_1.fastq.gz  
│   ├── runC\_2.fastq.gz  
│   ├── runD\_1.fastq.gz  
│   ├── runD\_2.fastq.gz  
├── sample3\_se/  
│   ├── runE.fastq.gz  
│   ├── runF.fastq.gz  

---

## Contamination Description Table

Prepare a tab-separated values (TSV) file listing the contaminated samples to generate, with the following columns:

| Column                     | Description                                                                |
| -------------------------- | -------------------------------------------------------------------------- |
| `contaminated_sample_name` | Name of the contaminated sample to be created                              |
| `source_sample_name`       | Name of the sample used as contamination source (must exist in input data) |
| `sink_sample_name`         | Name of the sample acting as the sink (must exist in input data)           |
| `num_reads_source`    | Number of reads in the source                                    |
| `num_reads_contaminated_sample`      | Number of reads in the contaminated sample                                      |
| `contamination_rate`           | Contamination rate (between 0 and 1)   |

📄 Example table: [conta\_desc\_example.tsv](conta_desc_example.tsv)

---

### Basic command

```bash
nextflow run simulate.nf \
  --project_name <output_files_prefix> \
  --gene_count_table <path_to_gene_count> \
  --conta_desc_table <path_to_description_table.tsv> \
  --output_dir <path_to_output_dir>
 ```
### Optional Parameters

You can use the same optional parameters as the Meteor workflow.

---

## Output Files

The output directory contains the following main file:
* `<project_name>.contaminated_samples.meteor2_species_ab_profiles.tsv`
  Species abundance table based on estimated genome coverage (rows: MSPs, columns: samples). This file is compatible with CroCoDeEL.  

As well as the Meteor workflow output files.

---

