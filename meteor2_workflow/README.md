
# Nextflow Workflow for metagenomic profiling with Meteor2

## Overview

This is a Nextflow ([https://www.nextflow.io/](https://www.nextflow.io/)) workflow to run **Meteor2** ([https://github.com/metagenopolis/Meteor](https://github.com/metagenopolis/Meteor)) on short-read, paired-end metagenomic sequencing data.
It produces gene and species abundance tables among other outputs, enabling downstream microbial profiling and comparative analysis.

> ⚠️ Note: This workflow uses a custom Meteor2 Conda package ([https://anaconda.org/fplazaonate/Meteor](https://anaconda.org/fplazaonate/Meteor)) to ensure full compatibility with Nextflow.

---

## Requirements

* Nextflow ([https://www.nextflow.io/](https://www.nextflow.io/))
* Conda (Anaconda or Miniconda) ([https://www.anaconda.com/download/](https://www.anaconda.com/download/))

Nextflow will automatically create the appropriate Conda environment and install Meteor2 on the first run.

### Tested with

* Nextflow 25.04.4
* Conda 23.3.1

---

## Input Data Structure

Organize your metagenomic data as follows:

* Store all sequencing data in a single directory.
* Each sample must be placed in its own subdirectory.
* Each subdirectory should contain one or more FASTQ files.

> ⚠️ File naming must follow the convention with `_1` and `_2` indicating forward and reverse reads, respectively.

Example:

sequencing\_data/
├── sample1/
│   ├── runA\_1.fastq.gz
│   ├── runA\_2.fastq.gz
│   ├── runB\_1.fastq.gz
│   ├── runB\_2.fastq.gz
├── sample2/
│   ├── runC\_1.fastq.gz
│   ├── runC\_2.fastq.gz
│   ├── runD\_1.fastq.gz
│   ├── runD\_2.fastq.gz

---

## Run the Workflow

The workflow `meteor.nf` performs the following steps:

* Downloads the human gut gene catalogue (`meteor_download`)
* Indexes the FASTQ files for downstream analysis (`index_fastq_files`)
* Maps reads to the gene catalogue (`meteor mapping`)
* Generates gene and species abundance tables (`meteor profile`)
* Merges results across all processed samples (`meteor merge`)

### Basic command

```bash
nextflow run meteor.nf \
  --project_name <output_files_prefix> \
  --sequencing_data_dir <path_to_sequencing_data> \
  --output_dir <output_dir>
 ```
### Optional Parameters

You can specify a different gene catalogue and FASTQ file extension if needed.
Below, we use the mouse gut gene catalogue and uncompressed FASTQ files.

```bash
nextflow run meteor.nf \
  --project_name <output_files_prefix> \
  --sequencing_data_dir <path_to_sequencing_data> \
  --output_dir <output_dir>
  --catalogue_name mm_5_0_gut \
  --fastq_extension '.fastq'
 ```

#### Defaults:

* \--catalogue\_name: hs\_10\_4\_gut
* \--fastq\_extension: .fastq.gz

---

## Output Files

The output directory contains the following main files:

* `<project_name>_vs_<catalogue_name>.meteor2_raw_gene_profiles.tsv`
  Gene count table (rows: gene IDs, columns: samples). Useful for generating simulated contaminated samples by rarefaction.

* `<project_name>_vs_<catalogue_name>.meteor2_species_ab_profiles.tsv`
  Species abundance table based on estimated genome coverage (rows: MSPs, columns: samples). This file is compatible with CroCoDeEL.

* `<project_name>_vs_<catalogue_name>.meteor2_species_tax.tsv`
  GTDB taxonomy table for the identified species.

---

