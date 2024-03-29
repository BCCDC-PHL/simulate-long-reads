# Simulate Long Reads
This pipeline is designed to simulate long sequencing reads from a set of reference genomes, for a variety of data quality and quantity.
The intention of this pipeline is to be used as a tool to create datasets that can be used as part of a validation experiment for other
microbial genomics analysis pipelines. When used in that application, the simulated datasets should be integrated into a larger validation plan
that incorporates real sequencing data that is derived from real sequencing runs on instuments that are identical or as similar as possible to those
that will be used in production.

After simulating the reads, they are mapped back onto the input genomes using [`minimap2`](https://github.com/lh3/minimap2). QC statistics are collected
from the alignment using [`qualimap`](https://github.com/scchess/Qualimap) and [`samtools`](https://github.com/samtools/samtools).

Some basic quality statistics are also collected on the simulated reads using [`fastp`](https://github.com/OpenGene/fastp).

The pipeline currently supports oxford nanopore reads, simulated using the [Badread](https://github.com/rrwick/Badread)
simulator, using its default error profile.

## Dependencies
Dependencies are listed in [environments.yml](environments/environment.yml). When using the `conda` profile, they will be installed automatically.
Other runtime environments such as docker, podman or singularity are not currently supported.

## Usage
To prepare a pipeline run, first collect your reference genome assemblies together in a directory. Assemblies should be in FASTA format and will
only be recognized for input if they use one of these filename extensions: `.fa`, `.fasta`, `.fna`.

To run the pipeline with default parameters:

```
nextflow run BCCDC-PHL/simulate-long-reads \
  -profile conda \
  --cache ~/.conda/envs \
  --assembly_input </path/to/assemblies> \
  --outdir </path/to/outdir>
```

Several parameters are available to control the quality and quantity of the simulated reads. The pipeline will automatically create simulated outputs at
eight mean coverage depths: 5x, 15x, 25x, 35x, 50x, 75x, 100x, 150x.

| Parameter               | Description                                                                                                                | Default |
|:------------------------|:---------------------------------------------------------------------------------------------------------------------------|--------:|
| `mean_read_length`      | Mean length of the reads to be simulated (bases).                                                                          | 15000   |
| `stdev_read_length`     | Standard deviation of the length of the reads to be simulated (bases).                                                     | 13000   |
| `percent_junk_reads`    | This percentage of reads will be low-complexity junk.                                                                      | 1       |
| `percent_random_reads`  | This percentage of reads will be random sequence.                                                                          | 1       |
| `percent_chimeras`      | Percentage at which separate fragments join together.                                                                      | 1       |
| `replicates`            | Number of replicates to generate for each genome, at each depth.                                                           | 1       |

### Introducing Contamination
Contamination can be introduced using the `--contaminants` parameter. The flag takes a `.csv` formatted file as an argument, with the following fields:

```
ID
ASSEMBLY
PROPORTION
```

...where `ID` is an identifier for the contaminant (avoid using spaces in the identifier), `ASSEMBLY` is a path to a FASTA-formatted genome that contaminant reads will be simulated from, and `PROPORTION` is a floating point number between 0 and 1 that determines the proportion of the sample that is contaminated by that genome. The sum of the values in the `PROPOTION` field should not be greater than 1.

## Outputs
The pipeline creates one output directory per set of simulated reads, below the directory provided for the `--outdir` parameter. The output directories are named
using the filename of the input reference genomes (excluding the file extension), with the additon of a 4-character string that is derived from the MD5 checksum of
the sample ID and the parameters used to generate the simulated reads.

For example, if we were to run this pipeline on a reference genome called [`ATCC-BAA-2787.fasta`](https://genomes.atcc.org/genomes/680bf0f0947a443c), including three contaminant genomes ([`ATCC-BAA-3053`](https://genomes.atcc.org/genomes/0c1563977c244589), [`ATCC-BAA-3038`](https://genomes.atcc.org/genomes/67f3a5f5558b4da1), [`ATCC-BAA-710`](https://genomes.atcc.org/genomes/100b594ab4114233)), and output the files  into a directory named `output`, the pipeline would produce the following outputs:

```
output
├── ATCC-BAA-2787-083b
│   ├── ATCC-BAA-2787-083b.bam
│   ├── ATCC-BAA-2787-083b.bam.bai
│   ├── ATCC-BAA-2787-083b_bamqc
│   ├── ATCC-BAA-2787-083b_fastp.csv
│   ├── ATCC-BAA-2787-083b_fastp.json
│   ├── ATCC-BAA-2787-083b_qualimap_bamqc_genome_results.csv
│   ├── ATCC-BAA-2787-083b_L.fastq.gz
│   ├── ATCC-BAA-2787-083b_read_simulation_parameters.csv
│   ├── ATCC-BAA-2787-083b_samtools_stats_summary.txt
│   └── contaminants
│       ├── ATCC-BAA-2787-083b-ATCC-BAA-3038_num_contaminant_reads.csv
│       ├── ATCC-BAA-2787-083b-ATCC-BAA-3053_num_contaminant_reads.csv
│       ├── ATCC-BAA-2787-083b-ATCC-BAA-710_num_contaminant_reads.csv
│       ├── ATCC-BAA-3038_contaminant_L.fastq.gz
│       ├── ATCC-BAA-3053_contaminant_L.fastq.gz
│       └── ATCC-BAA-710_contaminant_L.fastq.gz
├── ATCC-BAA-2787-096a
│   ├── ATCC-BAA-2787-096a.bam
│   ├── ATCC-BAA-2787-096a.bam.bai
│   ├── ATCC-BAA-2787-096a_bamqc
│   ├── ATCC-BAA-2787-096a_fastp.csv
│   ├── ATCC-BAA-2787-096a_fastp.json
│   ├── ATCC-BAA-2787-096a_qualimap_bamqc_genome_results.csv
│   ├── ATCC-BAA-2787-096a_L.fastq.gz
│   ├── ATCC-BAA-2787-096a_read_simulation_parameters.csv
│   ├── ATCC-BAA-2787-096a_samtools_stats_summary.txt
│   └── contaminants
│       ├── ATCC-BAA-2787-096a-ATCC-BAA-3038_num_contaminant_reads.csv
│       ├── ATCC-BAA-2787-096a-ATCC-BAA-3053_num_contaminant_reads.csv
│       ├── ATCC-BAA-2787-096a-ATCC-BAA-710_num_contaminant_reads.csv
│       ├── ATCC-BAA-3038_contaminant_L.fastq.gz
│       ├── ATCC-BAA-3053_contaminant_L.fastq.gz
│       └── ATCC-BAA-710_contaminant_L.fastq.gz
...
```

For both `fastp` and `qualimap`, the original outputs are parsed to create simplified `.csv`-formatted outputs.

The headers for the `_fastp.csv` output are:

```
sample_id
total_reads
total_bases
q20_bases
q30_bases
q20_rate
q30_rate
```

The headers for the `_qualimap_bamqc_genome_results.csv` output are:

```
sample_id
median_insert_size
mean_coverage
stdev_coverage
proportion_genome_covered_over_5x
proportion_genome_covered_over_10x
proportion_genome_covered_over_20x
proportion_genome_covered_over_30x
proportion_genome_covered_over_40x
proportion_genome_covered_over_50x
```

For each set of simulated reads, a record is made of which parameters were used to generate the reads. This record is stored in the `_read_simulation_parameters.csv` file and has
the headers:

```
sample_id
replicate
random_seed
fold_coverage
mean_read_length
stdev_read_length
junk_reads_percent
random_reads_percent
chimeras_percent
```

If contaminants are included, then a `contaminants` sub-directory will be created within each output directory. That directory includes files named with the simulated library ID and the contaminant ID, followed by `_num_contaminant_read_pairs.csv`. Those files have the following headers:

```
sample_id
contaminant_id
num_simulated_reads
num_contaminant_reads
target_contaminant_proportion
```

In addition to the `_num_contaminant_reads.csv` files, the reads that were used as contaminants are also included in the `contaminants` directory.
