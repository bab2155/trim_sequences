# summarize_reads.py

Script that BLASTs dataset against primer and adapter sequences, creates summary statistics, and trims primer and adapter sequences.

## Setup

Requires `python 3.6.1` or greater

Install dependencies:

`pip install -r requirements.txt`

## Usage

`python ./summarize_reads.py <path to config.ini file>`

A sample config.ini file was included that was used to generate the files in the `output` folder.

## Config file

A config file is required to define relevant paths, as well as the primer and adapter sequences. The following template shows required section and variable names, with a short description.

**Template**

```
[INPUT]

# Primer sequence
PRIMER_SEQUENCE = CGCCGTTTCCCAGTAGGTCTC

# Adapter sequence
ADAPTER_SEQUENCE = ACTGAGTGGGAGGCAAGGCACACAGGGGATAGG

# Path to FASTA file of dataset
FASTA_PATH = ./input/test.fna

# Path to QUAL file of dataset
QUAL_PATH = ./input/test.qual

[BLAST]

# Path to local blastn executable
BLASTN_PATH = tools/ncbi-blast-2.8.1+/bin/blastn

# Path to local makeblastb executable
MAKEDB_PATH = tools/ncbi-blast-2.8.1+/bin/makeblastdb

[OUTPUT]

# Output directory
OUTPUT_DIR = ./output/

# tmp directory
TMP_DIR = ./tmp/
```

## Output

The script creates 4 output files

1. `blast_output.m8`

Blast output file in m8 format

2. `output.fasta`

Fasta file containing reads greater than 100bp, average read quality scores greater than 20, primers and adaptors trimmed. Primers are only trimmed if they exist at the start of the read. Adapters are trimmed if they exist anywhere in the read.

3. `output.tsv`

Tab de-limited text file containing the read identifiers along with the starting and end positions of the primer or adaptor sequences. Only reads that contain either the **whole** primer or **whole** adapter sequence are included.

4. `summary.tsv`

Provides the following summary statistics

    1) Total number of reads in the dataset.
    2) Total number of reads greater than 100 bp.
    3) Total number of reads with average quality scores greater than 20.
    4) Total number of reads with primer sequences.
    5) Total number of reads with adaptor sequences.
    6) Total number of reads with both primer and adaptor sequences.

NOTE: Reads that contain primer or adapter sequences anywhere in the read were counted towards these numbers.

## Notes on trimming

Primers were trimmed only if they were found in the beginning of the read.

Adapters found anywhere in the read were trimmed from the first base of the adapter to the end of the read.

Both primers and adapters were only trimmed if the **full, uninterrupted** sequence was found in the read.

