import os
import sys
import argparse
import configparser
from typing import List
from statistics import mean
from collections import namedtuple

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Blast.Applications import NcbiblastnCommandline, NcbimakeblastdbCommandline


###############
# START FUNCTIONS, CLASSES, AND SUBCLASSES
###############

# Header of BLAST output in m8 format
M8_HEADER = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
             'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']

# namedtuple to handle each line of BLAST output
BlastEntry = namedtuple('BlastEntry', M8_HEADER)


def make_records(file_path: str, file_format: str) -> List[SeqIO.SeqRecord]:
    '''Reads FASTA or QUAL files into structured format

    Input parameters:
        file_path (str): path to FASTA or QUAL file
        file_format (str): file format descriptor. Either 'fasta' or 'qual'

    Returns:
        (List[SeqIO.SeqRecord]): list of records in structured format
    '''

    out = []
    for record in SeqIO.parse(file_path, file_format):
        out.append(record)

    return out


def count_reads(records: List[SeqIO.SeqRecord]) -> int:
    '''Count all unique identifiers from FASTA records

    Input params:
        records (List[SeqIO.SeqRecord]): List of records from FASTA file

    Returns:
        (int): read count from records
    '''

    return len(set([record.id for record in records]))


def filter_reads_by_length(records: List[SeqIO.SeqRecord], length: int) -> List[SeqIO.SeqRecord]:
    '''Get length of each read string and keep only those greater than length

    Input params:
        records (List[SeqIO.SeqRecord]: list of records from FASTA file
        length (int): read length threshold which at or below reads will be removed

    Returns
        (List[SeqIO.SeqRecord]): list of records with read lengths above length
    '''

    return [record for record in records if len(record) > length]


def filter_quality_scores(qual: List[SeqIO.SeqRecord], threshold: float) -> List[SeqIO.SeqRecord]:
    '''Calculate average quality score for each read and keep those greater than threshold

    Input parameters:
        qual (List[SeqIO.SeqRecord]): list of records from QUAL file
        threshold (float): threshold of quality scores which at or below will be filtered out

    Returns:
         (List[SeqIO.SeqRecord]): QUAL records with average quality scores above threshold
    '''

    def get_average_score(qual_record: SeqIO.SeqRecord) -> float:
        '''Calculates average quality score from QUAL record.

        Input paramters:
            qual_record (SeqIO.SeqRecord): single record from QUAL file

        Returns:
            (float): average quality score of the read.
                     -1 if record does not contain "phred_quality" key.
        '''

        return mean(qual_record.letter_annotations.get('phred_quality', [-1]))

    return [record for record in qual if get_average_score(record) > threshold]


def get_reads_containing_entire_sequence(blast: List[BlastEntry], query_id: str,
                                         query_sequence: str) -> List[BlastEntry]:
    '''Returns reads that contain the entire query_sequence anywhere within the read.

    Input parameters:
        blast (List[BlastEntry]): list of BLAST entries
        query_id (str): subject sequence ID from BLAST output
        query_sequence (str): sequence being BLASTed against

    Returns:
        (List[BlastEntry]): list of BLAST entries that contain the entire query_sequence anywhere in the read.
    '''

    passing_reads = [entry for entry in blast if
                     entry.sseqid == query_id and
                     entry.length == str(len(query_sequence)) and
                     float(entry.pident) == 100]

    return passing_reads


def parse_blast_output(blast_out: str) -> List[BlastEntry]:
    '''Parse BLAST output into structured format

    Input parameters:
        blast_output (str): raw blast output string in m8 format

    Returns:
        out (List[BlastEntry]): list of namedtuples representing each line of BLAST output

    '''

    out = blast_out.split('\n')
    out = [s.split('\t') for s in out]
    out = [BlastEntry._make(s) for s in out if len(s) == 12]

    return out


def intersect_fasta_qual(fasta: List[SeqIO.SeqRecord], qual: List[SeqIO.SeqRecord]) -> List[SeqIO.SeqRecord]:
    '''Filter FASTA records by intersecting IDs from FASTA records and QUAL records

    Input parameters:
        fasta (List[SeqIO.SeqRecord]): FASTA records
        qual (List[SeqIO.SeqRecord]): QUAL records

    Returns:
        out (List[SeqIO.SeqRecord]): FASTA records that have IDs in both input FASTA and QUAL records
    '''

    passing_read_ids = set([read.id for read in fasta])
    passing_qual_score_ids = set([record.id for record in qual])
    filtered_ids = passing_read_ids.intersection(passing_qual_score_ids)

    # Subset FASTA records by passing ids
    out = [record for record in fasta_records if record.id in filtered_ids]

    return out

###############
# END FUNCTIONS, CLASSES, AND SUBCLASSES
###############


###############
# Set up argument parser.
###############

# Because there is only one parameter (path to config file), it's only intention was to handle
# displaying usage summary to the user.
parser = argparse.ArgumentParser(description='Summarize various metrics from FASTA and QUAL files. Used particularly to '
                                             'answer specific questions posed.')
parser.add_argument('config_path', help='path to config file')
args = parser.parse_args()

###############
# Parse config file
###############

# Load in config file
config_path = sys.argv[1]
config = configparser.ConfigParser()
config.read(config_path)

# Define output and tmp directories. Make the directories if they don't exist.
OUT_DIR = config['OUTPUT'].get('OUTPUT_DIR')
if not os.path.exists(OUT_DIR):
    os.makedirs(OUT_DIR)

TMP_DIR = config['OUTPUT'].get('TMP_DIR', './tmp/')
if not os.path.exists(TMP_DIR):
    os.makedirs(TMP_DIR)

PRIMER_SEQUENCE = config['INPUT'].get('PRIMER_SEQUENCE') # Primer sequence
ADAPTER_SEQUENCE = config['INPUT'].get('ADAPTER_SEQUENCE') # Adapter sequence
FASTA_PATH = config['INPUT'].get('FASTA_PATH') # Path to input FASTA file
QUAL_PATH = config['INPUT'].get('QUAL_PATH') # Path to corresponding input QUAL file
SEQUENCE_FASTA_PATH = os.path.join(TMP_DIR, 'sequences.fasta') # FASTA file made from primer and adapter sequences.
BLASTN_PATH = config['BLAST'].get('BLASTN_PATH') # Path to blastn executable
MAKEDB_PATH = config['BLAST'].get('MAKEDB_PATH') # Path to makeblastdb executable
BLAST_OUTPUT_PATH = os.path.join(OUT_DIR, 'blast_output.m8') # BLAST m8 output path
OUTPUT_FASTA_PATH = os.path.join(OUT_DIR, 'output.fasta') # Output FASTA path
OUT_TABFILE_PATH = os.path.join(OUT_DIR, 'output.tsv') # Output tab-delimited file path
SUMMARY_PATH = os.path.join(OUT_DIR, 'summary.tsv') # Path to answered questions

###############
# Load FASTA and QUAL files
###############

print('Loading FASTA and QUAL files...', end='')

fasta_records = make_records(FASTA_PATH, 'fasta')
qual_records = make_records(QUAL_PATH, 'qual')

print('Done!')


###############
# Run BLAST
###############

print('Running BLAST...', end='')

# Convert primer and adapter sequences to FASTA file
primer_record = SeqIO.SeqRecord(Seq(PRIMER_SEQUENCE), id='primer', description='')
adapter_record = SeqIO.SeqRecord(Seq(ADAPTER_SEQUENCE), id='adapter', description='')
sequence_records = [primer_record, adapter_record]
SeqIO.write(sequence_records, SEQUENCE_FASTA_PATH, 'fasta')

# Make db from primer/adapter FASTA file
makedb = NcbimakeblastdbCommandline(cmd=MAKEDB_PATH, input_file=SEQUENCE_FASTA_PATH, input_type='fasta', parse_seqids=True, dbtype='nucl', out=SEQUENCE_FASTA_PATH)()

# Run blastn
blast_output = NcbiblastnCommandline(cmd=BLASTN_PATH, db=SEQUENCE_FASTA_PATH, query=FASTA_PATH, outfmt=6, task='blastn', perc_identity=100)()[0]

# Parse BLAST lines into structured format
blast_entries = parse_blast_output(blast_output)

print('Done!')

###############
# Answer supplied questions
###############

print('Generating answers to questions...', end='')

# Initialize list containing all question and answer tuples. This will later be written to a file.
output = []


# 1) Total number of reads in the dataset.
total_reads = count_reads(fasta_records)
output.append(('Total number of reads in the dataset', total_reads))


# 2) Total number of reads greater than 100 bp.
READ_LENGTH_THRESHOLD = 100

reads_gt100 = filter_reads_by_length(fasta_records, READ_LENGTH_THRESHOLD)
output.append(('Total number of reads greater than 100 bp', len(reads_gt100)))


# 3) Total number of reads with average quality scores greater than 20.
QUALITY_SCORE_THRESHOLD = float(20)

filtered_quality_records = filter_quality_scores(qual_records, QUALITY_SCORE_THRESHOLD)
output.append(('Total number of reads with average quality scores greater than 20', len(filtered_quality_records)))


# 4) Total number of reads with primer sequences.
reads_with_primer_sequence = get_reads_containing_entire_sequence(blast_entries, 'primer', PRIMER_SEQUENCE)
output.append(('Total number of reads with primer sequences', len(reads_with_primer_sequence)))


# 5) Total number of reads with adaptor sequences.
reads_with_adapter_sequence = get_reads_containing_entire_sequence(blast_entries, 'adapter', ADAPTER_SEQUENCE)
output.append(('Total number of reads with adaptor sequences', len(reads_with_adapter_sequence)))


# 6) Total number of reads with both primer and adaptor sequences.
if any([not reads_with_primer_sequence, not reads_with_adapter_sequence]):
    # If either read list is empty, return empty
    reads_with_both_sequences = []
else:
    # Otherwise, find intersection of read IDs from the set of read IDs containing primer sequences and the set of read
    # IDs containing adapter sequences
    reads_with_both_sequences = set([read.qseqid for read in reads_with_primer_sequence]).intersection(
        set([read.qseqid for read in reads_with_adapter_sequence]))

output.append(('Total number of reads with both primer and adaptor sequences', len(reads_with_both_sequences)))


# Write answers to all questions to file
with open(SUMMARY_PATH, 'w') as f:
    for line in output:
        f.write('{}\t{}\n'.format(line[0], line[1]))

print('Done!')

###############
# 1) Blast output file in m8 format.
###############

print('Writing BLAST m8 output...', end='')

with open(BLAST_OUTPUT_PATH, 'w') as f:
    f.write(blast_output)

print('Done!')

###############
# 2) Fasta file containing reads greater than 100bp, average read quality scores greater than 20, primers and adaptors
# trimmed.
###############

print('Filtering FASTA file by read length and quality scores...', end='')

filtered_fasta = intersect_fasta_qual(reads_gt100, filtered_quality_records)

print('Done!')

print('Trimming primer and adapter sequences...', end='')

# Trim primers and adapters
for read in filtered_fasta:
    if read.seq.startswith(PRIMER_SEQUENCE):
        read.seq = read.seq[len(PRIMER_SEQUENCE):]

    index = read.seq.find(ADAPTER_SEQUENCE)
    if index != -1:
        read.seq = read.seq[:index]

SeqIO.write(filtered_fasta, OUTPUT_FASTA_PATH, 'fasta')

print('Done!')

###############
# 3) Tab de-limited text file containing the read identifiers along with the starting and end positions of the primer
# or adaptor sequences.
###############

print('Generating tab-delimited output file of query sequence positions...', end='')

tab_output = []
reads_with_any_sequence = reads_with_primer_sequence + reads_with_adapter_sequence
for entry in reads_with_any_sequence:
        tab_output.append((entry.qseqid, entry.sseqid, entry.qstart, entry.qend))

with open(OUT_TABFILE_PATH, 'w') as f:
    for line in tab_output:
        f.write('{}\n'.format('\t'.join(line)))

print('Done!')

# Print output summary
print('\nOutput paths:')
print('================')
print('Answers to provided questions: {}'.format(SUMMARY_PATH))
print('BLAST output in m8 format: {}'.format(BLAST_OUTPUT_PATH))
print('Filtered FASTA file: {}'.format(OUTPUT_FASTA_PATH))
print('Tab-delimited summary output: {}'.format(OUT_TABFILE_PATH))
