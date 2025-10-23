#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# Author: Skyler Kuhn

# Standard Library
from __future__ import print_function
from textwrap import dedent
from datetime import datetime
import argparse, gzip, os, sys

# Constants
# Usage and help section
_NAME = os.path.basename(sys.argv[0])
_HELP = dedent("""
@Usage:
    $ ./{0} [-h] [--version] \\
            [--split-strain-fasta-ids]
            --input-blast-result PER_GENE_BLASTN_RESULTS \\
            --input-strain-fasta BLASTN_FASTA_DB \\
            --input-strain-ids STRAIN_SEQUENCE_IDS \\
            --output-strain-fasta OUTPUT_FASTA

@About:
    Given a parsed output and filtered blastn output
    file for a given gene, a file containing all the
    genome acession IDs for a strain of interest (i.e
    MSSA, MRSA-mecI-negative, MRSA-mecI-positive),
    and a FASTA file containing all the sequences in
    the blastn database, this script will create a
    strain-specific FASTA file containing the top
    hit's alignment sequence.

@Required:
    -i, --input-blast-result PER_GENE_BLASTN_RESULTS
        Per-gene, filtered input BLASTN results.
    -f, --input-strain-fasta BLASTN_FASTA_DB
        FASTA file used to create the BLASTN database.
        This contains all the genomic sequences for
        Staphylococcus aureus that were pulled from NCBI.
    -s, --input-strain-ids STRAIN_SEQUENCE_IDS
        File containing a single column with all
        of the genomic accession IDs associated
        with the strain of interest. The strains
        of interest are:
         • MSSA (mecA - strains)
         • MRSA (mecA + and mecI + strains)
         • MRSA (mecA + and mecI - strains)
        This can be parsed from the gene_dection sheet
        in the deliveried XLSX spreadsheet. The IDs in this
        file should match sequence IDs in the BLASTN_FASTA_DB.
    -o, --output-strain-fasta OUTPUT_FASTA
        Output FASTA file to write the gene aligment
        sequences for the strain of interest.
@Options:
    -d, --split-strain-fasta-ids
                   Split the sequence IDs in the input
                   strain FASTA file on the '|' character.
                   This is useful if the sequence IDs in
                   the FASTA file contain additional info.
    -h, --help     Shows this help message and exits.
    -v, --version  Prints the version and exits.

@Example:
    $ ./{0} \\
          --split-strain-fasta-ids \\
          --input-blast-result MecA_blastn_results.tsv \\
          --input-strain-ids MSSA_seqids.txt \\
          --input-strain-fasta staphylococcus_aureus_combined.fa \\
          --output-strain-fasta MSSA_mecA_aligned_sequence.fa
""".format(_NAME)
)

# Semantic version
_VERISON = '1.0.0'


# Classes and helper functions
class InvalidBasePairError(Exception):
    """Raised when a nucleotide cannot be mapped to a known complement base pair.
    This may occur when the provided sequence contains unknown/invalid base pairs.
    In this scenario, the corresponding complementary basepair cannot be determined.

    @attributes:
        sequence: Input DNA sequence which caused the error
        codon:    The exact nucleotide/basepair that cannot be mapped.
    """
    def __init__(self, sequence, nucleotide):
        self.sequence = sequence
        self.nucleotide = nucleotide
        self.message = """Error: Invalid nucleotide, '{}', within DNA sequence to reverse complement!
            └── Please view the provided the provided DNA sequence to complement:
                > {}""".format(self.nucleotide, self.sequence)
        super(Exception, self).__init__(self.message)

    def __str__(self):
        return "{} -> {}".format(self.nucleotide, self.message)


def err(*message, **kwargs):
    """Prints any provided args to standard error.
    kwargs can be provided to modify print functions
    behavior.
    @param message <any>:
        Values printed to standard error
    @params kwargs <print()>
        Key words to modify print function behavior
    """
    print(*message, file=sys.stderr, **kwargs)



def fatal(*message, **kwargs):
    """Prints any provided args to standard error
    and exits with an exit code of 1.
    @param message <any>:
        Values printed to standard error
    @params kwargs <print()>
        Key words to modify print function behavior
    """
    err(*message, **kwargs)
    sys.exit(1)


def timestamp(format="%Y-%m-%d %H:%M:%S"):
    """Returns a formatted timestamp string
    for the current time.
    @param format <str>:
        Format string for the timestamp, default:
        "%Y-%m-%d %H:%M:%S" which is equivalent to
        "2023-10-01 12:00:00" for example.
    @return <str>:
        Formatted timestamp string, i.e. "2023-10-01 12:00:00"
    """
    return datetime.now().strftime(format)


def log(*message):
    """Logs a message to standard output with a timestamp.
    @param message <any>:
        Values printed to log
    """
    print("[{0}] {1}".format(
        timestamp(),
        " ".join([str(m) for m in message]))
    )


def check_permissions(parser, path, *args, **kwargs):
    """Checks permissions using os.access() to see the
    user is authorized to access a file/directory. Checks
    for existence, read, write and execute via args:
        • os.F_OK (tests existence)
        • os.R_OK (tests read)
        • os.W_OK (tests write)
        • os.X_OK (tests exec)
    @param parser <argparse.ArgumentParser() object>:
        Argparse parser object
    @param path <str>:
        Name of path to check
    @param args <any>:
        Positional args to pass to os.access()
    @param kwargs <any>:
        Named kwargs to pass to os.access()
    @return path <str>:
        Returns absolute path if it exists and the
        checked permssions are setup are correct.
    """
    if not os.path.exists(path):
        parser.error(
            "Path '{}' does not exists! Failed to provide vaild input.".format(path)
        )
    if not os.access(path, *args, **kwargs):
        parser.error(
            "Path '{}' exists, but cannot read path due to permissions!".format(path)
        )
    return os.path.abspath(path)


def stripped(s):
    """Cleans string to remove quotes
    @param s <str>:
        String to remove quotes or clean
    @return s <str>:
        Cleaned string with quotes removed
    """
    return s.strip('"').strip("'").lstrip().rstrip()


def parse_cli_arguments():
    """Parses command line arguments and returns
    an argparse.parse_args object.
    @return <argparse.parse_args()>:
        Parsed command line arguments
    """
    parser = argparse.ArgumentParser(
        add_help=False,
        description=_HELP,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage = argparse.SUPPRESS,
    )
    # Input Per-gene, filtered (i.e alignment length 
    # fraction and pidentity) BLASTN file to parse
    parser.add_argument(
        '-i', '--input-blast-result',
        type = lambda file: check_permissions(parser, file, os.R_OK), 
        required=True,
        help=argparse.SUPPRESS
    )
    # FASTA file used to build the BLAST database,
    # this should contain all the genomic sequences
    # for Staphylococcus aureus pulled from NCBI.
    # This sequence IDs in this file should match
    # those in the --input-strain-ids file.
    parser.add_argument(
        '-f', '--input-strain-fasta',
        type = lambda file: check_permissions(parser, file, os.R_OK), 
        required=True,
        help=argparse.SUPPRESS
    )
    # Genomic accessions identifers for the strain
    # of interest, file with single column of IDs
    # to parse from the --input-strain-fasta
    parser.add_argument(
        '-s', '--input-strain-ids',
        type = lambda file: check_permissions(parser, file, os.R_OK), 
        required=True,
        help=argparse.SUPPRESS
    )
    # Output FASTA file to write the strain-specific
    # gene alignment sequences
    parser.add_argument(
        '-o', '--output-strain-fasta',
        type=str, required=True,
        help=argparse.SUPPRESS
    )
    # Option to split the sequence IDs in the input
    # strain FASTA file on the '|' character.
    parser.add_argument(
        '-d', '--split-strain-fasta-ids',
        action='store_true',
        required = False,
        default = False,
        help=argparse.SUPPRESS
    )
    # Get version information
    parser.add_argument(
        '-v', '--version',
        action='version',
        help = argparse.SUPPRESS,
        version='%(prog)s {0}'.format(_VERISON)
    )
    # Add custom help message
    parser.add_argument(
        '-h', '--help',
        action='help',
        help=argparse.SUPPRESS
    )
    return parser.parse_args()


def fasta(filename):
    """
    Reads in a FASTA file and yields each of its entries.
    The generator yields each sequence identifier and its 
    corresponding sequence to ensure a low memory profile. 
    If a sequence occurs over multiple lines, the yielded 
    sequence is concatenated. This function can read both
    plain text and gzip-compressed FASTA files.
    @param filename <str>:
        Path of FASTA file to read and parse
    @yield chrom, sequence <str>, <str>:
        Yields each seq id and seq in the FASTA file
    """
    open_func = gzip.open if filename.endswith('.gz') else open
    with open_func(filename, 'rt') as file:
        sequence, chrom = '', ''
        for line in file:
            line = line.strip()
            if line.startswith('>') and sequence:
                # base case for additional entries
                yield chrom, sequence
                chrom = line[1:] # remove the > symbol
                sequence = ''
            elif line.startswith('>'):
                # base case for first entry in fasta file
                chrom = line[1:] # remove the > symbol
            else:
                # concatenate multi-line sequences
                sequence += line
        else:
            yield chrom, sequence


def fold(sequence, max_length=80):
    """Folds or formats a long sequence so it does not exceed `max_length`.
    If a `sequence` does exceeed the limit, it will wrap sequenced onto a new line.
    @param sequence <str>:
        Coding DNA reference sequence to translate (transcript sequence or CDS sequence)
    @returns folded sequence <str>:
        Folded or formatted sequence that does not exceed `max_length`
    """
    # Clean sequence prior to conversion
    sequence = sequence.strip()
    # Format sequence to max length 
    chunks = []
    for i in range(0, len(sequence)-1, max_length):
        chunk = sequence[i:i+max_length]
        chunks.append(chunk)
    return "\n".join(chunks)


def reverse_complement(sequence):
    """Returns the reverse complement of a DNA sequence.
    @param sequence <str>:
        DNA sequence to find the reverse complement of
    @return rev_comp_sequence <str>:
        Reverse complement of the input DNA sequence
    """
    # Mappings for complementary bases
    # under the IUPAC system:
    # https://en.wikipedia.org/wiki/Nucleic_acid_notation
    COMPLEMENT = {
        # Canonical base pairs
        'A':'T', 'T':'A', 'C':'G', 'G':'C', 
        'a':'t', 't':'a', 'c':'g', 'g':'c',
        'U':'A', 'u':'a', 'N':'N', 'n':'n',
        # Ambiguous base pairs
        'W':'S', 'S':'W', 'M':'K', 'K':'M',
        'w':'s', 's':'w', 'm':'k', 'k':'m',
        'R':'Y', 'Y':'R', 'B':'V', 'V':'B',
        'r':'y', 'y':'r', 'b':'v', 'v':'b',
        'D':'H', 'd':'h', 'H':'D', 'h':'d',
        # Gaps
        '-':'-'
    }
    # Find the reverse complement DNA sequence
    revcomp = ""
    for bp in sequence[::-1]:
        try:
            revcomp += COMPLEMENT[bp]
        except KeyError:
            raise InvalidBasePairError(sequence, bp)
    return revcomp


def index_header(file_header, filename, require=[]):
    """Returns the index of each column_name
    as a dictionary.
    @param file_header <str>:
        First line of a file, containing column names
    @param filename <str>:
        Name of the file the header was parsed from,
        used for error reporting.
    @param require <list[str]>:
        List of column names that must be present
        in the header, otherwise a fatal error
        will be raised.
    @return idx <dict[str]=int>:
        Column name to index dictionary
    """
    idx = {}
    tokens = [
        stripped(c.strip()) \
            for c in file_header.strip().split('\t')
    ]
    # Create column name to index mapping
    for i,c in enumerate(tokens):
        idx[c]=i
    # Check that required columns exist
    missing = [r for r in require if r not in idx ]
    if missing: 
        fatal("Fatal error: Missing required columns in input file '{0}': {1}".format(filename, ",".join(missing)))
    return idx


def split_fasta_id(fasta_id, char='|', field_idx=0):
    """Splits a FASTA sequence ID on the character and field index.
    @param fasta_id <str>:
        FASTA sequence ID to split
    @param char <str>:
        Character to split the FASTA ID on, default '|'
    @param field_idx <int>:
        Field index to return after splitting, default 0
    @return <str>:
        Parsed field of the split FASTA sequence ID.
    """
    return fasta_id.split(char)[field_idx]


def main():
    """
    Main method, entry point that runs when program is directly invoked.
    """
    # Parse command line arguments
    args = parse_cli_arguments()
    log("Running script with the following options: ", args)

    # Create output directory
    output_dir = os.path.abspath(os.path.dirname(args.output_strain_fasta))
    if not os.path.exists(output_dir):
        try: os.makedirs(output_dir)
        except OSError as e: 
            fatal("Fatal error: Failed to create output directory: '{0}'\n{1}".format(output_dir, e))

    # Read in strain sequence IDs
    strain_seqids = set()
    log("Started reading in strain sequence IDs from: ", args.input_strain_ids)
    with open(args.input_strain_ids, 'r') as fh:
        for line in fh:
            line = stripped(line)
            if line:
                strain_seqids.add(line)
    log("Finished reading in strain sequence IDs, parsed {0} strain IDs.".format(len(strain_seqids)))

    # Read in and parse the per-gene filtered
    # BLASTN results file
    blastn_seqids_to_alignment_region = {}
    # Backup default column indices if
    # the BLASTN output file does not
    # contain header with column names
    BLAST_SEQID_IDX = 1          # Column index for subject sequence ID in BLASTN output
    BLAST_SUBJECT_START_IDX = 8  # Column index for subject start in BLASTN output
    BLAST_SUBJECT_END_IDX = 9    # Column index for subject end in BLASTN output
    # Note: if the subject, start is greater than the subject end
    # then the alignment is on the reverse strand. You need to find
    # the reverse complement of the sequence in that case.
    log("Started parsing input BLASTN results file: ", args.input_blast_result)
    blastn_lines_parsed = 0
    with open(args.input_blast_result, 'r') as fh:
        # Try to parse header for column indices,
        # use default indices if header not present
        col2idx = index_header(next(fh), filename = args.input_blast_result)
        BLAST_SEQID_IDX = col2idx.get('sseqid', BLAST_SEQID_IDX)
        BLAST_SUBJECT_START_IDX = col2idx.get('sstart', BLAST_SUBJECT_START_IDX)
        BLAST_SUBJECT_END_IDX = col2idx.get('send', BLAST_SUBJECT_END_IDX)
        for line in fh:
            line = stripped(line)
            if line and not line.startswith('#'):
                fields = line.split('\t')
                subject_seqid = stripped(fields[BLAST_SEQID_IDX])
                subject_start = int(stripped(fields[BLAST_SUBJECT_START_IDX]))
                subject_end = int(stripped(fields[BLAST_SUBJECT_END_IDX]))
                # Store the alignment region for the subject sequence ID
                if subject_seqid in strain_seqids:
                    blastn_lines_parsed += 1
                    blastn_seqids_to_alignment_region[subject_seqid] = (subject_start, subject_end)
    log("Finished parsing input BLASTN results file, saved {0} regions.".format(blastn_lines_parsed))

    # Parse the input staphylococcus aureus FASTA file
    # containing the genomic sequences for all strains,
    # This file was used as input to build the BLASTN DB.
    # It is worth noting this file is massive (346 GB).
    # It contains all the Staphylococcus aureus genomes
    # on NCBI, which is more than 130K genomes.
    log("Starting parsing input Staphylococcus aureus FASTA file: ", args.input_strain_fasta)
    parsed_alignments_from_fasta = 0
    with open(args.output_strain_fasta, 'w') as out_fh:
        for seqid, sequence in fasta(args.input_strain_fasta):
            seqid = stripped(seqid)
            # Split sequence identifer on pipe char,
            # this contains the same accession ID
            # in the blast results
            seqid = seqid if not args.split_strain_fasta_ids \
            else split_fasta_id(
                seqid,
                char = '|',
                field_idx = 0
            )
            if seqid in blastn_seqids_to_alignment_region:
                parsed_alignments_from_fasta += 1
                subject_start, subject_end = blastn_seqids_to_alignment_region[seqid]
                if subject_start < subject_end:
                    # Forward strand alignment,
                    # subtract 1 from start for
                    # 0-based indexing of string
                    alignment_sequence = sequence[subject_start-1:subject_end]
                else:
                    # Reverse strand alignment,
                    # Flip the start and end,
                    # subtract 1 from start for
                    # 0-based indexing of string,
                    # and find the reverse
                    # complement of the sequence
                    alignment_sequence = sequence[subject_end-1:subject_start]
                    # Get reverse complement
                    alignment_sequence = reverse_complement(alignment_sequence)
                # Write alignment sequence to output FASTA file
                out_fh.write(
                    ">{0}\n{1}\n".format(
                        seqid,
                        alignment_sequence
                    )
                )
    log("Finished parsing input Staphylococcus aureus FASTA file, parsed {0} sequences.".format(parsed_alignments_from_fasta))
    log("Finished writing strain-specific gene alignment sequences to output FASTA file: ", args.output_strain_fasta)

if __name__ == '__main__':
    # Call main
    main()