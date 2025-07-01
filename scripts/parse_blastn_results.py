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
            [--length-filter 0.95] \\
            [--percent-identity-filter 95.0] \\
            --input BLASTN_RESULTS \\
            --blastdb-seqs BLASTDB_SEQIDS \\
            --output-prefix OUTPUT_PREFIX

@About:
    Given an output file from blastn, this script
    will parse the results to collapse the results
    on the accession id of the strain sequence,
    and it will report presence/absence of each
    gene and the top hit for each accession id.

@Required:
    -i, --input BLASTN_RESULTS
        Input BLASTN results file to parse.
    -s, --blastdb-seqs BLASTDB_SEQIDS
        File containing seqeunce IDs of all
        sequences in the blastn databases.
        This can be parsed from the FASTA
        file that was used to create the
        blastn database via the following
        command below:
          $ grep '^>' file.fa \\
              | sed 's/^>//g' > blastn_seqids.txt
    -o, --output-prefix OUTPUT_PREFIX
        Output prefix to writing the collapsed
        results. This script will produce more than
        one file as output.
@Options:
    -l, --length-filter 0.95
        Filters any blastn results whose length is
        less than this value as a function to the
        length of the gene of interest. A value of
        1.00 would represent a blastn result whose
        length is same as the gene of interest.
        Range: 0.0-1.0
        Default: 0.95
    -p, --percent-identity-filter 95.0
        Filters any blastn results whose percent
        identity score is less than this value. A
        value of 100 would represent a perfect match.
        Range: 0.0-100.0
        Default: 95.0
    -h, --help     Shows this help message and exits.
    -v, --version  Prints the version and exits.

@Example:
    $ ./{0} \\
            --input blastn_results-with-header.tsv \\
            --blastdb-seqs blastn_seqids.txt
            --output blastn_filtered
""".format(_NAME)
)

# Semantic version
_VERISON = '1.0.0'


# Helper functions
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
    # Input BLASTN file to parse
    parser.add_argument(
        '-i', '--input',
        type=str, required=True,
        help=argparse.SUPPRESS
    )
    # BlastDB sequence identifers 
    parser.add_argument(
        '-s', '--blastdb-seqs',
        type=str, required=True,
        help=argparse.SUPPRESS
    )
    # Output prefix for results
    parser.add_argument(
        '-o', '--output-prefix',
        type=str, required=True,
        help=argparse.SUPPRESS
    )
    # Length filter
    parser.add_argument(
        '-l', '--length-filter',
        type=float, required=False,
        default=0.95,
        help=argparse.SUPPRESS
    )
    # Percent identity filter
    parser.add_argument(
        '-p', '--percent-identity-filter',
        type=float, required=False,
        default=95.00,
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


def stripped(s):
    """Cleans string to remove quotes
    @param s <str>:
        String to remove quotes or clean
    @return s <str>:
        Cleaned string with quotes removed
    """
    return s.strip('"').strip("'")


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


def index_header(file_header):
    """Returns the index of each column_name
    as a dictionary.
    @param file_header <str>:
        First line of a file, containing column names
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
    return idx


def get_with_default(line_list, column_name_idx_dict, column_name, default_value="NA"):
    """Get a value from a list using the column name index
    dictionary. If the column name does not exist in the
    dictionary, return the default value (i.e "NA").
    @param line_list <list[str]>:
        List of values from a line in a file. This is the
        list that get are retrieving information from.
        This function is used to return a value (with a
        default value if missing) from within a list or
        dictionary comprehension.
    @param column_name_idx_dict <dict[str]=int>:
        Dictionary mapping column names to their index
    @param column_name <str>:
        Column name to look up in the dictionary
    @param default_value <str>:
        Default value to return if the column name does not
        exist in the dictionary. Defaults to "NA".
    @return value <str>:
        Value from the list at the index of the column name,
        or the default value if the column name does not exist.
    """
    # Default value to return if column name DNE
    parsed_value = default_value
    # Try to get the index of the column name
    # This can result in a KeyError/IndexError
    # if the column name does not exist in the
    # dict. 
    if column_name in column_name_idx_dict:
        # Get the index of the column name
        # from the dictionary.
        # This will raise a KeyError if the
        # column name does not exist in the dict.
        list_idx = column_name_idx_dict[column_name]
        if list_idx < len(line_list):
            # If the index is within the bounds of the list,
            # return the value at that index.
            parsed_value = line_list[list_idx]
            # Remove quotes from the value
            parsed_value = stripped(parsed_value)
    return parsed_value


def index_file(
        file,
        first_key="sseqid",
        second_key="qseqid",
        multi_value_columns=["pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore"]
    ):
    """Parses and indexes a blastn results file into a dictionary
    for quick lookups later. This file is indexed by keys representing
    the accession_id and gene_id, and it contains value representing a
    list of dictionary where each dictionary contains blastn results
    for a contig (1:M relationship between accession_id and different
    contig results):
        • pident
        • length
        • mismatch
        • gapopen
        • qstart
        • qend
        • sstart
        • send
        • evalue
        • bitscore
    i.e indexed_file[sseqid][qseqid] = [{...},{...}, ...]
    @param file <str>:
        File to parse and index. Must contain a header with
        the columns listed in keys and values. The index of
        these columns will be automatically resolved.
    @param first_key <str>:
        Key to use for the first level of the index. Defaults
        to "sseqid". This is the first key in the nested 
        dictionary. Please note: the first key will be
        split '|' where the first element is parsed/used.
    @param second_key <str>:
        Key to use for the second level of the index. Defaults
        to "qseqid". This is the second key in the
        nested dictionary.
    @param multi_value_columns <list[str]>:
        List of columns that contain multiple values per key.
        These values are stored as a list of dicts in the
        nested dictionary under [sseqid][qseqid].
        This represents the blastn results for a given
        accession_id for a given gene. Each element in the
        returned list is a dict containing the following
        information:
            • pident
            • length
            • mismatch
            • gapopen
            • qstart
            • qend
            • sstart
            • send
            • evalue
            • bitscore
    @return file_idx <dict[str][str]=list[dict,dict,...]>:
        Nested dictionary where,
            • 1st_key = sseqid
            • 2nd_key = qseqid
    @note: file_idx[sseqid][qseqid] returns a list of dicts,
    where each element in the list contain blastn results.
    """
    log("Started indexing input file: {0}".format(file))
    file_idx = {}
    # Handler for opening files, i.e.
    # uncompressed or gzip files
    open_func = gzip.open if file.endswith('.gz') else open
    line_number = 0  # Used for error reporting
    with open_func(file, 'rt') as fh:
        header = next(fh)
        col_idx = index_header(header)
        for line in fh:
            # Increment line number
            line_number += 1
            # Split the line into columns
            tokens = line.strip().split('\t')
            # Concatenate mutiple keys into
            # a single key separated by the
            # key_delim character
            _k1 = tokens[col_idx[first_key]].split("|")[0]   # represents accession_id
            _ds  = tokens[col_idx[first_key]].split("|")[1]  # description
            _k2 = tokens[col_idx[second_key]]                # reprsents gene_id
            if _k1 not in file_idx:
                file_idx[_k1] = {}
            if _k2 not in file_idx[_k1]:
                file_idx[_k1][_k2] = []
            # Parse multi-value columns
            mv_attrs = {v: get_with_default(tokens,col_idx,v) for v in multi_value_columns}
            mv_attrs["description"] = _ds
            file_idx[_k1][_k2].append(mv_attrs)
    log("Finished indexing input file: {0} ({1} lines)".format(file, line_number))
    return file_idx


if __name__ == '__main__':
    # Parse command line arguments
    args = parse_cli_arguments()

    # Sanity check for usage
    if len(sys.argv) == 1:
        # Nothing was provided
        fatal('Invalid usage: {0} [-h] ...'.format(os.path.basename(sys.argv[0])))

    log("Running parse blastn script with the following options: ", args)
    # Create output directory if
    # it does not exist
    output_dir = os.path.abspath(os.path.dirname(args.output_prefix))
    if not os.path.exists(output_dir):
        try: os.makedirs(output_dir)
        except OSError as e:
            fatal(
                "Fatal error: Failed to create output directory: {0}\n{1}".format(
                    output_dir, e
                )
            )

    # Gene lengths
    GENE_LENGTHS = {
        "BlaZ": 846,
        "MecA": 2007,
        "MecI": 372,
        "Stk1": 1995,
        "Stp1": 744,
        "BlaR1": 1758,
        "BlaI": 381
    }

    # Parse and collapse blastn results
    FIRST_KEY = "sseqid"
    SECOND_KEY = "qseqid"
    PARSE_1toM_COLUMNS = [
        "pident",
        "length",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "sstart",
        "send",
        "evalue",
        "bitscore"
    ]
    blastn_dict = index_file(
        file=args.input,
        first_key = FIRST_KEY,
        second_key=SECOND_KEY,
        multi_value_columns=PARSE_1toM_COLUMNS,
    )
    log('Started filtering blastn results based on percent identity and length filters')
    # Filter results based on percent identity
    # and length filters
    for k,v in blastn_dict.items():
        for gene, results_list in blastn_dict[k].items():
            # Get the gene length
            gene_length = GENE_LENGTHS.get(gene, 0)
            # Filter results based on percent identity and length
            blastn_dict[k][gene] = [
                result for result in results_list
                if (float(result["pident"]) >= args.percent_identity_filter and
                    float(result["length"]) >= (args.length_filter * gene_length))
            ]
    log('Finished filtering blastn results based on percent identity and length filters')

    log("Started collapsing blastn results to a single top hit per gene per accession_id")
    # Filter 1:M blastn results to
    # get the top hit for genes 
    # aligning to different scaffolds,
    # sort results as a function of
    # length and percent identity
    # and collapse results to a
    # single result.
    for k,v in blastn_dict.items():
        for gene, results_list in blastn_dict[k].items():
            # Sort results by length and percent identity
            results_list.sort(
                key=lambda x: (float(x["length"]), float(x["pident"])),
                reverse=True
            )
            # Collapse results to a single result
            if results_list:
                blastn_dict[k][gene] = results_list[0]  # dict
            else:
                blastn_dict[k][gene] = []
    log("Finished collapsing blastn results to a single top hit per gene per accession_id")
    
    # Write collapsed results to output file
    log("Started writing collapsed results to file: {0}_collapsed_results.tsv".format(args.output_prefix))
    WRITE_1toM_COLUMNS = [
        "pident",
        "length",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "sstart",
        "send",
        "evalue",
        "bitscore",
        "description"
    ]
    with open(args.output_prefix + "_collapsed_results.tsv", 'w') as out_fh:
        # Write header
        header = [
            SECOND_KEY, FIRST_KEY,
            *WRITE_1toM_COLUMNS
        ]
        out_fh.write("\t".join(header) + "\n")
        # Write results
        for k,v in blastn_dict.items():
            for gene, result in v.items():
                if isinstance(result, list) and not result:
                    continue  # Skip empty results
                # Flatten the result dict
                result_values = [result.get(col, "NA") for col in WRITE_1toM_COLUMNS]
                out_fh.write("\t".join([gene, k] + result_values) + "\n")
    log("Finished writing collapsed results to file: {0}_collapsed_results.tsv".format(args.output_prefix))
    
    log("Started in blastn sequence identifter file")
    blastdb_sequences = []
    with open(args.blastdb_seqs, 'r') as file:
        # Read sequence identifiers
        blastdb_sequences = [stripped(line.strip()) for line in file if line.strip()]
    log("Finished in blastn sequence identifter file")

    log("Started creating gene presence/absence output file: {0}_gene_presence_absence.tsv".format(args.output_prefix))    
    # Create gene presence/absence output file
    # ACCESSION_ID  BlaZ    MecA    MecI   Stk1   Stp1   BlaR1   BlaI   NGenes
    # ID1           0       1       0      1      1      0       1      4
    # ID2           1       0       1      0      1      1       1      5
    with open(args.output_prefix + "_gene_presence_absence.tsv", 'w') as out_fh:
        # Write header
        header = ["ACCESSION_ID"] + list(GENE_LENGTHS.keys()) + ["NGenes"]
        out_fh.write("\t".join(header) + "\n")
        # Write results
        for seq_id in blastdb_sequences:
            # Initialize presence/absence dictionary
            presence_absence = {gene: 0 for gene in GENE_LENGTHS.keys()}
            # Check if the sequence ID exists in the blastn_dict
            if seq_id in blastn_dict:
                for gene, result in blastn_dict[seq_id].items():
                    if isinstance(result, dict) and result:
                        # If the result is a non-empty dict,
                        # then the gene is present
                        presence_absence[gene] = 1
            # Count number of genes present
            n_genes_present = sum(presence_absence.values())
            # Write presence/absence results
            presence_absence_values = [str(presence_absence[gene]) for gene in GENE_LENGTHS.keys()]
            out_fh.write("\t".join([seq_id] + presence_absence_values + [str(n_genes_present)]) + "\n")
    log("Finished creating gene presence/absence output file: {0}_gene_presence_absence.tsv".format(args.output_prefix))
