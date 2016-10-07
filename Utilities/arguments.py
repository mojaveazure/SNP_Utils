#!/usr/bin/env python3

import sys
if sys.version_info.major is not 3:
    sys.exit("Please use Python 3 for this module")

import argparse
import os


BLAST_CONFIG = os.getcwd() + '/blast_config.ini'

#   Make an argument parser
def make_argument_parser():
    parser = argparse.ArgumentParser(add_help=True)
    #   Subparser for BLAST- vs Alignment-based SNP prediction
    subparsers = parser.add_subparsers(
        title='Subroutine',
        description='Choose a subroutine',
        dest='method',
        help="'BLAST' means we are running a BLAST search to find SNPs, 'SAM' means we are using a SAM file, 'CONFIG' is used to configure a BLAST search"
    )
    #   BLAST subparser
    blast = subparsers.add_parser('BLAST')
    blast_modes = blast.add_mutually_exclusive_group(required=True)
    blast_modes.add_argument(
        '-x',
        '--xml',
        dest='xml',
        type=argparse.FileType('r'),
        default=None,
        required=False,
        metavar='BLAST XML RESULTS',
        help="Optional BLAST XML results already found, incompatible with '-c | --config'"
    )
    blast_modes.add_argument(
        '-c',
        '--config',
        dest='config',
        type=str,
        default=None,
        metavar='BLAST CONFIG FILE',
        help="Configuration file for running BLAST, incompatible with '-x | --xml'"
    )
    blast.add_argument(
        '-l',
        '--lookup',
        dest='lookup',
        type=str,
        default=None,
        required=True,
        metavar='ILLUMINA LOOKUP TABLE',
        help="SNP lookup table in Illumina format"
    )
    blast.add_argument(
        '-o',
        '--outname',
        dest='outname',
        type=str,
        default='output',
        required=False,
        metavar='OUTPUT NAME',
        help="Name of output file, without suffix. Defaults to 'output'"
    )
    blast.add_argument(
        '-r',
        '--rank-snps',
        dest='rank',
        action='store_const',
        const=True,
        default=False,
        required=False,
        help="Do we rank the SNPs and take the lowest e-value/highest scoring BLAST results? Pass '-r | --rank-snps' to say yes"
    )
    #   Alignment subparser
    sam = subparsers.add_parser('SAM')
    sam.add_argument(
        '-s',
        '--sam-file',
        dest='samfile',
        type=argparse.FileType('r'),
        default=None,
        required=True,
        metavar='SAM FILE',
        help="SAM File with aligned SNP contextual sequences"
    )
    sam.add_argument(
        '-r',
        '--reference',
        dest='reference',
        type=str,
        default=None,
        required=True,
        metavar='REFERENCE SEQUENCE',
        help="Reference sequence in FASTA format"
    )
    sam.add_argument(
        '-l',
        '--lookup',
        dest='lookup',
        type=str,
        default=None,
        required=True,
        metavar='ILLUMINA LOOKUP TABLE',
        help="SNP lookup table in Illumina format"
    )
    sam.add_argument(
        '-o',
        '--outname',
        dest='outname',
        type=str,
        default='output',
        required=False,
        metavar='OUTPUT NAME',
        help="Name of output file, without suffix. Defaults to 'output'"
    )
    #   Config subparser
    config = subparsers.add_parser('CONFIG')
    ref_db = config.add_mutually_exclusive_group(required=True)
    ref_db.add_argument(
        '-r',
        '--reference',
        dest='subject',
        type=str,
        default=None,
        metavar='REFERNCE SEQUENCE',
        help="Reference sequence in FASTA format, incompatible with '-d | --database'"
    )
    ref_db.add_argument(
        '-d',
        '--database',
        dest='database',
        type=str,
        default=None,
        metavar='REFERENCE DATABASE',
        help="Reference BLAST nucleotide database, incompatible with '-r | --reference'"
    )
    config.add_argument(
        '-e',
        '--evalue',
        dest='evalue',
        type=float,
        default=1e-1,
        required=False,
        metavar='E-VALUE',
        help="E-value threshold, defaults to '1e-1'"
    )
    config.add_argument(
        '-s',
        '--max-seqs',
        dest='max_seqs',
        type=int,
        default=3,
        required=False,
        metavar='MAX SEQS',
        help="Max seqs to keep per query, defaults to '3'"
    )
    config.add_argument(
        '-m',
        '--max-hsps',
        dest='max_hsps',
        type=int,
        default=3,
        required=False,
        metavar='MAX HSPS',
        help="Max hsps to keep per hit, defaults to '3'"
    )
    config.add_argument(
        '-i',
        '--identity',
        dest='identity',
        type=float,
        default=99,
        required=False,
        metavar='PERCENT IDENTITY',
        help="Percent identity to match, defaults to '99'"
    )
    config.add_argument(
        '-k',
        '--keep-query',
        dest='keep_query',
        action='store_const',
        const=True,
        default=False,
        required=False,
        help="Do we keep the query FASTA file? Pass '-k | --keep-query' to say yes"
    )
    config.add_argument(
        '-c',
        '--config',
        dest='outfile',
        default=BLAST_CONFIG,
        required=False,
        metavar='BLAST CONFIG FILE',
        help="Name of BLAST config file to write to, defaults to '" + BLAST_CONFIG + "'"
    )
    return parser
