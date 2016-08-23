#!/usr/bin/env python3

import sys
if sys.version_info.major is not 3:
    sys.exit("Please use Python 3 for this module: " + __name__)

import os

from Objects import snp

try:
    from overload import overload
    from Bio.Blast.Applications import NcbiblastnCommandline
except ImportError as error:
    sys.exit("Failed to find  " + error.name)


#   A custom error
class BLASTFailedError(Exception):
    """BLAST seems to have failed..."""


#   A function to run BLASTn
@overload
def run_blastn(query, subject, evalue, max_seqs, max_hsps, keep_query):
    """Run BLASTn"""
    try:
        assert isinstance(query, str)
        assert isinstance(subject, str)
        assert isinstance(evalue, float)
        assert isinstance(max_seqs, int)
        assert isinstance(max_hsps, int)
        assert isinstance(keep_query, bool)
    except AssertionError:
        raise
    #   Create an output name
    query_base = os.path.basename(os.path.splitext(query)[0])
    subject_base = os.path.basename(os.path.splitext(subject)[0])
    blast_out = os.getcwd() + '/' + query_base + '_' + subject_base + '_BLAST.xml'
    #   Setup BLASTn
    blastn = NcbiblastnCommandline(
        query=query,
        subject=subject,
        evalue=evalue,
        outfmt=5,
        max_target_seqs=max_seqs,
        max_hsps=max_hsps,
        out=blast_out
    )
    #   Run BLASTn
    print(blastn, file=sys.stderr)
    blastn()
    if not os.path.exists(blast_out):
        raise BLASTFailedError
    if not keep_query:
        os.remove(query)
    return blast_out

@run_blastn.add
def run_blastn(bconf):
    """Run BLASTn"""
    try:
        assert isinstance(bconf, dict)
        query = bconf['query']
        subject = bconf['subject']
        evalue = bconf['evalue']
        max_seqs = bconf['max_seqs']
        max_hsps = bconf['max_hsps']
        keep_query = bconf['keep_query']
        blast_out = run_blastn(
            query=query,
            subject=subject,
            evalue=evalue,
            max_seqs=max_seqs,
            max_hsps=max_hsps,
            keep_query=keep_query
        )
        return blast_out
    except AssertionError:
        raise
    except KeyError as error:
        sys.exit("Malformed BLAST config! Missing entry: " + str(error.args))


#   A function to write a FASTA file from a lookup dictionary
def write_fasta(lookup_dict, query):
    """Write a FASTA file from the lookups in a dictionary"""
    try:
        assert isinstance(lookup_dict, dict)
        for lookup in lookup_dict.values():
            assert isinstance(lookup, snp.Lookup)
    except AssertionError:
        raise
    query_base = os.path.basename(os.path.splitext(query)[0])
    fasta_file = os.getcwd() + '/' + query_base + '.fasta'
    with open(fasta_file, 'w') as f:
        for lookup in lookup_dict.values():
            f.write(lookup.format_fasta())
            f.write('\n')
    return fasta_file
