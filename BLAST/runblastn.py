#!/usr/bin/env python3

import sys
if sys.version_info.major is not 3:
    sys.exit("Please use Python 3 for this module: " + __name__)

import os
import re

from Objects import snp

try:
    from overload import overload
    from Bio.Blast.Applications import NcbiblastnCommandline
except ImportError as error:
    sys.exit("Failed to find  " + error.name)


#   A custom error
class BLASTFailedError(Exception):
    """BLAST seems to have failed..."""


#   Another custom error
class MalformedConfigError(Exception):
    """The BLAST config file isn't formed correctly"""


#   Ensure we have our BLAST database
def validate_db(db_path):
    """Find a BLAST nucleotide database"""
    try:
        assert isinstance(db_path, str)
    except AssertionError:
        raise TypeError
    db_name = os.path.basename(db_path)
    nhr = re.compile(r'(%s\.[0-9\.]*nhr)' % db_name).search
    nin = re.compile(r'(%s\.[0-9\.]*nin)' % db_name).search
    nsq = re.compile(r'(%s\.[0-9\.]*nsq)' % db_name).search
    nal = re.compile(r'(%s\.*nal)' % db_name).search
    db_directory = os.path.abspath(os.path.realpath(os.path.dirname(db_path)))
    if not db_directory:
        db_directory = os.getcwd()
    print('Searching for proper database files for', db_name, 'in', db_directory, file=sys.stderr)
    db_contents = '\n'.join(os.listdir(db_directory))
    if not nhr(db_contents) and not nin(db_contents) and not nsq(db_contents) and not nal(db_contents):
        raise FileNotFoundError("Failed to find the BLAST nucleotide database at " + db_path)


#   A function to run BLASTn
@overload
def run_blastn(cline, keep_query=True):
    """Run BLASTn"""
    try:
        assert isinstance(cline, NcbiblastnCommandline)
        assert isinstance(keep_query, bool)
    except AssertionError:
        raise TypeError
    print(cline, file=sys.stderr)
    cline()
    if not os.path.exists(cline.out):
        raise BLASTFailedError
    if not keep_query:
        os.remove(cline.query)
    return cline.out


@run_blastn.add
def run_blastn(query, subject, evalue, max_hits, max_hsps, identity, keep_query):
    """Run BLASTn"""
    try:
        assert isinstance(query, str)
        assert isinstance(subject, str)
        assert isinstance(evalue, float)
        assert isinstance(max_hits, int)
        assert isinstance(max_hsps, int)
        assert isinstance(identity, float)
        assert isinstance(keep_query, bool)
    except AssertionError:
        raise TypeError
    #   Create an output name
    print("Running BLAST against subject:", subject, file=sys.stderr)
    query_base = os.path.basename(os.path.splitext(query)[0])
    db_base = os.path.basename(os.path.splitext(subject)[0])
    blast_out = os.getcwd() + '/' + query_base + '_' + db_base + '_BLAST.xml'
    #   Setup BLASTn
    blastn = NcbiblastnCommandline(
        query=query,
        subject=subject,
        evalue=evalue,
        outfmt=5,
        max_target_seqs=max_hits,
        max_hsps=max_hsps,
        perc_identity=identity,
        out=blast_out
    )
    #   Run BLASTn
    outfile = run_blastn(cline=blastn, keep_query=keep_query)
    return outfile


@run_blastn.add
def run_blastn(query, database, evalue, max_hits, max_hsps, identity, keep_query):
    """Run BLASTn"""
    try:
        assert isinstance(query, str)
        assert isinstance(database, str)
        assert isinstance(evalue, float)
        assert isinstance(max_hits, int)
        assert isinstance(max_hsps, int)
        assert isinstance(identity, float)
        assert isinstance(keep_query, bool)
    except AssertionError:
        raise TypeError
    #   Create an output name
    print("Running BLAST with database:", database, file=sys.stderr)
    query_base = os.path.basename(os.path.splitext(query)[0])
    db_base = os.path.basename(os.path.splitext(database)[0])
    blast_out = os.getcwd() + '/' + query_base + '_' + db_base + '_BLAST.xml'
    #   Setup BLASTn
    blastn = NcbiblastnCommandline(
        query=query,
        db=database,
        evalue=evalue,
        outfmt=5,
        max_target_seqs=max_hits,
        max_hsps=max_hsps,
        perc_identity=identity,
        out=blast_out
    )
    #   Run BLASTn
    outfile = run_blastn(cline=blastn, keep_query=keep_query)
    return outfile


@run_blastn.add
def run_blastn(bconf):
    """Run BLASTn"""
    try:
        assert isinstance(bconf, dict)
        query = bconf['query']
        evalue = bconf['evalue']
        max_hits = bconf['max_hits']
        max_hsps = bconf['max_hsps']
        identity = bconf['identity']
        keep_query = bconf['keep_query']
    except AssertionError:
        raise TypeError
    except KeyError:
        raise MalformedConfigError
    try:
        if 'database' in bconf.keys():
            database = bconf['database']
            validate_db(database)
            blast_out = run_blastn(
                query=query,
                database=database,
                evalue=evalue,
                max_hits=max_hits,
                max_hsps=max_hsps,
                identity=identity,
                keep_query=keep_query
            )
        elif 'subject' in bconf.keys():
            subject = bconf['subject']
            blast_out = run_blastn(
                query=query,
                subject=subject,
                evalue=evalue,
                max_hits=max_hits,
                max_hsps=max_hsps,
                identity=identity,
                keep_query=keep_query
            )
        else:
            raise MalformedConfigError
        return blast_out
    except (FileNotFoundError, MalformedConfigError):
        raise


#   A function to write a FASTA file from a lookup dictionary
def write_fasta(lookup_dict, query):
    """Write a FASTA file from the lookups in a dictionary"""
    try:
        assert isinstance(lookup_dict, dict)
        for lookup in lookup_dict.values():
            assert isinstance(lookup, snp.Lookup)
    except AssertionError:
        raise TypeError
    query_base = os.path.basename(os.path.splitext(query)[0])
    fasta_file = os.getcwd() + '/' + query_base + '.fasta'
    with open(fasta_file, 'w') as f:
        for lookup in lookup_dict.values():
            f.write(lookup.format_fasta())
            f.write('\n')
    return fasta_file
