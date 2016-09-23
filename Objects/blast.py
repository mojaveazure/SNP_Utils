#!/usr/bin/env python3
"""This module holds two classes and one error"""

import sys
if sys.version_info.major is not 3:
    sys.exit("Please use Python 3 for this module: " + __name__)

from itertools import repeat

from . import snp

try:
    from overload import overload
    from bs4 import element
except ImportError as error:
    sys.exit("Please install " + error.name)


#   An error I probably overuse...
class NoSNPError(Exception):
    """A SNP has not been found"""


#   A class definition for a BLAST Hsp
class Hsp(object):
    """This is a class for a BLAST Hsp
    It continas the following information:
        Chromosome name
        Query name
        Hsp e-value
        Query sequence
        Subject sequence
        Hsp start relative to subject
        Hsp end relative to subject
        Subject strand (forward or reverse)
        SNP Position relative to subject (once calcualted with Hsp.get_snp_position())
        """
    def __init__(self, chrom, name, evalue, qseq, hseq, hstart, hend, hstrand):
        try:
            assert isinstance(chrom, str)
            assert isinstance(name, str)
            assert isinstance(evalue, float)
            assert isinstance(qseq, str)
            assert isinstance(hseq, str)
            assert isinstance(hstart, int)
            assert isinstance(hend, int)
            assert isinstance(hstrand, int)
            assert hstrand == 1 or hstrand == -1
        except AssertionError:
            raise
        self._chrom = chrom
        self._name = name
        self._evalue = evalue
        self._query = qseq
        self._subject = hseq
        self._start = hstart
        self._end = hend
        self._hstrand = hstrand
        self._snp_pos = None

    def __repr__(self):
        return self._name + ":" + str(self._evalue)

    def get_chrom(self):
        """Get the chromosome that the hsp matched to"""
        return self._chrom

    def get_name(self):
        """Get the query name of the hsp"""
        return self._name

    def get_rc(self):
        """Did the query align to the forward (False) or reverse (True) strand"""
        return self._hstrand == -1

    def get_snp_position(self, query_snp, expected):
        """Calculate the SNP position relative to reference"""
        try:
            assert len(query_snp) is 1
            assert isinstance(expected, int)
        except AssertionError:
            raise
        #   If we don't find the SNP in this hsp
        if self._query.find(query_snp) is -1:
            raise NoSNPError # Error out
        #   If the hsp sequence is less than our expected value
        elif len(self._query) < expected or (self._query.find(query_snp) < expected and self._query.count(query_snp) == 1):
            #   Do some guesswork as to where the SNP is
            self._snp_pos = self._query.find(query_snp)
            if self.get_rc():
                return self._start - self._snp_pos
            else:
                return self._start + self._snp_pos
        #   Otherwise
        else:
            #   Start with the expected position
            self._snp_pos = expected
            #   Get a substring up to this position
            qsubstr = self._query[:self._snp_pos + 1]
            #   While we haven't found our query SNP
            while qsubstr[-1] != query_snp:
                old = self._snp_pos # Old position
                self._snp_pos = expected + qsubstr.count('-') # Add on any indels to our position
                if old is self._snp_pos: # If we didn't change our position guess
                    raise NoSNPError # Error out
                qsubstr = self._query[:self._snp_pos + 1] # Replace our substring
            #   Calculate insertions
            hsubstr = self._subject[:self._snp_pos + 1]
            insertion = hsubstr.count('-')
            if self.get_rc():
                return self._start - self._snp_pos - insertion
            else:
                return self._start + self._snp_pos - insertion

    @overload
    def get_subject_allele(self):
        """Get the reference allele if the SNP has been found"""
        try:
            assert self._snp_pos is not None
            ref = self._subject[self._snp_pos]
            assert ref.upper() in 'ACGTN'
            return ref
        except AssertionError:
            raise NoSNPError
        # if self._snp_pos is None:
        #     raise NoSNPError
        # return self._subject[self._snp_pos]

    @get_subject_allele.add
    def get_subject_allele(self, query_snp, expected):
        """Get the reference allele if the SNP has been found"""
        genomic_position = self.get_snp_position(query_snp, expected)
        return (genomic_position, self.get_subject_allele())


#   A class definition for holding Iterations
class SNPIteration(object):
    """This is a class for holding BLAST iterations
    It holds the SNP name and a list of hsps"""
    @staticmethod
    def GET_VALUE(tag, value):
        try:
            assert isinstance(tag, element.Tag)
        except AssertionError:
            raise
        return tag.findChild(value).text

    _VALS = ['Hsp_evalue', 'Hsp_hit-from', 'Hsp_hit-to', 'Hsp_hit-frame', 'Hsp_qseq', 'Hsp_hseq']
    def __init__(self, iteration):
        try:
            assert isinstance(iteration, element.Tag)
        except AssertionError:
            raise
        self._snpid = SNPIteration.GET_VALUE(iteration, 'Iteration_query-def')
        self._hsps = []
        self._fail = False
        #   Start parsing hits and hsps
        for hit in iteration.findAll('Hit'):
            self._parse_hit(hit)
        # If we don't have any hsps, set 'self._fail' to True
        if len(self._hsps) < 1:
            self._fail = True

    def __repr__(self):
        return self._snpid + '(' + str(len(self._hsps)) + ' hsp(s))'

    def _parse_hsp(self, hsp):
        """Parse the hsp section of a BLAST XML"""
        try:
            assert isinstance(hsp, element.Tag)
        except AssertionError:
            raise
        #   If there are no gaps in this hsp
        if hsp.findChild('Hsp_midline').text.count(' ') < 1:
            #   Skip
            raise NoSNPError
        #   Collect all values
        hsp_vals = map(SNPIteration.GET_VALUE, repeat(hsp, len(self._VALS)), self._VALS)
        #   Return as a tuple
        return tuple(hsp_vals)

    def _parse_hit(self, hit):
        """Parse the hit section of a BLAST XML"""
        try:
            assert isinstance(hit, element.Tag)
        except AssertionError:
            raise
        chrom = SNPIteration.GET_VALUE(hit, 'Hit_def') # Get the chromosome information
        hsps = [] # A list to hold hsps
        for hsp in hit.findAll('Hsp'):
            try:
                #   Try to parse the hsp
                hsp_vals = self._parse_hsp(hsp)
                hsps.append(hsp_vals)
            except NoSNPError: # If there's no gap, skip
                continue
        for hsp in hsps:
            #   Unpack our tuple
            (evalue, hsp_start, hsp_end, strand, query, reference) = hsp
            #   Make a Hit
            hsp = Hsp(
                chrom=chrom,
                name=self._snpid,
                evalue=float(evalue),
                qseq=query,
                hseq=reference,
                hstart=int(hsp_start),
                hend=int(hsp_end),
                hstrand=int(strand)
            )
            #   Add our hsp to the list of hsp
            self._hsps.append(hsp)

    def get_snpid(self):
        """Get the SNP ID for this iteration"""
        return self._snpid

    def check_fail(self):
        """See if we are lacking any hits"""
        if self._fail:
            raise NoSNPError

    def hit_snps(self, lookup):
        """Find SNPs for our Hits"""
        try:
            assert isinstance(lookup, snp.Lookup)
            self.check_fail()
        except AssertionError:
            raise
        except NoSNPError:
            raise
        #   Some holding lists
        no_snp = []
        snp_list = []
        try:
            for hit in self._hsps:
                #   Make a hit out of every SNP
                s = snp.SNP(lookup, hit)
                snp_list.append(s)
                # snp_list = [snp.SNP(lookup, hit) for hit in self._hits]
        except NoSNPError:
            no_snp.append(lookup.get_snpid())
        return(snp_list, no_snp)


