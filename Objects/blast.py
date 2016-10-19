#!/usr/bin/env python3
"""This module holds two classes and one error"""

import sys
if sys.version_info.major is not 3:
    sys.exit("Please use Python 3 for this module: " + __name__)


from itertools import repeat

from . import snp
from Utilities.utilities import rank_remove

try:
    from overload import overload
    from bs4 import element
except ImportError as error:
    sys.exit("Please install " + error.name)


#   A tuple of values that need to be parsed for each HSP
VALS = (
    'Hsp_bit-score',
    'Hsp_evalue',
    'Hsp_hit-from',
    'Hsp_hit-to',
    'Hsp_hit-frame',
    'Hsp_identity',
    'Hsp_align-len',
    'Hsp_qseq',
    'Hsp_hseq'
)

#   The message found in a BLAST iteration if no hits were found
NO_HIT_MESSAGE = 'No hits found'

#   An error I probably overuse...
class NoSNPError(Exception):
    """A SNP has not been found"""


#   A class definition for a BLAST Hsp
class Hsp(object):
    """This is a class for a BLAST Hsp
    It continas the following information:
        Chromosome name                     (str)
        Query name                          (str)
        Hsp e-value                         (float)
        Query sequence                      (str)
        Subject sequence                    (str)
        Hsp start relative to subject       (int)
        Hsp end relative to subject         (int)
        Subject strand (forward or reverse) [1 | -1]
        Hsp bit-score                       (int or float)
        Identity                            (int)
        Aligned length                      (int)
        SNP Position relative to subject (once calcualted with Hsp.get_snp_position())
        """
    def __init__(self, chrom, name, evalue, qseq, hseq, hstart, hend, hstrand, bit_score, identity, aligned_length):
        try:
            assert isinstance(chrom, str)
            assert isinstance(name, str)
            assert isinstance(evalue, float)
            assert isinstance(qseq, str)
            assert isinstance(hseq, str)
            assert isinstance(hstart, int)
            assert isinstance(hend, int)
            assert isinstance(hstrand, int)
            assert isinstance(bit_score, int) or isinstance(bit_score, float)
            assert isinstance(identity, int)
            assert isinstance(aligned_length, int)
        except AssertionError:
            raise TypeError
        try:
            assert hstrand == 1 or hstrand == -1
            assert identity <= aligned_length
        except AssertionError:
            raise ValueError
        self._chrom = chrom
        self._name = name
        self._evalue = evalue
        self._query = qseq
        self._subject = hseq
        self._start = hstart
        self._end = hend
        self._hstrand = hstrand
        self._bits = bit_score
        self._identity = identity
        self._alength = aligned_length
        self._snp_pos = None
        self._snp = None

    def __repr__(self):
        return self._name + ":" + str(self._evalue)

    def __eq__(self, other):
        if isinstance(other, Hsp):
            name_bool = self._name == other._name
            eval_bool = self._evalue == other._evalue
            score_bool = self._bits == other._bits
            return name_bool and eval_bool and score_bool
        elif isinstance(other, str):
            return self._name == other
        elif isinstance(other, float):
            return self._evalue == other
        else:
            return NotImplemented

    def __lt__(self, other):
        if isinstance(other, Hsp):
            if self._evalue == other._evalue:
                #   A small evalue should go with a larger score
                return self._bits > other._bits
            else:
                return self._evalue < other._evalue
        elif isinstance(other, float):
            return self._evalue < other
        else:
            return NotImplemented

    def __le__(self, other):
        if isinstance(other, Hsp):
            if self._evalue <= other._evalue:
                #   A small evalue should go with a larger score
                return self._bits >= other._bits
            else:
                return False
        elif isinstance(other, float):
            return self._evalue <= other
        else:
            return NotImplemented

    def __bool__(self):
        return bool(self._snp_pos)

    def __hash__(self):
        return hash(self._name)

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
        try: # Type checking
            assert isinstance(expected, int)
        except AssertionError:
            raise TypeError
        try: # Value checking
            assert len(query_snp) is 1
        except AssertionError:
            raise ValueError
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
        """Get the reference allele"""
        try:
            assert self._snp_pos is not None
            ref = self._subject[self._snp_pos]
            assert ref.upper() in 'ACGTN'
            return ref
        except AssertionError:
            raise NoSNPError

    @get_subject_allele.add
    def get_subject_allele(self, query_snp, expected):
        """Get the reference allele"""
        genomic_position = self.get_snp_position(query_snp, expected)
        return (genomic_position, self.get_subject_allele())

    @overload
    def add_snp(self, this_snp):
        """Add a SNP to our Hsp"""
        try: # Type checking
            assert isinstance(this_snp, snp.SNP)
        except AssertionError:
            raise TypeError
        self._snp = this_snp

    @add_snp.add
    def add_snp(self, lookup):
        """Add a SNP to our HSP"""
        try: # Type checking
            assert isinstance(lookup, snp.Lookup)
        except AssertionError:
            raise TypeError
        try: # Value checking
            assert self._name == lookup.get_snpid()
            s = snp.SNP(lookup=lookup, hsp=self)
            self.add_snp(this_snp=s)
        except AssertionError:
            raise ValueError
        except:
            raise

    def get_snp(self):
        """Return the SNP associated with this Hsp"""
        if not self._snp:
            raise NoSNPError
        return self._snp


#   A function to get the value from a tag
def get_value(tag, value):
    """Get the value from an element.Tag object"""
    try:
        assert isinstance(tag, element.Tag)
        return tag.findChild(value).text
    except AssertionError:
        raise TypeError
    except:
        raise


#   A function to parse the HSP section of a BLAST XML file
def parse_hsp(hsp):
    """Parse the HSP section of a BLAST XML result"""
    try: # Type checking
        assert isinstance(hsp, element.Tag)
    except AssertionError:
        raise TypeError
    #   Ensure that our HSP has at least one mismatch
    if hsp.findChild('Hsp_midline').text.count(' ') < 1:
        raise NoSNPError
    hsp_vals = map(get_value, repeat(hsp, len(VALS)), VALS)
    return tuple(hsp_vals)


#   A function to parse the Hit section of a BLAST XML file
def parse_hit(snpid, hit):
    """Parse the Hit section of a BLAST XML result"""
    try:
        assert isinstance(snpid, str)
        assert isinstance(hit, element.Tag)
    except AssertionError:
        raise TypeError
    chrom = get_value(hit, 'Hit_def')
    vals = list()
    hsps = list()
    for hsp in hit.findAll('Hsp'):
        try:
            hsp_vals = parse_hsp(hsp)
            vals.append(hsp_vals)
        except NoSNPError:
            continue
    for val in vals:
        (bit_score, evalue, hsp_start, hsp_end, strand, identity, align_length, query, reference) = val
        hsp = Hsp(
            chrom=chrom,
            name=snpid,
            evalue=float(evalue),
            qseq=query,
            hseq=reference,
            hstart=int(hsp_start),
            hend=int(hsp_end),
            hstrand=int(strand),
            bit_score=float(bit_score),
            identity=int(identity),
            aligned_length=int(align_length)
        )
        hsps.append(hsp)
    if hsps:
        return hsps
    else:
        return False


#   A function to rank HSPs
def rank_hsps(hsps):
    """Rank HSPs and take the highest scoring (smallest) HSP(s)"""
    try:
        assert isinstance(hsps, list) or isinstance(hsps, tuple)
        for h in hsps:
            assert isinstance(h, Hsp)
    except AssertionError:
        raise TypeError
    #   Sort our HSPs into a dictionary
    hsp_dict = {}
    for h in hsps:
        if h.get_name() in hsp_dict.keys():
            hsp_dict[h.get_name()].append(h)
        else:
            hsp_dict[h.get_name()] = [h]
    #   Rank-remove our HSPs
    hsps_ranked = {name : rank_remove(hsp_list, lowest=True) for name, hsp_list in hsp_dict.items()}
    return hsps_ranked
