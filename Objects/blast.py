#!/usr/bin/env python3
"""This module holds two classes and one error"""

import sys
if sys.version_info.major is not 3:
    sys.exit("Please use Python 3 for this module: " + __name__)


from Utilities.utilities import rank_remove
from . import snp

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

#   The message if no definition line, used to decide if we use accessions instead of deflines
NO_DEF_MESSAGE = 'No definition line'

#   An error I probably overuse...
class NoSNPError(Exception):
    """A SNP has not been found"""
    def __init__(self, message): # Make this error infinitely more useful
        super(self.__class__, self).__init__(message)
        self.message = str(message)


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
        try: # Type checking
            assert isinstance(chrom, str)
            assert isinstance(name, str)
            assert isinstance(evalue, float)
            assert isinstance(qseq, str)
            assert isinstance(hseq, str)
            assert isinstance(hstart, int)
            assert isinstance(hend, int)
            assert isinstance(hstrand, int)
            assert isinstance(bit_score, (int, float))
            assert isinstance(identity, int)
            assert isinstance(aligned_length, int)
        except AssertionError:
            raise TypeError
        try: # Value checking
            assert hstrand == 1 or hstrand == -1
            assert identity <= aligned_length
        except AssertionError:
            raise ValueError
        #   Hsp information
        self._chrom = chrom # Hit_def/Hit_accession
        self._name = name # Iteration_query-def
        self._evalue = evalue # Hsp_evalue
        self._query = qseq # Hsp_qseq
        self._subject = hseq # Hsp_hseq
        self._start = hstart # Hsp_hit-from
        self._end = hend # Hsp_hit-to
        self._hstrand = hstrand #Hsp_hit-strand
        self._bits = bit_score # Hsp_bit-score
        self._identity = identity # Hsp_identity
        self._alength = aligned_length # Hsp_align-len
        #   SNP information
        self._snp_pos = None
        self._snp = None

    def __repr__(self):
        return self._name + ":" + str(self._evalue)

    def __eq__(self, other):
        if isinstance(other, Hsp): # Check equality with other Hsp objects
            name_bool = self._name == other._name # Ensure our name is equal to the other name
            eval_bool = self._evalue == other._evalue # Is our evalue equal to the other evalue
            score_bool = self._bits == other._bits # Is our score equal to the other score
            return name_bool and eval_bool and score_bool # Return whether all of these are true
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

    @overload
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
            print('noquery', self, file=sys.stderr)
            raise NoSNPError # Error out
        #   If the hsp sequence is less than our expected value
        elif (len(self._query) < expected) or (self._query.find(query_snp) <= expected and self._query.count(query_snp) == 1):
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
                self._snp_pos = expected + qsubstr.count('-') # Add any deletions to our position
                if old is self._snp_pos: # If we didn't change our position guess
                    print("loop error", self, file=sys.stderr)
                    raise NoSNPError # Error out
                qsubstr = self._query[:self._snp_pos + 1] # Replace our substring
            #   Calculate insertions
            hsubstr = self._subject[:self._snp_pos + 1]
            insertion = hsubstr.count('-')
            if self.get_rc():
                return self._start - self._snp_pos - insertion
            else:
                return self._start + self._snp_pos - insertion

    @get_snp_position.add
    def get_snp_position(self, lookup):
        """Calculate the SNP position relative to reference"""
        try:
            assert isinstance(lookup, snp.Lookup)
        except AssertionError:
            raise TypeError
        try:
            assert lookup.get_snpid() == self._name
        except AssertionError:
            raise ValueError
        iupac = lookup.get_sequence(iupac=True).upper()
        query = self._query.replace('-', '').upper()
        if self._alength >= lookup.get_length() or iupac.find(query) is 0:
            return self.get_snp_position(
                query_snp=lookup.get_code(),
                expected=lookup.get_forward_position()
            )
        elif self._alength < lookup.get_length():
            expected = iupac.find(lookup.get_code(), lookup.get_forward_position())
            aligned = iupac.find(query)
            if expected > aligned:
                adj = expected - aligned
            elif expected < aligned:
                adj = lookup.get_forward_position() + (aligned - expected) + query.find(lookup.get_code())
            return self.get_snp_position(
                query_snp=lookup.get_code(),
                expected=adj
            )
        else:
            raise ValueError("Something happened...")

    @overload
    def get_subject_allele(self):
        """Get the reference allele"""
        try:
            assert self._snp_pos is not None
            ref = self._subject[self._snp_pos]
            assert ref.upper() in 'ACGTN'
            return ref
        except AssertionError:
            print('noref', self, file=sys.stderr)
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
        print('SNP!', self, file=sys.stderr)
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
            # s = snp.SNP(lookup=lookup, hsp=self)
            # self.add_snp(this_snp=s)
            self.add_snp(this_snp=snp.SNP(lookup=lookup, hsp=self))
        except AssertionError:
            raise ValueError
        except:
            print('hmm', self, file=sys.stderr)
            raise

    def get_snp(self):
        """Return the SNP associated with this Hsp"""
        if not self._snp:
            raise NoSNPError('No snp for me: ' + self)
        return self._snp


#   A function to get the value from a tag
def get_value(tag, value):
    """Get the value from an element.Tag object"""
    try:
        assert isinstance(tag, element.Tag)
        assert isinstance(value, str)
        return tag.findChild(value).text
    except AssertionError:
        raise TypeError("'tag' must be of type 'element.tag'; 'value' must be of type 'str'")
    except AttributeError:
        return ''
    except:
        raise


#   A function to parse the HSP section of a BLAST XML file
def parse_hsp(hsp):
    """Parse the HSP section of a BLAST XML result"""
    try: # Type checking
        assert isinstance(hsp, element.Tag)
    except AssertionError:
        raise TypeError("'hsp' must be of type 'element.Tag'")
    #   Ensure that our HSP has at least one mismatch
    if hsp.findChild('Hsp_midline').text.count(' ') < 1:
        raise NoSNPError('No mismatches found for this Hsp')
    return tuple(get_value(tag=hsp, value=val) for val in VALS)


#   A function to parse the Hit section of a BLAST XML file
def parse_hit(snpid, hit):
    """Parse the Hit section of a BLAST XML result"""
    try: # Type checking
        assert isinstance(snpid, str)
        assert isinstance(hit, element.Tag)
    except AssertionError:
        raise TypeError("'snpid' must be of type 'str'; 'hit' must be of type 'element.Tag'")
    #   Get and check the chromosome information
    chrom = get_value(tag=hit, value='Hit_def')
    if chrom == NO_DEF_MESSAGE:
        chrom = get_value(tag=hit, value='Hit_accession')
    hsps = list() # Create a holding list
    for hsp in hit.findAll('Hsp'):
        try:
            hsp_vals = parse_hsp(hsp)
        except NoSNPError as nosnp:
            print(snpid + ':', nosnp.message)
        else:
            (bit_score, evalue, hsp_start, hsp_end, strand, identity, align_length, query, reference) = hsp_vals # Unpack our tuple
            hsp = Hsp( # Create an Hsp object
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
            hsps.append(hsp) # Add to our list
    if hsps:
        return hsps
    else:
        return None


#   A function to rank HSPs
def rank_hsps(hsps):
    """Rank HSPs and take the highest scoring (smallest) HSP(s)"""
    try:
        assert isinstance(hsps, (list, tuple))
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
