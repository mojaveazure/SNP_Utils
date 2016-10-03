#!/usr/bin/env python3
"""A module for SNP stuff"""

import sys
if sys.version_info.major is not 3:
    sys.exit("Please install Python 3 for this module: " + __name__)

import re

from . alignment import Alignment
from . import blast
from . blast import NoSNPError

try:
    from overload import overload
except ImportError as error:
    print("Please install " + error.name)


#   An error to say we weren't given a base
class NotABaseError(Exception):
    """You have not provided a single-character: [ACGTacgt]"""


#   An error to say that our given objects do not match
class NoMatchError(Exception):
    """The objects you have provided do not match"""


#   A class to hold SNP information
class SNP(object):
    """A class to hold VCF information about an individual SNP. Stores the
    following information:
        SNP ID
        Contig
        Reference Position
        Reference Base
        Alternate Base
    """
    @staticmethod
    def reverse_complement(base):
        """Get the reverse complement of a nucleotide"""
        try:
            assert isinstance(base, str)
            assert len(base) is 1
            assert base in 'ACGTNacgtn'
            rc = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
            return base.translate(rc)
        except AssertionError:
            raise NotABaseError

    @overload
    def __init__(self, lookup):
        try:
            assert isinstance(lookup, Lookup)
        except AssertionError:
            raise TypeError
        self._snpid = lookup.get_snpid()
        #   Everything else gets made in a bit
        self._contig = None
        self._position = None
        self._reference = None
        self._alternate = None

    @__init__.add
    def __init__(self, lookup, alignment, reference):
        try:
            assert isinstance(alignment, Alignment)
            assert isinstance(reference, dict)
            self.__init__(lookup)
        except:
            raise TypeError
        if self._snpid != alignment.get_name():
            raise NoMatchError
        #   The contig is found in the alignment
        self._contig = alignment.get_contig()
        #   Get the rest of the information
        self._calculate_position(lookup, alignment) # True SNP position
        self._find_states(lookup, alignment, reference) # Reference and alternate states

    @__init__.add
    def __init__(self, lookup, hsp):
        try:
            assert isinstance(hsp, blast.Hsp)
            self.__init__(lookup)
        except:
            raise TypeError
        if self._snpid != hsp.get_name():
            raise NoMatchError
        #   The contig and SNP position are found in the hsp
        self._contig = hsp.get_chrom()
        self._position = hsp.get_snp_position(lookup.get_code(), lookup.get_forward_position())
        #   Get the rest of the information
        self._find_states(lookup, hsp) # Reference and alternate states

    def __repr__(self):
        return self._snpid

    def __eq__(self, other):
        if isinstance(other, SNP):
            return self._snpid == other._snpid
        elif isinstance(other, Lookup):
            return self._snpid == other._snpid and self._alternate == other.get_alternate(self._reference) and self._alternate != 'N'
        elif isinstance(other, str):
            return self._snpid == other
        elif isinstance(other, int):
            return self._position == other
        else:
            return NotImplemented

    def __hash__(self):
        return hash(self._snpid)

    def _calculate_position(self, lookup, alignment):
        """Calculate the position of the SNP in the reference sequence"""
        index = 0 # Index of our split CIGAR string
        if alignment.get_rc() or lookup.get_rc(): # If we're reverse complementing
            qpos = lookup.get_reverse_position() - 1 # Start with the reverse position of the SNP, must subtract one
        else: # Otherwise
            qpos = lookup.get_forward_position() # Start with the forward posittion
        while True: # Endless loop to do weird things...
            try: # While we have a CIGAR string to parse
                old = qpos # Store our previously calculated SNP position
                #   Seach the CIGAR string as a list, starting with index 0, for indels
                if re.search('M', alignment.get_cigar()[index]): # If we have a perfect match
                    if qpos < int(''.join(re.findall(r'\d+', alignment.get_cigar()[index]))): # If our SNP is in the perfect match
                        break # Exit the loop, we have our position
                if re.search('D', alignment.get_cigar()[index]): # If we have a deletion relative to reference
                    qpos += int(''.join(re.findall(r'\d+', alignment.get_cigar()[index]))) # Add the deletion to our SNP position
                if re.search('[IS]', alignment.get_cigar()[index]): # If we have an insertion relative to reference
                    qpos -= int(''.join(re.findall(r'\d+', alignment.get_cigar()[index]))) # Subtract the insertion from our SNP postion
                index += 1 # Increase the index
                if qpos <= 0 or qpos >= lookup.get_length(): # If we've gone beyond the scope of our lookup: 0 is before the sequence, lookup.get_length() is after
                    qpos = old # Go back to our previously calculated SNP postion
                    break # Exit the loop, we have our position
            except IndexError: # If we run out of CIGAR string codes
                break # Exit the loop, we have our position
        self._position = alignment.get_position() + qpos # Our SNP position is at the mapping position plus the SNP position

    @overload
    def _find_states(self, lookup, alignment, reference):
        """Get the reference and alternate alleles"""
        #   Get the reference allele, given our contig and position found above
        self._reference = reference[self._contig][self._position - 1] # Subtract one as FASTA is 1-based and Python is 0-based
        if alignment.get_rc(): # If we're reverse complement
            alt = lookup.get_alternate(self.reverse_complement(self._reference))
            self._alternate = self.reverse_complement(alt)
        else:
            self._alternate = lookup.get_alternate(self._reference) # An 'N' will be returned if the reference allele doesn't match with our IUPAC code

    @_find_states.add
    def _find_states(self, lookup, hsp):
        """Get the reference and alternate alleles"""
        try:
            assert isinstance(lookup, Lookup)
            assert isinstance(hsp, blast.Hsp)
            self._reference = hsp.get_subject_allele() # Get the reference allele from the Hit
            self._alternate = lookup.get_alternate(self._reference) # Get the alternate from the Lookup
            if hsp.get_rc():
                self._reference = self.reverse_complement(self._reference)
                self._alternate = self.reverse_complement(self._alternate)
        except AssertionError:
            raise TypeError
        except NoSNPError:
            raise
        except NotABaseError:
            raise

    def get_snpid(self):
        """Get the SNP ID"""
        return self._snpid

    def get_chrom(self):
        """Get the chromosome/contig"""
        return self._contig

    def get_position(self):
        """Get the position of the SNP"""
        return self._position

    def check_masked(self):
        """Check to see if our alternate allele is masked"""
        return self._alternate == 'N'

    def format_vcf(self):
        """Format the information in VCF style"""
        #   Create a list of information for a VCF file
        vcf_line = [
            self._contig,
            str(self._position),
            self._snpid,
            self._reference,
            self._alternate,
            '.',
            '.',
            's'
            ]
        return '\t'.join(vcf_line) # Join everything together with a tab


#   A class definition for a lookup sequence in Illumina format
class Lookup(object):
    """This is a class for a SNP lookup sequence in Illumina format
    It contains the following information:
        SNP ID
        Sequence with Illumina syntax
        Sequence in IUPAC codes
        Sequence Length
        SNP Position from forward
        SNP Position from reverse
    """

    # A dictionary of IUPAC codes for SNPs
    _IUPAC_CODES = {
        'R' : 'AG',
        'Y' : 'CT',
        'S' : 'CG',
        'W' : 'AT',
        'K' : 'GT',
        'M' : 'AC'
    }

    def __init__(self, snpid, sequence):
        #   We're given the SNP ID and sequence when making the object, everything else 
        #   can be made with _capture_snp() and _find_iupac()
        self._snpid = snpid
        self._sequence = sequence
        self._forward_position = None
        self._reverse_position = None
        self._snp = None
        self._code = None
        self._iupac = None
        #   Get the rest of the information we need for our lookup
        self._capture_snp()
        self._find_iupac()

    def __repr__(self):
        return self._snpid + ':' + self._code

    def __eq__(self, other):
        if isinstance(other, Lookup):
            return self._snp == other._snpid and self._code == other._code
        elif isinstance(other, str):
            if len(other) is 1:
                return self._code == other.upper()
            elif len(other) > 1:
                return self._snpid == other
            else:
                return False
        else:
            return NotImplemented

    def _capture_snp(self):
        """Capture the SNP and it's position from the start and end of the sequence"""
        #   Get the forward position
        self._forward_position = self._sequence.find('[')
        #   Get the reverse position
        self._reverse_position = len(self._sequence) - self._sequence.find(']')
        #   Get the SNP
        self._snp = self._sequence[self._forward_position:self._sequence.find(']') + 1]

    def _find_iupac(self):
        """Create an IUPAC version of the sequence and calculate it's length"""
        #   Create a string of the two states of the SNP in alphabetical order
        ordered_snp = ''.join(sorted(re.findall('[ACGT]', self._snp)))
        #   Find the IUPAC code for the SNP
        self._code = ''.join([c for c, o in self._IUPAC_CODES.items() if ordered_snp == o])
        #   Create the IUPAC version of the sequence
        self._iupac = re.sub(r'\[%s\]' % self._snp[1:-1], self._code, self._sequence)

    def get_snpid(self):
        """Get the SNP ID"""
        return self._snpid

    def get_forward_position(self):
        """Get the SNP position in the forward direction"""
        return self._forward_position

    def get_reverse_position(self):
        """Get the SNP position in the reverse direction"""
        return self._reverse_position

    def get_code(self):
        """Get the IUPAC code"""
        return self._code

    def get_length(self):
        """Get the length of the IUPAC sequence"""
        return len(self._iupac)

    #   Search the IUPAC codes for an alternate allele of a SNP
    def get_alternate(self, reference):
        """Get the alternate allele given an IUPAC code and reference allele"""
        ref = re.compile(u'(%s)' % reference) # Regex to ensure that our found reference allele is covered by the IUPAC code
        alt = re.compile(u'([^%s])' % reference) # Regex to find the alternate allele
        if ref.search(self._IUPAC_CODES[self._code]): # If our reference allele is plausible given our IUPCA code
            alternate = alt.search(self._IUPAC_CODES[self._code]).group() # Get the alternate
            return alternate
        else: # Otherwise, give an 'N'
            return 'N'

    def format_fasta(self):
        """Return this lookup in FASTA format"""
        fasta = [
            '>' + self._snpid,
            self._iupac
        ]
        return '\n'.join(fasta)
