#!/usr/bin/env python3

import sys
if sys.version_info.major is not 3:
    sys.exit("Please use Python 3 for this module: " + __name__)

import re

#   A class definition for a SAM alignment
class Alignment(object):
    """A class to hold a SAM Alignment
    It contains the following information:
        Query Name
        Bitwise Flag
        Reference sequence name
        1-based leftmost mapping position
        Cigar string
    """
    _CIGAR = re.compile(u'([0-9]+[A-Z]+)') # Regex to break the CIGAR string into component codes using re.findall()
    def __init__(self, line):
        try:
            assert isinstance(line, str)
        except AssertionError:
            raise TypeError
        #   There's only some information that we need for our alignment, everything else is forgotten
        split_line = line.strip().split() # Remove leading and trailing whitespace, then split the line by column
        self._qname = split_line[0] # First column in a SAM file
        self._flag = int(split_line[1]) # Second column
        self._rname = split_line[2] # Third column
        self._pos = int(split_line[3]) # Fourth column, should be an int
        self._cigar = self._CIGAR.findall(split_line[5]) # Sixth column, after breaking up the CIGAR string into a list of component codes

    def __repr__(self):
        return self._rname + ':' + self._qname

    def get_rc(self):
        """Check to see if we're reverse complementing our sequence"""
        #   If the 16th bit is set, it's reverse complement
        return self._flag is 16

    def get_name(self):
        """Get the alignment name"""
        return self._qname

    def get_position(self):
        """Get the alignment position"""
        return self._pos

    def get_contig(self):
        """Get the reference sequence name"""
        return self._rname

    def get_cigar(self):
        """Get the CIGAR string as a list"""
        return self._cigar

    def check_flag(self):
        """Make sure we don't have extraneous alignments"""
        return self._flag is 0 or self._flag is 16

