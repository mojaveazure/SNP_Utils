#!/usr/bin/env python3
"""Extra utilties"""

import sys
if sys.version_info.major is not 3:
    sys.exit("Please use Python 3 for this module: " + __name__)


import os
from datetime import datetime

#   A class to hold values for an ##INFO declaration
class INFO(object):
    """A class for a VCF ##INFO declaration
    It holds the following information:
        INFO ID             (str)
        INFO Number         (int or str)
        INFO Type           (str)
        INFO Description    (str)
    """

    _VALID_TYPES = ['Integer', 'Float', 'Flag', 'Character', 'String']
    _CHAR_NUMS = 'AGR.'

    def __init__(self, infoid, number, infotype, description):
        try: # Type checking
            assert isinstance(infoid, str)
            assert isinstance(number, (int, str))
            assert isinstance(infotype, str)
            assert isinstance(description, str)
        except AssertionError:
            raise TypeError
        try: # Value checking
            assert infotype in self._VALID_TYPES
            if infotype == 'Flag':
                assert int(number) is 0
            if isinstance(number, str):
                assert len(number) is 1
                assert number.upper() in self._CHAR_NUMS
        except AssertionError:
            raise ValueError
        self._infoid = infoid
        self._number = str(number).upper()
        self._type = infotype
        self._description = description

    def __call__(self):
        """Format the VCF ##INFO declaration for printing"""
        vals = [ # Assemble the meat of the INFO declaration
            'ID=' + self._infoid,
            'Number=' + self._number,
            'Type=' + self._type,
            'Description="' + self._description + '"'
        ]
        info = '##INFO=<' + ','.join(vals) + '>' # Create the full declaration
        return info

    def get_id(self):
        """Get the ID for this INFO"""
        return self._infoid

    def get_number(self):
        """Get the number for this INFO"""
        return self._number

    def get_type(self):
        """Get the type for this INFO"""
        return self._type

    def get_description(self):
        """Get the description for this INFO"""
        return self._description


#   A function to remove duplicates from a list
def deduplicate_list(with_dups, key):
    """Deduplicate a list, or tuple given a list, set, or tuple of values that should not be in the first list
        with_dups is the list with values found in key"""
    try:
        assert isinstance(with_dups, (list, tuple))
        assert isinstance(key, list) or isinstance(key, set) or isinstance(key, tuple)
    except AssertionError:
        raise
    deduped = set()
    for value in with_dups:
        if value not in key:
            deduped.add(value)
    print("Removed", len(with_dups) - len(deduped), "duplicates", file=sys.stderr)
    return deduped


#   A function to rank values in a list or tuple and keep only the highest/lowest value(s)
def rank_remove(data, lowest=False):
    """Rank values of a list or tuple and take only the highest (or lowest)"""
    try:
        #   Ensure our data is a list or tuple
        assert isinstance(data, (list, tuple))
        assert isinstance(lowest, bool)
    except AssertionError:
        raise TypeError("'data' must be of type 'list' or 'tuple', 'lowest' must be of type 'bool'")
    #   Create our first result
    try:
        result = [data[0]]
    except KeyError: # There must be at least one data entry
        raise ValueError("'data' must have at least one value")
    try:
        for value in data[1:]: # For every other data point
            if lowest: # If we're looking for the lowest rank
                if value < result[0]: # If lower
                    result = [value] # Replace
                elif value == result[0]: # If equal
                    result.append(value) # Add
                else: # Otherwise
                    continue # Pass
            else: # If we're looking for the highest rank
                if value > result[0]: # If greater
                    result = [value] # Replace
                elif value == result[0]: # If equal
                    result.append(value) # Add
                else: # Otherwise
                    continue # Pass
    except KeyError: # If there is only one data point (taken above)
        pass # Skip this, we return the one point
    except: # Other exceptions
        raise # Let someone else deal with it
    return result # Return our list


#   A function to make a VCF header
def vcf_header(ref, info=None):
    """Make a header for a VCF file given a reference genome and an optional list of INFO objects"""
    try: # Type checking
        assert isinstance(ref, str)
        assert isinstance(info, list) or info is None
        for i in info:
            assert isinstance(i, INFO)
    except TypeError:
        pass # Ignore TypeErrors if 'info' is None
    except AssertionError: # Catch assertions and raise as TypeError
        raise TypeError
    #   The first several header lines
    fileformat = '##fileformat=VCFv4.2'
    filedate = '##fileDate=' + datetime.now().strftime("%Y%m%d")
    source = '##source=' + os.path.basename(sys.argv[0])
    reference = '##reference=' + os.path.basename(ref)
    colnames = [
        '#CHROM',
        'POS',
        'ID',
        'REF',
        'ALT',
        'QUAL',
        'FILTER',
        'INFO'
    ]
    header = [
        fileformat,
        filedate,
        source,
        reference
    ]
    #   Add lines from INFO objects
    if info:
        header.extend([i() for i in info])
    header.append('\t'.join(colnames))
    return '\n'.join(header)
