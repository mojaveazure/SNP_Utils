#!/usr/bin/env python3
"""Extra utilties"""

import sys
if sys.version_info.major is not 3:
    sys.exit("Please use Python 3 for this module: " + __name__)


#   A function to remove duplicates from a list
def deduplicate_list(with_dups, key):
    """Deduplicate a list, set, or tuple given a list, set, or tuple of values that should not be in the first list
        with_dups is the list with values found in key"""
    try:
        assert isinstance(with_dups, list) or isinstance(with_dups, set) or isinstance(with_dups, tuple)
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
        assert isinstance(data, list) or isinstance(data, tuple)
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
