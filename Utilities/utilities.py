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
