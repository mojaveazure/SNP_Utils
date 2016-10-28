#!/usr/bin/env python3

"""This module contains various functions to filter potential SNPs in hopes of removing duplicate potential SNPs"""

import sys
if sys.version_info.major is not 3:
    sys.exit("Please use Python 3 for this module: " + __name__)


from copy import deepcopy

from Objects.snp import SNP, read_map, filter_snps
from Objects.wrappers import NcbiblastdbcmdCommandline
from . utilities import INFO, closest_value

try:
    from Bio import SeqIO
    from Bio.Application import ApplicationError
except ImportError as error:
    sys.exit("Please install " + error.name)

#   Two INFO objects for alternate SNP positions information
ALTCHR = INFO( # Create an INFO object for alternate chromosomes
    infoid='ALTCHR',
    number='.',
    infotype='String',
    description='Alternate chromosomal positions'
)
ALTPOS = INFO( # Create an info object for alternate positions
    infoid='ALTPOS',
    number='.',
    infotype='Integer',
    description='Alternate positions'
)

#   A function to sort a list of SNPs by ID and optionally chromosome/contig
def _sort_snps(snplist, bychrom=False):
    """Sort SNPs by SNP ID and, optionally, chromosome/contig"""
    try:
        assert isinstance(snplist, (list, tuple))
        for snp in snplist:
            assert isinstance(snp, SNP)
        assert isinstance(bychrom, bool)
    except AssertionError:
        raise TypeError
    snp_dict = dict()
    for snp in snplist:
        if bychrom:
            key = (snp.get_snpid(), snp.get_chrom())
        else:
            key = snp.get_snpid()
        if key in snp_dict:
            snp_dict[key].append(snp)
        else:
            snp_dict[key] = [snp]
    return snp_dict


#   A function to filter potential SNPs by chromosomal/contig
def chrom_filter(args, vcfinfo, propersnps):
    """Filter potential SNPs by chromosomes/contigs using a genetic map"""
    #   Create holding collections
    mapped_snps = list() # List of mapped SNPs
    try: # Type checking
        assert isinstance(args, dict)
        assert isinstance(vcfinfo, list)
        for info in vcfinfo:
            assert isinstance(info, INFO)
        assert isinstance(propersnps, (filter, list, tuple))
        #   Sort our SNPs into several lists of the same marker
        snp_dict = _sort_snps(snplist=propersnps)
    except AssertionError:
        raise TypeError
    except:
        raise
    try: # Value checking
        assert 'map' in args.keys()
    except AssertionError:
        raise ValueError
    print("Filtering SNPs by hit chromsome/contig", file=sys.stderr)
    #   Add our VCF header
    infoheader = set(deepcopy(vcfinfo)) # NO SIDE EFFECTS
    infoheader.add(ALTCHR)
    infoheader.add(ALTPOS)
    #   Open our map file
    map_dict = read_map(args['map'])
    #   Create our list of mapped snps
    for snpid, snplist in snp_dict.items():
        #   Only filter by map if there are multiple SNPs for the same marker and we have a map position for the marker
        if len(snplist) > 1 and snpid in map_dict.keys():
            mapped_snps += filter_snps(
                snp_list=snplist,
                genetic_map=map_dict[snpid],
                altchr=ALTCHR.get_id(),
                altpos=ALTPOS.get_id()
            )
        else: # Otherwise, add the SNP list to our list of mapped SNPs
            mapped_snps += snplist
    return (mapped_snps, list(infoheader))


#   A function to filter by a distance threshold
def distance_filter(args, vcfinfo, propersnps):
    """Filter potential SNPs by a distance threshold"""
    #   Create holding collections
    filtered_snps = list() # List of filtered SNPs
    # mapped_snps = list() # List of mapped SNPs
    try: # Type checking
        assert isinstance(args, dict)
        assert isinstance(vcfinfo, list)
        for info in vcfinfo:
            assert isinstance(info, INFO)
        assert isinstance(propersnps, (filter, list, tuple))
        #   Sort our SNPs into several lists of the same marker
        snp_dict = _sort_snps(snplist=propersnps, bychrom=True)
    except AssertionError:
        raise TypeError
    except:
        raise
    try: # Value checking
        assert 'threshold' in args.keys()
        assert isinstance(args['threshold'], int) and args['threshold'] >= 0
    except AssertionError:
        raise ValueError
    print("Filtering SNPs with a minimum distance threshold of", args['threshold'], file=sys.stderr)
    #   Add our VCF header
    infoheader = set(deepcopy(vcfinfo)) # NO SIDE EFFECTS
    infoheader.add(ALTCHR)
    infoheader.add(ALTPOS)
    #   Start filtering by distance threshold
    for snplist in snp_dict.values():
        #   Sort our SNP list (SNPs have a __lt__ attribute based on position so this works as is)
        sorted_snps = sorted(snplist)
        left_snps = [sorted_snps.pop(0)] # First SNP in sorted list is the leftmost SNP
        while sorted_snps: # The while loop stops when we're done with our list of sorted SNPs
            left_index = len(left_snps) - 1 # Take the furthest right of our left SNPs
            next_snp = sorted_snps.pop(0) # Get the next SNP
            if next_snp - left_snps[left_index] < args['threshold']: # If the distance between the two SNPs are less than our minimum distance threshold
                #   Add the next SNP as an alternate position to our left SNP
                left_snps[left_index].add_info(
                    key=ALTCHR.get_id(),
                    value=next_snp.get_chrom()
                )
                left_snps[left_index].add_info(
                    key=ALTPOS.get_id(),
                    value=next_snp.get_position()
                )
            else: # Otherwise, add the next SNP to our list of left-most SNPs
                left_snps.append(next_snp)
        filtered_snps += left_snps # Add our list of left-most SNPs to the overall SNP list
    return (filtered_snps, list(infoheader))


#   Use the genetic map to finish the filtering
def map_filter(args, vcfinfo, propersnps, refgen, method, bconf=None):
    """Filters potatoes"""
    #   Create holding collections
    chrom_lengths = dict() # Dictionary to hold chromosomal sequence lengths
    closest_snps = list() # List of filtered SNPs
    try: # Type checking
        assert isinstance(args, dict)
        assert isinstance(vcfinfo, list)
        for info in vcfinfo:
            assert isinstance(info, INFO)
        assert isinstance(propersnps, (filter, list, tuple))
        #   Sort our SNPs into several lists of the same marker
        snp_dict = _sort_snps(snplist=propersnps, bychrom=True)
        assert isinstance(refgen, str)
        assert isinstance(method, str)
        assert isinstance(bconf, dict) or bconf is None
    except AssertionError:
        raise TypeError
    try: # Value checking
        assert 'map' in args
        assert method in {'BLAST', 'SAM'}
        if method == 'BLAST':
            assert bconf is not None
        if bconf:
            assert 'subject' in bconf or 'database' in bconf
    except AssertionError:
        raise ValueError
    print("Filtering SNPs by relative location on the genetic map", file=sys.stderr)
    #   Add our VCF header
    infoheader = set(deepcopy(vcfinfo)) # NO SIDE EFFECTS
    infoheader.add(ALTCHR)
    infoheader.add(ALTPOS)
    #   Read in the genetic map
    map_dict = read_map(args['map'])
    #   Get the length of each chromosome in the reference
    if method == 'SAM' or 'subject' in bconf:
        try:
            print("Reading in reference FASTA", refgen, "This might take a while...", file=sys.stderr)
            reference = SeqIO.to_dict(SeqIO.parse(refgen, 'fasta'))
            for chrom, seq in reference.items():
                chrom_lengths[chrom] = len(seq)
        except FileNotFoundError as error:
            sys.exit("Failed to find " + error.filename)
    elif 'database' in bconf:
        try:
            print("Parsing blast database", bconf['database'], file=sys.stderr)
            accessioncmd = NcbiblastdbcmdCommandline(
                db=bconf['database'],
                dbtype='nucl',
                entry='all',
                outfmt='%a'
            )
            accout = accessioncmd()[0]
            accessions = (acc for acc in accout.strip().split())
            for acc in accessions:
                lengthcmd = NcbiblastdbcmdCommandline(
                    db=bconf['database'],
                    dbtype='nucl',
                    entry=acc,
                    outfmt='%l'
                )
                length = lengthcmd()[0]
                chrom_lengths[acc] = int(length)
        except ApplicationError as error:
            sys.exit(error)
    else:
        raise NotImplementedError
    print("Found", len(chrom_lengths), 'chromosomes', file=sys.stderr)
    print("Filtering", len(snp_dict), "SNP IDs", file=sys.stderr)
    for key, snplist in snp_dict.items():
        (snpid, chrom) = key
        try:
            assert len(snplist) > 1
            these_snps = deepcopy(snplist)
            map_length = max([m.get_map_distance() for m in map_dict.values() if m.get_chrom() == chrom])
            map_proportion = map_dict[snpid].get_map_distance() / map_length
            length_dict = {snp / chrom_lengths[chrom] : snp for snp in these_snps}
            closest = closest_value(iterable=length_dict, value=map_proportion, return_item=True)
            these_snps.remove(closest)
            for remaining in these_snps:
                closest.add_info(
                    key=ALTCHR.get_id(),
                    value=remaining.get_chrom()
                )
                closest.add_info(
                    key=ALTPOS.get_id(),
                    value=remaining.get_position()
                )
            closest_snps.append(closest)
        except AssertionError: # Don't need to do all this crap and spit out error messages
            closest_snps.extend(snplist)
        except KeyError:
            print("No map found for", snpid, file=sys.stderr)
            closest_snps.extend(snplist)
        except ValueError:
            print("No genetic map for chromosome", chrom, file=sys.stderr)
            closest_snps.extend(snplist)
    return (closest_snps, list(infoheader))







