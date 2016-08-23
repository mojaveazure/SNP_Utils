#!/usr/bin/env python3
"""A script to find SNP positions from a BLAST search"""

import sys
if sys.version_info.major is not 3:
    sys.exit("Please use Python 3 for this script")

from Utilities import arguments
from Objects import snp
from Objects import blast
from Objects.blast import NoSNPError
from BLAST import runblastn
from BLAST import configure

try:
    from bs4 import BeautifulSoup
    from bs4 import FeatureNotFound
    from Bio import SeqIO
except ImportError as error:
    sys.exit("Please install " + error.name)


#   BLAST-based
def blast_based(args, lookup_dict):
    try:
        assert isinstance(args, dict)
        assert isinstance(lookup_dict, dict)
    except AssertionError:
        raise
    if not args['xml']:
        print("Running BLAST", file=sys.stderr)
        bconf = configure.parse_config(args['config'])
        bconf['query'] = runblastn.write_fasta(lookup_dict, args['lookup'])
        blast_out = runblastn.run_blastn(bconf)
        blast_xml = open(blast_out, 'r')
    else:
        print("Loading BLAST XML file", file=sys.stderr)
        blast_xml = args['xml']
    try:
        blast_soup = BeautifulSoup(blast_xml, 'xml')
    except FeatureNotFound:
        sys.exit("Please install 'lxml' to properly parse the BLAST results")
    iteration_dictionary = {}
    no_hit = []
    snp_list = []
    for query in blast_soup.findAll('Iteration'):
        iteration = blast.SNPIteration(query)
        iteration_dictionary[iteration.get_snpid()] = iteration
    for snpid, l in lookup_dict.items():
        if snpid in iteration_dictionary:
            try:
                snps, fail = iteration_dictionary[snpid].hit_snps(l)
                snp_list += snps
                no_hit += fail
            except NoSNPError:
                print("No hit found for", snpid, file=sys.stderr)
                no_hit.append(l.get_snpid())
    blast_xml.close()
    return(snp_list, no_hit)


#   Write the output files
def write_outputs(args, snp_filter, masked_filter, no_snps):
    """Write the output files"""
    try:
        assert isinstance(snp_filter, filter)
        assert isinstance(masked_filter, filter)
    except AssertionError:
        raise
    outfile = args['outname'] + '.vcf'
    maskedfile = args['outname'] + '_masked.vcf'
    failedfile = args['outname'] + '_failed.log'
    snp_list = list(snp_filter)
    if len(snp_list) < 1:
        print("Failed to find any SNPs", file=sys.stderr)
    else:
        print("Writing output to", outfile, file=sys.stderr)
        with open(outfile, 'w') as o:
            for s in snp_list:
                o.write(s.format_vcf())
                o.write('\n')
    masked_list = list(masked_filter)
    if len(masked_list) > 0:
        print("Writing", len(masked_list), "masked SNPs to", maskedfile, file=sys.stderr)
        with open(maskedfile, 'w') as m:
            for s in masked_list:
                m.write(s.format_vcf())
                m.write('\n')
    if len(no_snps) > 0:
        print("Failed to find", len(no_snps), "SNPs", file=sys.stderr)
        with open(failedfile, 'w') as f:
            for s in no_snps:
                f.write(s)
                f.write('\n')

#   Run the program
def main():
    """Find SNP positions"""
    parser = arguments.make_argument_parser()
    if not sys.argv[1:]:
        sys.exit(parser.print_help())
    args = vars(parser.parse_args())
    if args['method'] == 'CONFIG':
        configure.make_config(args)
    else:
        lookup_dict = {}
        with open(args['lookup'], 'r') as lk:
            for line in lk:
                split = line.strip().split()
                l = snp.Lookup(split[0], split[1])
                lookup_dict[l.get_snpid()] = l
        if args['method'] == 'BLAST':
            snp_list, no_snps = blast_based(args, lookup_dict)
        elif args['method'] == 'SAM':
            raise NotImplementedError("Finding SNPs using a SAM file is not yet implemented")
        masked = filter(lambda s: s.check_masked(), snp_list)
        proper_snps = filter(lambda s: not s.check_masked(), snp_list)
        write_outputs(args, proper_snps, masked, no_snps)


if __name__ == '__main__':
    main()
