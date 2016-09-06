#!/usr/bin/env python3
"""A script to find SNP positions from a BLAST search"""

import sys
if sys.version_info.major is not 3:
    sys.exit("Please use Python 3 for this script")

from Utilities import arguments
from Utilities import utilities
from Objects import snp
from Objects import blast
from Objects.snp import NotABaseError
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
    #   If we weren't provided a BLAST XML file
    if not args['xml']:
        #   Run BLASTn
        print("Running BLAST", file=sys.stderr)
        #   Configure BLAST
        bconf = configure.parse_config(args['config'])
        #   Make our FASTA file
        bconf['query'] = runblastn.write_fasta(lookup_dict, args['lookup'])
        #   Run BLAST
        blast_out = runblastn.run_blastn(bconf=bconf)
        #   Open the XML file
        blast_xml = open(blast_out, 'r')
    else:
        #   Read the XML file provided
        print("Loading BLAST XML file", file=sys.stderr)
        blast_xml = args['xml']
    try:
        #   Make soup out of the XML
        blast_soup = BeautifulSoup(blast_xml, 'xml')
    except FeatureNotFound:
        sys.exit("Please install 'lxml' to properly parse the BLAST results")
    iteration_dictionary = {} # Holds iterations stored by SNP ID
    no_hit = [] # List of SNP IDs
    snp_list = [] # List of snp.SNPs
    #   For every query in our BLAST soup
    for query in blast_soup.findAll('Iteration'):
        iteration = blast.SNPIteration(query) # Make an iteration
        iteration_dictionary[iteration.get_snpid()] = iteration # Add to our dictionary
    #   For every SNP in our lookup table
    for snpid, l in lookup_dict.items():
        try:
            if snpid in iteration_dictionary: # If the SNP was found in BLAST
                snps, fail = iteration_dictionary[snpid].hit_snps(l) # Get SNPs and failures for every hit for this iteration
                snp_list += snps # Add to our list of SNPs
                no_hit += fail # Add to our list of failures
            else:
                raise NoSNPError
        except NoSNPError:
            print("No hit found for", snpid, file=sys.stderr)
            no_hit.append(snpid)
        except NotABaseError:
            print("Something happened with", snpid, file=sys.stderr)
            no_hit.append(snpid)
    #   Close the XML file
    blast_xml.close()
    return(snp_list, no_hit)


#   Write the output files
def write_outputs(args, snp_filter, masked_filter, no_snps, method):
    """Write the output files"""
    try:
        assert isinstance(args, dict)
        assert isinstance(snp_filter, filter)
        assert isinstance(masked_filter, filter)
        assert isinstance(no_snps, list)
        assert isinstance(method, str)
    except AssertionError:
        raise
    header = '##fileformat=VCFv4.2\n##INFO<ID=s,Number=1,Type=Flag,Description="Variant is calculated from %s">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' % method
    outfile = args['outname'] + '.vcf'
    maskedfile = args['outname'] + '_masked.vcf'
    failedfile = args['outname'] + '_failed.log'
    snp_list = list(snp_filter)
    if len(snp_list) < 1:
        print("Failed to find any SNPs", file=sys.stderr)
    else:
        print("Writing", len(snp_list), "SNPs to", outfile, file=sys.stderr)
        with open(outfile, 'w') as o:
            o.write(header)
            for s in snp_list:
                o.write(s.format_vcf())
                o.write('\n')
    print("Removing masked SNPs that were actually found", file=sys.stderr)
    masked_list = utilities.deduplicate_list(list(masked_filter), snp_list)
    if len(masked_list) > 0:
        print("Writing", len(masked_list), "masked SNPs to", maskedfile, file=sys.stderr)
        with open(maskedfile, 'w') as m:
            m.write(header)
            for s in masked_list:
                m.write(s.format_vcf())
                m.write('\n')
    print("Removing failed SNPs that were actually found", file=sys.stderr)
    no_snps_deduped_parital = utilities.deduplicate_list(no_snps, [s.get_snpid() for s in snp_list])
    print("Removing failed SNPs that were masked", file=sys.stderr)
    no_snps_deduped = utilities.deduplicate_list(no_snps_deduped_parital, [s.get_snpid() for s in masked_list])
    if len(no_snps_deduped) > 0:
        print("Writing", len(no_snps_deduped), "failed SNPs to", failedfile, file=sys.stderr)
        with open(failedfile, 'w') as f:
            for s in no_snps_deduped:
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
            method = 'BLAST'
            snp_list, no_snps = blast_based(args, lookup_dict)
        elif args['method'] == 'SAM':
            method = 'SAM'
            raise NotImplementedError("Finding SNPs using a SAM file is not yet implemented")
        masked = filter(lambda s: s.check_masked(), snp_list)
        proper_snps = filter(lambda s: not s.check_masked(), snp_list)
        write_outputs(args, proper_snps, masked, no_snps, method)


if __name__ == '__main__':
    main()
