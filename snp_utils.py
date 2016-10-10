#!/usr/bin/env python3
"""A script to find SNP positions"""

import sys
if sys.version_info.major is not 3:
    sys.exit("Please use Python 3 for this script")


try:
    from Utilities import arguments
    from Utilities import utilities
    from Objects import snp
    from Objects import blast
    from Objects.snp import NotABaseError
    from Objects.snp import NoMatchError
    from Objects.blast import NoSNPError
    from BLAST import runblastn
    from BLAST import configure
except ImportError:
    sys.exit("Please make sure you are in the 'SNP_Utils' directory to load custom modules")


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
        raise TypeError
    #   If we weren't provided a BLAST XML file
    if 'xml' not in args.keys():
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
    no_hits = set() # A set to hold no hits
    snp_list = []
    hsps = [] # A list of HSPs from the XML file
    for query in blast_soup.findAll('Iteration'):
        snpid = blast.get_value(tag=query, value='Iteration_query-def')
        #   Ask if no hits were found
        try:
            if blast.get_value(tag=query, value='Iteration_message').capitalize == blast.NO_HIT_MESSAGE:
                print('No hit for', snpid, file=sys.stderr)
                no_hits.add(snpid)
                continue
        except AttributeError:
            pass
        #   For every hit found
        for hit in query.findAll('Hit'):
            hit_num = blast.get_value(tag=hit, value='Hit_num')
            this_hsps = blast.parse_hit(snpid=snpid, hit=hit)
            try:
                hsps += this_hsps
            except TypeError:
                print('No HSPs for', snpid, 'hit number:', hit_num, file=sys.stderr)
                no_hits.add(snpid)
                continue
    for hsp in hsps:
        snpid = hsp.get_name()
        lookup = lookup_dict[snpid]
        try:
            hsp.add_snp(lookup=lookup)
        except (NoSNPError, NotABaseError, NoMatchError):
            print('No SNP for', hsp, file=sys.stderr)
            no_hits.add(snpid)
            hsps.remove(hsp)
    #   Close the XML file
    blast_xml.close()
    #   Rank, if asked for
    if args['rank']:
        final_hsps = blast.rank_hsps(hsps=hsps)
    else:
        final_hsps = {h.get_name() : h for h in hsps}
    hit_snps = set(final_hsps.keys())
    for hsps in final_hsps.values():
        for hsp in hsps:
            try:
                snp_list.append(hsp.get_snp())
            except NoSNPError:
                print('No SNP for', hsp, file=sys.stderr)
                no_hits.add(hsp.get_name())
    no_hits -= hit_snps
    return(snp_list, no_hits)


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
        raise TypeError
    header = '##fileformat=VCFv4.2\n##INFO<ID=s,Number=1,Type=Flag,Description="Variant is calculated from %s">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' % method
    outfile = args['outname'] + '.vcf'
    maskedfile = args['outname'] + '_masked.vcf'
    failedfile = args['outname'] + '_failed.log'
    snp_list = list(snp_filter)
    snp_list.sort(key=lambda s : (s.get_chrom(), s.get_position()))
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
    args = {key : value for key, value in vars(parser.parse_args()).items() if value is not None}
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
        write_outputs(args, proper_snps, masked, list(no_snps), method)


if __name__ == '__main__':
    main()
