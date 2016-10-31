#!/usr/bin/env python3
"""A script to find SNP positions"""

import sys
if sys.version_info.major is not 3:
    sys.exit("Please use Python 3 for this script")


try:
    from Utilities import filters
    from Utilities import arguments
    from Utilities import utilities
    from Utilities.utilities import INFO
    from Objects import snp
    from Objects import blast
    from Objects import alignment
    from Objects.snp import NotABaseError
    from Objects.snp import NoMatchError
    from Objects.blast import NoSNPError
    from BLAST import runblastn
    from BLAST import configure
except ImportError:
    sys.exit("Please make sure you are in the 'SNP_Utils' directory to load custom modules")


from os.path import basename

try:
    from bs4 import BeautifulSoup
    from bs4 import FeatureNotFound
    from Bio import SeqIO
except ImportError as error:
    sys.exit("Please install " + error.name)


REFERENCE_DEFAULT = basename(sys.argv[0])

#   BLAST-based
def blast_based(args, lookup_dict):
    """Run SNP_Utils using a BLAST-based method"""
    try: # Type checking
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
        bconf = dict() # In case of filtering
    try:
        #   Make soup out of the XML
        blast_soup = BeautifulSoup(blast_xml, 'xml')
    except FeatureNotFound:
        sys.exit("Please install 'lxml' to properly parse the BLAST results")
    no_hits = set() # A set to hold no hits
    snp_list = [] # A list to hold all SNPs found
    hsps = [] # A list of HSPs from the XML file
    ref_gen = blast.get_value(tag=blast_soup, value='BlastOutput_db')
    bconf['database'] = ref_gen
    if not ref_gen:
        ref_gen = REFERENCE_DEFAULT
        bconf = None
    for query in blast_soup.findAll('Iteration'):
        snpid = blast.get_value(tag=query, value='Iteration_query-def')
        #   Ask if no hits were found
        try:
            if blast.get_value(tag=query, value='Iteration_message').capitalize() == blast.NO_HIT_MESSAGE:
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
    if 'rank' in args.keys():
        final_hsps = blast.rank_hsps(hsps=hsps)
    else:
        final_hsps = {h.get_name() : h for h in hsps}
    for hsps in final_hsps.values():
        for hsp in hsps:
            try:
                snp_list.append(hsp.get_snp())
            except NoSNPError:
                print('No SNP for', hsp, file=sys.stderr)
                no_hits.add(hsp.get_name())
    hit_snps = {s.get_snpid() for s in snp_list}
    no_hits -= hit_snps
    return(snp_list, no_hits, ref_gen, bconf)


#   Alignment-blast_based
def alignment_based(args, lookup_dict):
    """Run SNP_Utils using a SAM file"""
    try: # Type checking
        assert isinstance(args, dict)
        assert isinstance(lookup_dict, dict)
    except AssertionError:
        raise TypeError
    alignment_dict = dict()
    no_hits = set()
    snp_list = list()
    try:
        #   Read in the reference FASTA file
        print("Reading in reference FASTA", args['reference'], "This might take a while...", file=sys.stderr)
        reference = SeqIO.to_dict(SeqIO.parse(args['reference'], 'fasta'))
        #   Read in the SAM alignment
        print("Reading in SAM file", args['samfile'], file=sys.stderr)
        with open(args['samfile'], 'r') as s:
            for line in s:
                if line.startswith('@'):
                    continue
                else:
                    a = alignment.Alignment(line)
                    if a.check_flag():
                        alignment_dict[a.get_name()] = a
    except FileNotFoundError as error:
        sys.exit("Failed to find " + error.filename)
    for snpid in lookup_dict:
        if snpid in alignment_dict.keys() and alignment_dict[snpid].get_contig() is not '*':
            s = snp.SNP(
                lookup=lookup_dict[snpid],
                alignment=alignment_dict[snpid],
                reference=reference
            )
            snp_list.append(s)
        else:
            print("No map for", snpid)
            no_hits.add(snpid)
    if len(snp_list) < 1:
        sys.exit("Failed to find any SNPs in " + args['samfile'] + " that were listed in " + args['lookup'])
    no_hits -= {s.get_snpid() for s in snp_list}
    return (snp_list, no_hits)


#   Write the output files
def write_outputs(args, snp_filter, masked_filter, no_snps, ref_gen, vcf_info):
    """Write the output files"""
    try:
        assert isinstance(args, dict)
        assert isinstance(snp_filter, (filter, list, tuple))
        assert isinstance(masked_filter, (filter, list, tuple))
        assert isinstance(no_snps, list)
        assert isinstance(ref_gen, str)
        assert isinstance(vcf_info, (list, tuple))
    except AssertionError:
        raise TypeError
    header = utilities.vcf_header(ref_gen, sorted(vcf_info))
    outfile = args['outname'] + '.vcf'
    maskedfile = args['outname'] + '_masked.vcf'
    failedfile = args['outname'] + '_failed.log'
    snp_list = list(snp_filter)
    snp_list.sort(key=lambda s: (s.get_chrom(), s.get_position()))
    if len(snp_list) < 1:
        print("Failed to find any SNPs", file=sys.stderr)
    else:
        print("Writing", len(snp_list), "SNPs to", outfile, file=sys.stderr)
        with open(outfile, 'w') as o:
            o.write(header)
            o.write('\n')
            for s in snp_list:
                o.write(s.format_vcf())
                o.write('\n')
    print("Removing masked SNPs that were actually found", file=sys.stderr)
    masked_list = list(masked_filter)
    if len(masked_list) > 0:
        print("Writing", len(masked_list), "masked SNPs to", maskedfile, file=sys.stderr)
        with open(maskedfile, 'w') as m:
            m.write(header)
            m.write('\n')
            for s in masked_list:
                m.write(s.format_vcf())
                m.write('\n')
    if len(no_snps) > 0:
        print("Writing", len(no_snps), "failed SNPs to", failedfile, file=sys.stderr)
        with open(failedfile, 'w') as f:
            for s in sorted(no_snps):
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
        arguments.validate_filters(args=args, parser=parser) # Validate filtering options
        #   Two holding collections
        lookup_dict = {}
        vcf_info = []
        #   Read in the lookup table
        with open(args['lookup'], 'r') as lk:
            print("Reading in lookup table", args['lookup'], file=sys.stderr)
            for line in lk:
                split = line.strip().split()
                l = snp.Lookup(split[0], split[1])
                lookup_dict[l.get_snpid()] = l
        #   Find SNPs given our method of BLAST or SAM
        if args['method'] == 'BLAST':
            method = 'BLAST'
            snp_list, no_snps, ref_gen, bconf = blast_based(args, lookup_dict)
            vcf_info.append(INFO(infoid='B', number=0, infotype='Flag', description='Variant Calculated from BLAST'))
        elif args['method'] == 'SAM':
            method = 'SAM'
            ref_gen = args['reference']
            bconf = None
            snp_list, no_snps = alignment_based(args, lookup_dict)
            vcf_info.append(INFO(infoid='S', number=0, infotype='Flag', description='Variant Calculated from SAM'))
        #   Start our filtering
        masked = filter(lambda s: s.check_masked(), snp_list)
        proper_snps = tuple(filter(lambda s: not s.check_masked(), snp_list))
        if 'map' in args and args['bychrom']:
            proper_snps, vcf_info = filters.chrom_filter(args=args, vcfinfo=vcf_info, propersnps=proper_snps)
        if 'threshold' in args.keys():
            proper_snps, vcf_info = filters.distance_filter(args=args, vcfinfo=vcf_info, propersnps=proper_snps)
        if 'map' in args and args['bydistance'] and ref_gen is not REFERENCE_DEFAULT:
            proper_snps, vcf_info = filters.map_filter(args=args, vcfinfo=vcf_info, propersnps=proper_snps, refgen=ref_gen, method=method, bconf=bconf)
        elif ref_gen is REFERENCE_DEFAULT:
            print("Could not determine reference genome, not filtering by genetic map distance", file=sys.stderr)
        write_outputs(args, proper_snps, masked, list(no_snps), ref_gen, vcf_info)


if __name__ == '__main__':
    main()
