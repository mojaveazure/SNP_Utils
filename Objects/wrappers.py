#!/usr/bin/env python3

"""Wrapper for various NCBI Blast+ applications not provided by BioPython"""

import sys
if sys.version_info.major is not 3:
    (sys.exit("Please use Python 3 for this modul, " + __name__))


from collections import OrderedDict

try:
    from Bio.Application import _Option, AbstractCommandline
except ImportError as error:
    sys.exit("Failed to find " + error.name)


#   Blastdbcmd partial wrapper
class NcbiblastdbcmdCommandline(AbstractCommandline):
    """Command line wrapper for using blastdbcmd"""

    #   Possible output formats
    _OUTFMTS = OrderedDict(
        [('%f', "Sequence in FASTA format"),
         ('%s', "Sequence data without defline"),
         ('%a', "Accession"),
         ('%g', "GI"),
         ('%o', "Original ID (OID)"),
         ('%t', "Sequence title"),
         ('%l', "Sequence length"),
         ('%h', "Sequence hash value"),
         ('%T', "Taxid"),
         ('%X', "Leaf-node taxids"),
         ('%e', "Memebership integer"),
         ('%L', "Common taxonomic name"),
         ('%C', "Common taxonomic name for leaf-node taxids"),
         ('%S', "Scientific name"),
         ('%N', "Scientific name for leaf-node taxids"),
         ('%B', "BLAST name"),
         ('%K', "Taxonomic super kingdom"),
         ('%P', "PIG"),
         ('%m', "Sequence masking data")]
    )

    def __init__(self, cmd='blastdbcmd', **kwargs):
        if cmd is not 'blastdbcmd':
            raise ValueError("This wrapper is specific for 'blastdbcmd'")
        outfmtmsg = "Output format where:\n"
        for item in self._OUTFMTS.items():
            outfmtmsg += '\t' + ' means '.join(item) + '\n'
        self.parameters = [
            _Option(
                names=['-db', 'db'],
                description="BLAST database name",
                is_required=True,
                filename=True,
                equate=False
            ),
            _Option(
                names=['-dbtype', 'dbtype'],
                description="Molecule type stored in BLAST database",
                checker_function=lambda value: value in ['guess', 'nucl', 'prot'],
                equate=False
            ),
            _Option(
                names=['-entry', 'entry'],
                description="Search string for sequence identifier",
                equate=False
            ),
            _Option(
                names=['-outfmt', 'outfmt'],
                description=outfmtmsg,
                checker_function=lambda value: value in self._OUTFMTS,
                equate=False
            )
        ]
        super(self.__class__, self).__init__(cmd, **kwargs)
