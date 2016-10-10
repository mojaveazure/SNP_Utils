# SNP_Utils

SNP\_Utils is a Python program that creates a VCF file for a list of SNPs given in an Illumina lookup table. There are three subroutines for SNP\_Utils: `CONFIG`, `BLAST`, and `SAM`

To get basic usage information, simply run the program without any arguments

```
$ ./snp_utils.py
usage: snp_utils.py [-h] {CONFIG,BLAST,SAM} ...

optional arguments:
  -h, --help          show this help message and exit

Subroutine:
  Choose a subroutine

  {CONFIG,BLAST,SAM}  'BLAST' means we are running a BLAST search to find
                      SNPs, 'SAM' means we are using a SAM file, 'CONFIG' is
                      used to configure a BLAST search
```

More detailed information can be found for each subroutine by passing the name of the subroutine and the `-h | --help` flag or by reading below

## Illumina Lookup Table

The Illumina lookup table is a two-column, headerless table that has a SNP ID and contextual sequence in with the SNP in brackets (`[A/B]`) shwogin the two states for the SNP (`A` and `B`). The two columns are *tab*-delimited.

Example

```
SNP_1   ACGTCACGATCGA[A/G]ACGTATGCGAAGTTCGCC
SNP_2   GCTAGACTACCAG[G/T]GTCACGATGCCGTCAGTC
```

## `CONFIG`

The `CONFIG` subroutine is used only when running BLAST within SNP\_Utils. This step is **required** before running BLAST, but not required before parsing a BLAST XML file. Options for `CONFIG` are as follows:

 - Choosing whether the reference is in FASTA or nucleotide BLAST database format and specifying the path to these files
 - E-value threshold
 - Maximum number of hits and hsps
 - Percent identity to keep
 - Whether or not to keep the FASTA file generated from the Illumina lookup table

The configuration file is written in [INI format](https://www.wikiwand.com/en/INI_file)

## `BLAST`

The `BLAST` subroutine is used to run and parse BLAST results to create a VCF file of SNPs from an Illumina lookup table. To run BLASTn within SNP\_Utils, you must configure using the `CONFIG` subroutine. `BLAST` can also rank SNPs; high-ranking SNPs are those that had a low e-value and high bit-score in BLAST. When ranking, only the highest hit for every SNP is kept, ties for highest include both SNPs. To parse a previously-generated BLAST XML file, you do not need to configure BLAST. Options for `BLAST` are as follows:

 - Choosing either a BLAST config or XML file as input
 - Setting the basename for the output
 - Choosing whether or not we rank SNPs

## `SAM`

The `SAM` subroutine is not yet implemented

## Dependencies
SNP\_Utils depends on the following:
 - [Python 3](https://www.python.org/downloads/)
 - [NCBI BLASTn](http://blast.ncbi.nlm.nih.gov/Blast.cgi)
 - [BioPython](http://biopython.org/wiki/Biopython)
 - [BeautifulSoup 4](https://www.crummy.com/software/BeautifulSoup/)
 - [Overload](https://pypi.python.org/pypi/overload)
 - [lxml](http://lxml.de/)

The last four all available on [PyPi](https://pypi.python.org/pypi) and can be downloaded using [pip3](https://pip.pypa.io/en/latest/installing/) (included with Python 3.4 or greater)
