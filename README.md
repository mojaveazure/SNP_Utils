# SNP_Utils

This depends on:
 - [Python 3](https://www.python.org/downloads/)
 - [NCBI BLASTn](http://blast.ncbi.nlm.nih.gov/Blast.cgi)
 - [BioPython](http://biopython.org/wiki/Biopython)
 - [BeautifulSoup 4](https://www.crummy.com/software/BeautifulSoup/)
 - [Overload](https://pypi.python.org/pypi/overload)
 - [lxml](http://lxml.de/)
 
To get help, utilize the wonderful `-h` option!

```python
./snp_utils.py -h
```

## Lookup Tables

The lookup table is a two-column table that has a SNP ID and contextual sequence in Illumina format (sequence with SNP as `[A/B]` where `A` and `B` are the two states for the SNP). The two columns are *tab*-delimited

Example

```
SNP_1   ACGTCACGATCGA[A/G]ACGTATGCGAAGTTCGCC
```
