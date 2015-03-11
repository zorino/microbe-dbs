# BacterialDB Fetcher

Utilities to fetch biological databases relevant for bacterial genomics analysis.

```
$ bash BacterialDB-loader.sh -h
 BacterialDB-loader.sh -d 2015_03
     option -g download ncbi genbank bacteria (ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/)
     option -p download ncbi plasmids (ftp://ftp.ncbi.nlm.nih.gov/genomes/Plasmids/)
     option -f download ebi phages (ftp://ftp.ebi.ac.uk/pub/databases/fastafiles/embl_genomes/genomes/Phage/)
     option -b create blast database (BLAST executables need to be in your path)
```

## Ray Communities for bacterial metagenomic profiling

This originate from the project [Ray Meta](https://github.com/sebhtml/ray) for the profiling of metagenomes, see also the [Paper in Genome Biology](http://dx.doi.org/doi:10.1186/gb-2012-13-12-r122).

### NCBI-taxonomy

It differs from NCBI-Genomes Bacteria because it fetches from :

ftp://ftp.ncbi.nlm.nih.gov/genomes/* VS ftp://ftp.ncbi.nlm.nih.gov/genbank/genomes/*

*Eventually to be merged together. Genbank subfolder seems more updated*

see the NCBI-taxonomy/README.md

### Gene Ontology (GO)

see the GO/README.md

