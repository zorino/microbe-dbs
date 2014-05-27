# BioDB Fetcher

Utilities to fetch biological databases relevant for biological annotations.


## UniprotKB

	bash Uniprot/Create-UniprotKB-release.sh -h

## NCBI-Genomes Bacteria

	bash NCBI-Bacterial-Genomes/Create-NCBI-Bacterial-Release.sh -h

## Ray Communities for bacterial profiling

This originate from the project [Ray Meta](https://github.com/sebhtml/ray) for the profiling of metagenomes, see also the [Paper in Genome Biology](http://dx.doi.org/doi:10.1186/gb-2012-13-12-r122).

### NCBI-taxonomy


It differs from NCBI-Genomes Bacteria because it fetches from :
ftp://ftp.ncbi.nlm.nih.gov/genomes/* VS ftp://ftp.ncbi.nlm.nih.gov/genbank/genomes/*

*Eventually to be merged together. Genbank subfolder seems more updated*

see the NCBI-taxonomy/README.md

### Gene Ontology (GO)

see the GO/README.md

