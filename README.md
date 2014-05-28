# BacterialDB Fetcher

Utilities to fetch biological databases relevant for bacterial genome or metagenome annotation.

**Dependency : GNU parallel**

## UniprotKB

```
$ bash Uniprot/Create-UniprotKB-release.sh -h

Create-UniprotKB-release.sh -d 2014_04
   optional :   [-r]  Download the RDF database
   optional :   [-b]  Create Blast Database
```

## NCBI-Genomes Bacteria

```
$ bash NCBI-Bacterial-Genomes/Create-NCBI-Bacterial-Release.sh -h

Create-NCBI-Bacterial-Release.sh -d 2014_04
   option -b create blast database
   option -c number of core [default 8] (for extraction processing)
   
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

