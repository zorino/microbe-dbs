# microbe-dbs #

Utilities to fetch biological databases relevant for microbial genomics analysis.

```
$ bash microbe-dbs.sh -h
 microbe-dbs.sh <databases> <output>
    
    databases :
      -> refseq_bacteria
      -> refseq_plasmid
      -> refseq_viral
      -> ncbi_ar
      -> vfdb
      -> ...
```

## Databases ##

### Bacteria ###

#### refseq_bacteria ####

* description: RefSeq: NCBI Reference Sequence Database
* release: ftp://ftp.ncbi.nlm.nih.gov/refseq/release/RELEASE_NUMBER
* url: ftp://ftp.ncbi.nlm.nih.gov/refseq/release/bacteria/
* license: NCBI places no restrictions on the use or distribution. However, submitters may claim copyright of the data they have submitted.

#### refseq_plasmid ####

* description: RefSeq: NCBI Reference Sequence Database
* release: ftp://ftp.ncbi.nlm.nih.gov/refseq/release/RELEASE_NUMBER
* url: ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plasmid/
* license: NCBI places no restrictions on the use or distribution. However, submitters may claim copyright of the data they have submitted.

#### cog ####

* description: Clusters of Orthologous Groups of proteins database
* release: 2014 update
* url: https://www.ncbi.nlm.nih.gov/COG/
* license: NCBI places no restrictions on the use or distribution. However, submitters may claim copyright of the data they have submitted.

#### pub_mlst ####

* description: PubMLST - Public databases for molecular typing and microbial genome diversity.
* release: by date (https://pubmlst.org/news.shtml)
* url: https://pubmlst.org/
* license: https://pubmlst.org/policy.shtml [7]

#### vfdb ####

* description: Virulence Factor Database
* release: by date
* url: http://www.mgc.ac.cn/VFs/download.htm
* license: none reported. Copyrighted to Key Laboratory of Systems Biology of Pathogens, Institue of Pathogen Biology, CAMS&PUMC, Bejing, China.

#### mibig ####

* description: The Minimum Information about a Biosynthetic Gene cluster (MIBiG) database
* release: 1.3 (September 3rd, 2016)
* url: http://mibig.secondarymetabolites.org
* license: Creative Commons Attribution 4.0 International License

#### BacDive ####

* description: The Bacterial Diversity Metadatabase
* release: by date
* url: https://www.bacdive.de/
* license: See Term of use and Copyright (https://bacdive.dsmz.de/about)


#### ncbi_arg ####

* description: NCBI - Bacterial Antimicrobial Resistance Reference Gene Database
* release: by date
* url: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA313047
* license: NCBI places no restrictions on the use or distribution. However, submitters may claim copyright of the data they have submitted.


### Virus ###

#### refseq_viral ####

* description: RefSeq: NCBI Reference Sequence Database
* release: ftp://ftp.ncbi.nlm.nih.gov/refseq/release/RELEASE_NUMBER
* url: ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/
* license: NCBI places no restrictions on the use or distribution. However, submitters may claim copyright of the data they have submitted.


### Other ###

#### ncbi_taxonomy ####

* description: NCBI Taxonomy
* release: by date
* url: ftp://ftp.ncbi.nih.gov/pub/taxonomy/
* license: NCBI places no restrictions on the use or distribution. However, submitters may claim copyright of the data they have submitted.

#### refseq_archaea ####

* description: RefSeq: NCBI Reference Sequence Database
* release: ftp://ftp.ncbi.nlm.nih.gov/refseq/release/RELEASE_NUMBER
* url: ftp://ftp.ncbi.nlm.nih.gov/refseq/release/archaea/
* license: NCBI places no restrictions on the use or distribution. However, submitters may claim copyright of the data they have submitted.

#### refseq_fungi ####

* description: RefSeq: NCBI Reference Sequence Database
* release: ftp://ftp.ncbi.nlm.nih.gov/refseq/release/RELEASE_NUMBER
* url: ftp://ftp.ncbi.nlm.nih.gov/refseq/release/fungi/
* license: NCBI places no restrictions on the use or distribution. However, submitters may claim copyright of the data they have submitted.

#### ebi_idmapping ####

* description: EBI - Uniprot ID mappings with different databases
* release: ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/RELEASE.metalink
* url: ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/
* license: Creative Commons Attribution-NoDerivs 3.0


### Ray meta ###

For Ray Meta DBs creation see ./ray-meta/
