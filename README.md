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
      -> 
```

## Databases ##

### Bacteria ###

#### refseq_bacteria ####

* url: ftp://ftp.ncbi.nlm.nih.gov/refseq/release/bacteria/
* license: NCBI places no restrictions on the use or distribution. However, submitters may claim copyright of the data they have submitted.

#### refseq_plasmid ####
* url: ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plasmid/
* license: NCBI places no restrictions on the use or distribution. However, submitters may claim copyright of the data they have submitted.

#### COG ####

* url: https://www.ncbi.nlm.nih.gov/COG/
* license: NCBI places no restrictions on the use or distribution. However, submitters may claim copyright of the data they have submitted.

#### vfdb ####

* url: http://www.mgc.ac.cn/VFs/download.htm
* license: none reported. Copyrighted to Key Laboratory of Systems Biology of Pathogens, Institue of Pathogen Biology, CAMS&PUMC, Bejing, China.

#### mibig ####

* url: http://mibig.secondarymetabolites.org
* license: Creative Commons Attribution 4.0 International License


### Virus ###

#### refseq_viral ####

* url: ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/
* license: NCBI places no restrictions on the use or distribution. However, submitters may claim copyright of the data they have submitted.


### Other ###

#### refseq_archaea ####

* url: ftp://ftp.ncbi.nlm.nih.gov/refseq/release/archaea/
* license: NCBI places no restrictions on the use or distribution. However, submitters may claim copyright of the data they have submitted.



### Ray meta ###

For Ray Meta DBs creation see ./ray-meta/
