# microbe-dbs #

Utilities to fetch biological databases relevant for microbial genomics analysis.

Some database require python to download and at least 2 libraries (biopython and requests).

You can install the environment with conda from the environment.yml file.


```shell
conda env create -f environment.yml
conda activate microbe-dbs
```

Use the bash CLI to list databases.

```shell
./microbe-dbs -h
 microbe-dbs <databases> <output>
    
    databases :
      -> uniprotkb_bacteria
      -> refseq_bacteria
      -> refseq_plasmid
      -> refseq_viral
      -> ncbi_arg
      -> vfdb
      -> ...
```

## Databases ##

### Bacteria ###

#### uniprotkb_bacteria ####

* description: UniprotKB (SwissProt/TrEMBL) - Bacteria
* release: ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/reldate.txt
* url: ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/
* license: Creative Commons Attribution (CC BY 4.0)

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

#### isfinder ####

* description: ISFinder Database
* release: by date
* url: https://github.com/thanhleviet/ISfinder-sequences
* license: none reported. Citation needed.

#### mibig ####

* description: The Minimum Information about a Biosynthetic Gene cluster (MIBiG) database
* release: Version 1.4 (August 6th, 2018)
* url: http://mibig.secondarymetabolites.org
* license: Creative Commons Attribution 4.0 International License

#### bacdive ####

* description: The Bacterial Diversity Metadatabase
* release: by date
* url: https://www.bacdive.de/
* license: See Term of use and Copyright (https://bacdive.dsmz.de/about)

#### resfinder_db ####

* description: ResFinder DB - ResFinder identifies acquired antimicrobial resistance genes and/or chromosomal mutations
* release: last git commit date
* url: https://bitbucket.org/genomicepidemiology/resfinder_db/
* license: Apache License, Version 2.0 (the "License")

#### ncbi_arg ####

* description: NCBI - Bacterial Antimicrobial Resistance Reference Gene Database
* release: by date
* url: https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/Data/
* license: Public Domain (see https://raw.githubusercontent.com/ncbi/amr/master/LICENSE)

#### ebi_mgnify_hg ####

* description: EBI MGnify Human Gut Proteins Catalogue 90 (high quality)
* release: by date
* url: http://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/2019_09/uhgp_catalogue/uhgp-90/
* license: See https://www.ebi.ac.uk/metagenomics/about


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

