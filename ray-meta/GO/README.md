
This is for building input files for gene ontology profiling
with Ray Ontologies.

	GO-Main.sh [output directory]

This will generate these files :
* EMBL_CDS_Sequences/
* 000-Annotations.txt
* 000-Ontologies.txt
* 000-Sequences.fasta
* Annotations.txt
* Associations.txt
* cds.fasta
* gene_association.goa_uniprot
* gene_ontology_ext.obo
* idmapping.dat
* OntologyTerms.txt
* Sequences-With-Annotations.fasta

Now, you can run Ray as usual (including Ray Méta plugins), but with
additional options to run Ray Communities plugins as well:

Configuration should look like this
(either in the submission file or in the configuration file) :


```
-k 31
-o Ray-Communities

-p SeqA_1.fastq SeqA_2.fastq
-p SeqB_1.fastq SeqB_2.fastq

-search GO/EMBL_CDS_Sequences

-gene-ontology GO/OntologyTerms.txt
	           GO/Annotations.txt

```	
