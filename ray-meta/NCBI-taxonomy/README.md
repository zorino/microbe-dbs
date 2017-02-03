
See also sebhtml [Paper-Replication-2012](https://github.com/sebhtml/Paper-Replication-2012)

[Ray](https://github.com/sebhtml/ray) can be utilized to classify k-mers in a taxonomy. To do so,
Ray needs a taxonomy. You can use anything for the taxonomy.
At our center, we are using Greengenes and NCBI.

See these documents for general documentation (in the Ray distribution)
 about graph coloring and taxonomic profiling
features (called Ray Communities):

- Documentation/Taxonomy.txt
- Documentation/BiologicalAbundances.txt

Run this:

	NCBI-taxonomy-Main.sh

Will fetch genome collection from NCBI ftp://ftp.ncbi.nlm.nih.gov/genomes/

And will generate these files:

Genome sequence files :
* NCBI-taxonomy/NCBI-Genomes-Bacteria
* NCBI-taxonomy/NCBI-Genomes-Bacteria_DRAFT
* NCBI-taxonomy/NCBI-Genomes-Viruses
* NCBI-taxonomy/NCBI-Genomes-Plasmids

Taxonomy files :
* NCBI-taxonomy/Genome-to-Taxon.tsv
* NCBI-taxonomy/TreeOfLife-Edges.tsv
* NCBI-taxonomy/Taxon-Names.tsv


Now, you can run Ray as usual (including Ray Méta plugins), but with
additional options to run Ray Communities plugins as well:

Configuration should look like this
(either in the submission file or in the configuration file) :


```
-k 31
-o Ray-Communities

-p SeqA_1.fastq SeqA_2.fastq
-p SeqB_1.fastq SeqB_2.fastq

-search NCBI-taxonomy/NCBI-Genomes-Bacteria
-search NCBI-taxonomy/NCBI-Genomes-Bacteria_DRAFT
-search NCBI-taxonomy/NCBI-Genomes-Viruses
-search NCBI-taxonomy/NCBI-Genomes-Plasmids

-with-taxonomy NCBI-taxonomy/Genome-to-Taxon.tsv
	           NCBI-taxonomy/TreeOfLife-Edges.tsv
	           NCBI-taxonomy/Taxon-Names.tsv

```	

So basically, the whole thing does a distributed de Bruijn graph really
fast (plugins for the distributed storage engine), assembles de novo the
data by distributed graph traversals (Ray Méta; plugin SeedExtender),
colors the graph with the reference genomes provided with the -search
option (Ray Communities, plugin Searcher), and computes taxonomic profiles
using the provided taxonomy (Ray Communities, -with-taxonomy, plugin PhylogenyViewer).


All that stuff is heavily distributed -- each Ray process has 32768 user-space threads
(workers) and you can throw as many Ray processes as you want to.


If you are running Ray on a buggy network (we had problems with Mellanox Infiniband MT26428,
revision a0), you can turn on virtual communications too.

