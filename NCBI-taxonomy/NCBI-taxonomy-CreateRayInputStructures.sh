#!/bin/bash
# Author: Sébastien Boisvert
# Edited by : Maxime Déraspe
# Date: 2014-05-27
#
# Script to Fetch NCBI Taxonomy and Genomes
# and create the files for RayCommunities
# 
# Datasets 
# ftp://ftp.ncbi.nih.gov/genomes/
# 	Bacteria
# 	Bacteria_DRAFT
# 	Plasmids
# 	
# ftp://ftp.ncbi.nih.gov/refseq/release/
# 	viral
#


program=$0
OutputDirectory=$1
waitingSeconds=1

curDir=`pwd`
programDir=$(dirname $program)

export PATH=$PATH:$curDir:$programDir

curDate=$(date +%Y-%m-%d)

if test "$OutputDirectory" = ""
then
	echo "Warning: no output directory was provided, will use NCBI-taxonomy_$curDate"
	OutputDirectory=NCBI-taxonomy_$curDate
fi

# From: http://stackoverflow.com/questions/592620/check-if-a-program-exists-from-a-bash-script
command -v wget >/dev/null 2>&1 || { echo >&2 "Error, needs wget but it's not installed.  Aborting."; exit 1; }

echo "Welcome to Ray technolologies !"
echo "How are you today ?"
echo "Good!"
echo ""
echo "This is the easy-to-use assistant for generating NCBI taxonomy files to use with Ray Communities"


echo "Bioinformatics operations will begin in $waitingSeconds seconds..."

sleep $waitingSeconds

echo "OutputDirectory is $OutputDirectory"


echo "We have wget !"



if test -d $OutputDirectory
then
	echo "Warning: directory $OutputDirectory already exists."
else
	mkdir $OutputDirectory
fi

cd $OutputDirectory

if test ! -d ftp.ncbi.nih.gov
then
	mkdir ftp.ncbi.nih.gov
	cd ftp.ncbi.nih.gov

        if test ! -f gi_taxid_nucl.dmp.gz
        then
    	    echo "Downloading gi_taxid_nucl.dmp.gz, please wait."
	    wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_nucl.dmp.gz
	    echo "Done."
        fi

        if test ! -f taxdump.tar.gz
        then
    	    echo "Downloading taxdump.tar.gz, please wait."
	    wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
	    echo "Done."
        fi

        if test ! -f bacteria.fna.tar.gz
        then
	    echo "Downloading all.fna.tar.gz, please wait."
	    wget -O bacteria.fna.tar.gz ftp://ftp.ncbi.nih.gov/genomes/Bacteria/all.fna.tar.gz
	    echo "Done."
        fi

        if test ! -f plasmids.fna.tar.gz
        then
    	    echo "Downloading plasmids.all.fna.tar.gz, please wait."
	    wget -O plasmids.fna.tar.gz ftp://ftp.ncbi.nih.gov/genomes/Plasmids/plasmids.all.fna.tar.gz
	    echo "Done."
        fi

        if test ! -f viruses.fna.tar.gz
        then
    	    echo "Downloading viruses.all.fna.tar.gz, please wait."
	    wget -O viruses.fna.tar.gz ftp://ftp.ncbi.nih.gov/genomes/Viruses/all.fna.tar.gz
	    echo "Done."
        fi


        if test ! -f viruses.fna.tar.gz
        then
    	    echo "Downloading RefSeq of Viruses, please wait."
            wget -r -A "*fna.gz" ftp://ftp.ncbi.nih.gov/refseq/release/viral/
	    echo "Done."
        fi


        if test ! -d Bacteria_DRAFT
        then
    	    mkdir Bacteria_DRAFT
            cd Bacteria_DRAFT
	    echo "Downloading all.fna.tar.gz, please wait."
    	    wget -r -nH --cut-dirs=2 --no-parent ftp://ftp.ncbi.nih.gov/genomes/Bacteria_DRAFT/
    	    echo "Done."
            cd ..
        fi


        cd ..
else
	echo "$OutputDirectory/ftp.ncbi.nih.gov already exists, skipping."
fi


if test ! -d uncompressed
then
	mkdir uncompressed
	cd uncompressed

	echo "Decompressing taxdump.tar.gz, please wait."
	mkdir taxdump
	cd taxdump
	cat ../../ftp.ncbi.nih.gov/taxdump.tar.gz|gunzip|tar -x
	cd ..
	echo "Done."

	echo "Decompressing bacteria.fna.tar.gz, please wait."
	mkdir bacteria.all.fna
	cd bacteria.all.fna
	cat ../../ftp.ncbi.nih.gov/bacteria.fna.tar.gz|gunzip|tar -x
	cd ..
	echo "Done."

	echo "Decompressing plasmids.fna.tar.gz, please wait."
	mkdir plasmids.all.fna
	cd plasmids.all.fna
	cat ../../ftp.ncbi.nih.gov/plasmids.fna.tar.gz|gunzip|tar -x
	cd ..
	echo "Done."

	echo "Decompressing viruses.fna.tar.gz, please wait."
	mkdir viruses.all.fna
	cd viruses.all.fna
	cat ../../ftp.ncbi.nih.gov/viruses.fna.tar.gz|gunzip|tar -x
	cd ..
	echo "Done."

	echo "Decompressing viruses.fna.tar.gz, please wait."
	mkdir viruses.refseq.all.fna
	cd viruses.refseq.all.fna
	cat ../../ftp.ncbi.nih.gov/refseq/release/viral/*.fna.gz|gunzip>sequences.fna
	cd ..
	echo "Done."


        echo "Decompressing Bacteria_DRAFT/*, please wait."
        mkdir bacteria_draft.all.fna
        cd bacteria_draft.all.fna
        for i in $(ls ../../ftp.ncbi.nih.gov/Bacteria_DRAFT/)
        do
            	mkdir $i
                cd $i
                rm -f *fna
                if test -f ../../../ftp.ncbi.nih.gov/Bacteria_DRAFT/$i/*.contig.fna.tgz
                then
                    tar xzf ../../../ftp.ncbi.nih.gov/Bacteria_DRAFT/$i/*.contig.fna.tgz
                elif test -f ../../../ftp.ncbi.nih.gov/Bacteria_DRAFT/$i/*.scaffold.fna.tgz
                then
                    tar xzf ../../../ftp.ncbi.nih.gov/Bacteria_DRAFT/$i/*.scaffold.fna.tgz
                else
                    echo "no fna file "
                fi
                cat *.fna > ../$i.fasta    
                cd ..
                rm -rf $i 
        done
        cd ..
        echo "Done."


	cd ..
fi


if test ! -f Genome-to-Taxon.tsv
then
	echo "Creating $OutputDirectory/Genome-to-Taxon.tsv, please wait."
	cat ftp.ncbi.nih.gov/gi_taxid_nucl.dmp.gz|gunzip > Genome-to-Taxon.tsv
	echo "Done."
fi

if test ! -f TreeOfLife-Edges.tsv
then
	echo "Creating $OutputDirectory/TreeOfLife-Edges.tsv, please wait."
	cat uncompressed/taxdump/nodes.dmp|awk '{print $3"\t"$1}' > TreeOfLife-Edges.tsv
	echo "Done."
fi


if test ! -f Taxon-Names.tsv
then
	echo "Creating $OutputDirectory/Taxon-Names.tsv, please wait."
	NCBI-taxonomy-Create-Taxon-Names.py uncompressed/taxdump/nodes.dmp uncompressed/taxdump/names.dmp Taxon-Names.tsv
	echo "Done."
fi

if test ! -d NCBI-Genomes-Bacteria
then
	echo "Creating $OutputDirectory/NCBI-Genomes-Bacteria, please wait."

	mkdir NCBI-Genomes-Bacteria
	cd NCBI-Genomes-Bacteria

	for i in $(ls ../uncompressed/bacteria.all.fna)
	do
		name=$(echo $i|sed 's/_uid/ /g'|awk '{print $1}')
	
		cat ../uncompressed/bacteria.all.fna/$i/*.fna > $name".fasta"
	done

	echo "Done."
	
	cd ..
fi

if test ! -d NCBI-Genomes-Plasmids
then
	echo "Creating $OutputDirectory/NCBI-Genomes-Plasmids, please wait."

	mkdir NCBI-Genomes-Plasmids
	cd NCBI-Genomes-Plasmids

        # plasmids.tar.gz has /am/ftp-genomes/Plasmids/fna/ path befores sequence files
	for i in $(ls ../uncompressed/plasmids.all.fna/am/ftp-genomes/Plasmids/fna/)
	do
            	header=`head -1 ../uncompressed/plasmids.all.fna/am/ftp-genomes/Plasmids/fna/$i`
	        name=`echo $header | \
        	    sed "{ 
		        s/>gi|.*|.*|.*| //g
        		s/, complete .*//g
	        	s/ /_/g
			s/\//-/g
	        	}"`
                c=1    
                while [ -f $name-$c.fasta ]
                do
                    let c=$(($c+1))
                done
                cat ../uncompressed/plasmids.all.fna/am/ftp-genomes/Plasmids/fna/$i > "$name-$c.fasta"
	done

	echo "Done."
	
	cd ..
fi


if test ! -d NCBI-Genomes-Viruses
then
	echo "Creating $OutputDirectory/NCBI-Genomes-Viruses, please wait."

	mkdir NCBI-Genomes-Viruses
	cd NCBI-Genomes-Viruses

	for i in $(ls ../uncompressed/viruses.all.fna)
	do
		name=$(echo $i|sed 's/_uid/ /g'|awk '{print $1}')
	
		cat ../uncompressed/viruses.all.fna/$i/*.fna > $name".fasta"
	done

	echo "Done."
	
	cd ..
fi

if test ! -d NCBI-Genomes-Bacteria_DRAFT
then
	echo "Creating $OutputDirectory/NCBI-Genomes-Bacteria_DRAFT, please wait."

	mkdir NCBI-Genomes-Bacteria_DRAFT
	cd NCBI-Genomes-Bacteria_DRAFT
	cp ../uncompressed/bacteria_draft.all.fna/*.fasta .
	echo "Done."
	cd ..
fi


echo ""
echo "Finished, checking files"
echo "Removing uncompressed data ... Please wait"

rm -rf uncompressed/

echo ""
echo "We thank you for your patience !"
echo ""
echo "NCBI taxonomy files are ready to be used with Ray Communities !"
echo "If you need more information, you can read these documents from the ray distribution:"
echo "- Documentation/NCBI-Taxonomy.txt (information about the NCBI taxonomy and Ray Communities)"
echo "- Documentation/Taxonomy.txt (information about the general taxonomy features in Ray Communities)"
echo "- Documentation/BiologicalAbundances.txt (information about biological profiling in Ray Communities)"
echo "- MANUAL_PAGE.txt (full list of options of Ray)"

echo "If you need support, send a email to denovoassembler-users@lists.sf.net"

echo ""
echo "Thank you for choosing Ray for your research."
echo "Happy open assembly and profiling to you !"

