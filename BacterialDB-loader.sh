#!/bin/bash
# author : maxime déraspe
# email : maxime@deraspe.net
#
# Script to update NCBI Genbank Bacterial Genomes
# fetch those files :
# ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/


if [[ $1 == "" ]];then
    echo "Need a Release Dates as argument (-d) and a DB to download (-g, -p, -f) !"
    echo " BacterialDB-loader.sh -d 2014_04"
    echo "     option -g download ncbi genbank bacteria (ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/)"
    echo "     option -p download ncbi plasmids (ftp://ftp.ncbi.nlm.nih.gov/genomes/Plasmids/)"
    echo "     option -f download ebi phages (ftp://ftp.ebi.ac.uk/pub/databases/fastafiles/embl_genomes/genomes/Phage/)"
    echo "     option -b create blast database (BLAST executables need to be in your path)"
    exit 0
else
    while getopts "d: g f p b h" opt
    do
	case $opt in
            d)
		date=$OPTARG
		;;
	    g)
		genome=1
		;;
	    f)
		phage=1
		;;
	    p)
		plasmid=1
		;;
            b)
		blast=1
		;;
            h)
		echo "Need a Release Dates as argument (-d) and a DB to download (-g, -p, -f) !"
		echo " BacterialDB-loader.sh -d 2014_04"
		echo "     option -g download ncbi genbank bacteria (ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/)"
		echo "     option -p download ncbi plasmids (ftp://ftp.ncbi.nlm.nih.gov/genomes/Plasmids/)"
		echo "     option -f download ebi phages (ftp://ftp.ebi.ac.uk/pub/databases/fastafiles/embl_genomes/genomes/Phage/)"
		echo "     option -b create blast database (BLAST executables need to be in your path)"
		exit 0
		;;
            \?)
		echo "Invalid Options !"
		exit 0
		;;
	esac
    done
fi


if [[ $date == "" ]]
then
    echo "Need a Release Dates as argument (-d) !"
    exit 0
else
    mkdir BacterialDB-$date
    cd BacterialDB-$date
fi


## DOES IT ALL - from scratch
# Fetch Bacteria Genomes
if [[ $genome -eq 1 ]]
then
    echo "#[NCBI Genbank Bacterial Genomes]#"
    echo "Listing NCBI Genbank Bacterial Genomes from ftp"
    wget -q --no-parent --no-remove-listing --spider ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/
    mv .listing Bacteria-NCBI-$date.listing
    nb_of_genomes=`wc -l < Bacteria-NCBI-$date.listing`
    iterator=1
    mkdir Bacteria-NCBI-$date
    cd Bacteria-NCBI-$date
    while read -r i
    do
	# dr-xr-xr-x   5 ftp      anonymous     4096 Mar 10 04:50 Abiotrophia_defectiva
	mkdir $i
	cd $i
	echo "Downloading $i assemblies [$iterator/$nb_of_genomes]"
	iterator=$(($iterator+1))
	wget -q ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/$i/assembly_summary.txt
	while read -r i
	do
	    wget -q -r -nH --cut-dirs=2 --retr-symlinks --no-remove-listing --retr-symlinks --no-remove-listing $i
	done < <(awk -F "\t" '{print $20}' assembly_summary.txt)
	cd ../
    done < <(awk '{print $9}' ../Bacteria-NCBI-$date.listing | sed "s/\\r//g")
fi


## Plasmids ##
if [[ $plasmid -eq 1 ]]
then
    mkdir Plasmids-NCBI-$date
    cd Plasmids-NCBI-$date
    echo " Fetching Bacterial Plasmids Now .."
    wget --no-parent --no-remove-listing --spider ftp://ftp.ncbi.nlm.nih.gov/genomes/Plasmids/fna/
fi

## Phages ##
if [[ $phages -eq 1 ]]
then
    mkdir Phages-EBI-$date
    cd Phages-EBI-$date
    echo " Fetching Bacterial Phages Now .."
    get -r -l1 --no-parent -nd ftp://ftp.ebi.ac.uk/pub/databases/fastafiles/embl_genomes/genomes/Phage/
    wget --no-parent --no-remove-listing --spider ftp://ftp.ncbi.nlm.nih.gov/genomes/Plasmids/fna/
fi


## Create Blast Database ##
if [[ $blast -eq 1 ]]
then
    # export PATH=/is1/commonPrograms/ncbi-blast-2.2.29+-src/c++/BUILD/bin:$PATH
    # Creating BlastDB
    mkdir BlastDB
    cd BlastDB
    # find ../ -name *.fna.gz | parallel -j $core "zcat {} >> all-contigs.fna"
    # makeblastdb -in uniprotKB_$date.fasta -input_type fasta -dbtype prot
    # mv uniprotKB_$date.fasta.* BlastDB
fi
