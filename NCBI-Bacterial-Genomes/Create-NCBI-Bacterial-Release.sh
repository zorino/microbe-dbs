#!/bin/bash
#
# Script to update NCBI Genbank Bacterial Genomes
# fetch those files :
# ftp://ftp.ncbi.nlm.nih.gov/genbank/genomes/
#    Bacteria_DRAFT/
#    Bacteria/
#    
# Dependency : GNU parallel


if [[ $1 == "" ]];then
    echo "Need a Release Dates as argument (-d) !"
    echo " Create-NCBI-Bacterial-Release.sh -d 2014_04"
    echo "     option -b create blast database"
    echo "     option -c number of core [default 8] (for extraction processing)"
    exit 0
else
    while getopts "d:c: r b h" opt
    do
	case $opt in
            d)
		date=$OPTARG
		;;
            c)
		core=$OPTARG
		;;
            r)
		rdf=1
		;;
            b)
		blast=1
		;;
            h)
		echo " Create-NCBI-Bacterial-Release.sh -d 2014_04"
		echo "     option -b create blast database"
		echo "     option -c number of core [default 8] (for extraction processing)"
		exit 0
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
    mkdir NCBI-Genbank-Bacteria-$date
    cd NCBI-Genbank-Bacteria-$date
fi

if [[ $core == "" ]]
then
    core=8
fi


## Bacteria_DRAFT ##

echo " Fetching Bacteria_DRAFT Genomes from NCBI ftp .."
echo " please wait .."
wget -A *.gbk* -A *.fna* -A *.faa* -r -nH --cut-dirs=2 --no-parent ftp://ftp.ncbi.nlm.nih.gov/genbank/genomes/Bacteria_DRAFT/
echo " Finished Bacteria_DRAFT fetching"

echo ""
echo " Processing the extraction of compressed file (.tgz) .."
echo " Using GNU parallel with -j $core"

cd Bacteria_DRAFT/
find . -name *.tgz | parallel -j $core "tar xvf {} -C {//}"
cd ../

echo " Bacteria_DRAFT Extraction Done !"
echo ""

echo " Fetching Bacteria Complete Genomes Now .."

## Bacteria ##
wget -A *.gbk* -A *.fna* -A *.faa* -r -nH --cut-dirs=2 --no-parent ftp://ftp.ncbi.nlm.nih.gov/genbank/genomes/Bacteria/


echo " Finished Bacteria fetching"

# mkdir Bacteria-All/
# cd Bacteria-All/

## Create Blast Database ##
if [[ $blast -eq 1 ]]
then
    export PATH=/is1/commonPrograms/ncbi-blast-2.2.29+-src/c++/BUILD/bin:$PATH
    mkdir BlastDB
    makeblastdb -in uniprotKB_$date.fasta -input_type fasta -dbtype prot
    mv uniprotKB_$date.fasta.* BlastDB
fi
