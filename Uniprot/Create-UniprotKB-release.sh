#!/bin/bash
#
# Script to update uniprot db
# fetch those files :
# ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/ ..
#    knowledgebase/complete/uniprot_trembl.fasta.gz
#    knowledgebase/complete/uniprot_sprot.fasta.gz
#    

if [[ $1 == "" ]];then
    echo "Need a Release Dates as argument (-d) !"
    echo " Create-UniprotKB-release.sh -d 2014_04 "
    echo "      optional :   [-r]  Download the RDF database"
    echo "      optional :   [-b]  Create Blast Database"
    exit 0
else
    while getopts "d: r b h" opt
    do
	case $opt in
            d)
		date=$OPTARG
		;;
            r)
		rdf=1
		;;
            b)
		blast=1
		;;
            h)
		echo " Create-UniprotKB-release.sh -d 2014_04 "
		echo "      optional :   [-r]  Download the RDF database"
		echo "      optional :   [-b]  Create Blast Database"
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
    mkdir release-$date
    cd release-$date
fi


## UniprotKB ##
mkdir UniprotKB
cd UniprotKB
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz
cd ../

## IDmapping ##
mkdir IDmapping
cd IDmapping
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz
cd ../

## Gunzip all files ##
gunzip UniprotKB/uniprot_sprot.fasta.gz &
gunzip UniprotKB/uniprot_trembl.fasta.gz &
gunzip IDmapping/idmapping.dat.gz
cat UniprotKB/uniprot_sprot.fasta > uniprotKB_$date.fasta
cat UniprotKB/uniprot_trembl.fasta >> uniprotKB_$date.fasta

## Create Blast Database ##
if [[ $blast -eq 1 ]]
then
    export PATH=/is1/commonPrograms/ncbi-blast-2.2.29+-src/c++/BUILD/bin:$PATH
    mkdir BlastDB
    makeblastdb -in uniprotKB_$date.fasta -input_type fasta -dbtype prot
    mv uniprotKB_$date.fasta.* BlastDB
fi

## RDF ##
if [[ $rdf -eq 1 ]]
then
    mkdir RDF
    cd RDF/
    wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/rdf/uniprot.rdf.gz
    wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/rdf/core.owl
    cd ../
fi
