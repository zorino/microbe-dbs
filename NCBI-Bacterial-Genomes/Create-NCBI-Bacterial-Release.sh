#!/bin/bash
# author : maxime déraspe
# email : maxime@deraspe.net
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

wget --no-parent --no-remove-listing --spider ftp://ftp.ncbi.nlm.nih.gov/genbank/genomes/Bacteria_DRAFT/
mv .listing Bacteria_DRAFT.listing

mkdir Bacteria_DRAFT
cd Bacteria_DRAFT

awk '{print $9}' ../Bacteria_DRAFT.listing | sed "s#\r##g" | parallel -j 1 \
"wget -a ../Bacteria_DRAFT.dwl.log -A *.gbk* -A *.fna* -A *.faa* -r -nH --cut-dirs=3 --no-parent ftp://ftp.ncbi.nlm.nih.gov/genbank/genomes/Bacteria_DRAFT/{/}"

# wget -A *.gbk* -A *.fna* -A *.faa* -r -nH --cut-dirs=2 --no-parent ftp://ftp.ncbi.nlm.nih.gov/genbank/genomes/Bacteria_DRAFT/

echo " Finished Bacteria_DRAFT fetching"

echo ""
echo " Processing the extraction of compressed file (.tgz) .."
echo " Using GNU parallel with -j $core"


find . -name *.tgz | parallel -j $core "tar xvf {} -C {//}"
cd ../

echo " Bacteria_DRAFT Extraction Done !"
echo ""

## Bacteria ##

echo " Fetching Bacteria Complete Genomes Now .."

wget --no-parent --no-remove-listing --spider ftp://ftp.ncbi.nlm.nih.gov/genbank/genomes/Bacteria/
mv .listing Bacteria.listing

mkdir Bacteria
cd Bacteria

awk '{print $9}' ../Bacteria.listing | sed "s#\r##g" | parallel -j 1 \
"wget -a ../Bacteria.dwl.log -A *.gbk* -A *.fna* -A *.faa* -r -nH --cut-dirs=3 --no-parent ftp://ftp.ncbi.nlm.nih.gov/genbank/genomes/Bacteria/{/}"

# wget -A *.gbk* -A *.fna* -A *.faa* -r -nH --cut-dirs=2 --no-parent ftp://ftp.ncbi.nlm.nih.gov/genbank/genomes/Bacteria/

cd ../

echo " Finished Bacteria fetching"

## Plasmids ##
echo " Fetching Bacterial Plasmids Now .."

wget --no-parent --no-remove-listing --spider ftp://ftp.ncbi.nlm.nih.gov/genomes/Plasmids/fna/
mv .listing Bacteria_Plasmids.listing

mkdir Bacteria_Plasmids
cd Bacteria_Plasmids

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Plasmids/plasmids.all.gbk.tar.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Plasmids/plasmids.all.fna.tar.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Plasmids/plasmids.all.faa.tar.gz

tar zxvf plasmids.all.gbk.tar.gz &
tar zxvf plasmids.all.fna.tar.gz &
tar zxvf plasmids.all.faa.tar.gz &

wait

for i in $(find . -name gbk -type d); do mv $i ../; done
for i in $(find . -name fna -type d); do mv $i ../; done
for i in $(find . -name faa -type d); do mv $i ../; done
rm -fr *
mv ../gbk ../fna ../faa .

cd ../

echo " Finished Bacterial Plasmids fetching"


echo " Creating Listing Files.."
find Bacteria* -name "*.faa" > All-Proteins.listing.txt
find Bacteria* -name "*.gbk" > All-Genbanks.listing.txt
find Bacteria* -name "*.fna" > All-Contigs.listing.txt

echo ""
echo " Creating Output Directories :"
echo "     All-Bacterial-Contigs/"
echo "     All-Bacterial-Proteins/"
echo "     All-Bacterial-Genbanks/"

mkdir All-Bacterial-Contigs
mkdir All-Bacterial-Proteins
mkdir All-Bacterial-Genbanks

cd All-Bacterial-Contigs
for i in `find ../ -name "*.fna"`
do
    ln -s $i .
done
cd ../

cd All-Bacterial-Proteins
for i in `find ../ -name "*.faa"`
do
    ln -s $i .
done
cd ../

cd All-Bacterial-Genbanks
for i in `find ../ -name "*.gbk"`
do
    ln -s $i .
done
cd ../


## Create Blast Database ##
if [[ $blast -eq 1 ]]
then
    export PATH=/is1/commonPrograms/ncbi-blast-2.2.29+-src/c++/BUILD/bin:$PATH
    mkdir BlastDB
    makeblastdb -in uniprotKB_$date.fasta -input_type fasta -dbtype prot
    mv uniprotKB_$date.fasta.* BlastDB
fi
