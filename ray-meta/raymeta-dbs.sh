#!/bin/bash
# author : maxime déraspe
# email : maximilien1er@gmail.com
#
# Script to update NCBI Genbank Bacterial Genomes
# fetch those files :
# ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/

date=$(date +%Y-%m-%d)

if [[ $1 == "" ]];then
    echo "Need at least one database to download (-g, -p, -f) !"
    echo " BacterialDB-.sh"
    echo "     option -u <directory> will update the databases you specify [-g,-p,-f] in their subdirectory"
    echo "     option -g download ncbi genbank bacteria (ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/)"
    echo "     option -p download ncbi plasmids (ftp://ftp.ncbi.nlm.nih.gov/genomes/Plasmids/)"
    echo "     option -f download ebi phages (ftp://ftp.ebi.ac.uk/pub/databases/fastafiles/embl_genomes/genomes/Phage/)"
    # echo "     option -b create blast database (BLAST executables need to be in your path)"
    exit 0
else
    while getopts "g f p b h u:" opt
    do
	      case $opt in
	          u)
		            update=$OPTARG
		            ;;
	          g)
		            genome=1
		            ;;
	          f)
		            phages=1
		            ;;
	          p)
		            plasmid=1
		            ;;
            b)
		            blast=1
		            ;;
            h)
		            echo "Need a Release Dates as argument (-d) and a DB to download (-g, -p, -f) !"
		            echo " BacterialDB-loader.sh -d <date>"
		            echo "     option -u will update the database you specify"
		            echo "     option -g download ncbi genbank bacteria (ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/)"
		            echo "     option -p download ncbi plasmids (ftp://ftp.ncbi.nlm.nih.gov/genomes/Plasmids/)"
		            echo "     option -f download ebi phages (ftp://ftp.ebi.ac.uk/pub/databases/fastafiles/embl_genomes/genomes/Phage/)"
		            # echo "     option -b create blast database (BLAST executables need to be in your path)"
		            exit 0
		            ;;
            \?)
		            echo "Invalid Options !"
		            exit 0
		            ;;
	      esac
    done
fi


if [[ $x$update != $x ]]
then
    if [[  ! -d "$update" ]]
    then
	      echo "The Directory to update $update doesn't exist !"
    else
	      cd $update
    fi
else
    mkdir BacterialDB_$date
    cd BacterialDB_$date
fi


## Fetch NCBI Bacterial Genomes ##
if [[ $genome -eq 1 ]]
then
    echo "#[NCBI Genbank Bacterial Genomes]#"
    echo "Listing NCBI Genbank Bacterial Genomes from ftp"
    wget -q --no-parent --no-remove-listing --spider ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/
    mv .listing Bacteria-NCBI-$date.listing
    nb_of_genomes=`wc -l < Bacteria-NCBI-$date.listing`
    iterator=1

    ## Go into Bacteria-NCBI root directory
    if [[ $x$update != $x ]]
    then
	      # Go into last Bacteria-NCBI directory
	      dir=$(find . -maxdepth 1 -type d -name  "Bacteria-NCBI-*")
	      if [[ $dir == "" ]]
	      then
	          echo "Can't Update Bacteria-NCBI.. are you sure you already downloaded it ?"
	          exit 0
	      else
	          cd $dir
	      fi
    else
	      mkdir Bacteria-NCBI-$date
	      cd Bacteria-NCBI-$date
    fi


    # Reading each species of root NCBI Bacterial from ftp structure
    while read -r i
    do
	      # dr-xr-xr-x   5 ftp      anonymous     4096 Mar 10 04:50 Abiotrophia_defectiva
	      if [ ! -d $i ]
	      then
	          mkdir $i
	          echo "NEW SPECIES: $i" >> ../NEWS_$date.logs
	      fi

	      cd $i

	      echo "Downloading/Checking $i assemblies [$iterator/$nb_of_genomes]"
	      iterator=$(($iterator+1))

	      break_iterator=5
	      while [ $break_iterator -gt 0 ]
	      do
	          wget -q ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/$i/assembly_summary.txt -O assembly_summary_$date.txt
	          if [ -f assembly_summary_$date.txt ]
	          then
		            break_iterator=0
	          else
		            break_iterator=$(($break_iterator-1))
		            sleep 3
	          fi
	      done

	      # Look if assembly_summary file has been downloaded
	      if [ ! -f assembly_summary_$date.txt ]
	      then
	          echo "ERROR $i: downloading assembly_summary_$date.txt" >> ../../DownloadErrors_$date.logs
	          cd ../
	          continue
	      fi

	      # Reading assembly_summary and downloading each assembly in it
	      while IFS=$";" read -r -a j
	      do
	          # assembly_accession header
	          if [[ "${j[0]}" == "# assembly_accession" ]]
	          then
		            continue
	          fi

	          assembly_version=$(echo ${j[15]} | sed "s/ /_/g" | sed "s#/#_#g")
	          assembly="${j[0]}_$assembly_version"

	          echo $assembly

	          if [ ! -d $assembly ]
	          then
		            echo "NEW ASSEMBLY: $i $assembly" >> ../../NEWS_$date.logs
		            mkdir $assembly
		            wget -q -r -nH --cut-dirs=2 --retr-symlinks \
		                 --no-remove-listing --retr-symlinks \
		                 --no-remove-listing --passive-ftp ${j[19]}
	          fi

	          cd $assembly

	          md5sum=$(md5sum -c md5checksums.txt 2> /dev/null | grep "$assembly_genomic" | grep -v OK$)

	          if [[ $md5sum != "" ]]
	          then
		            echo "ERROR $i: md5sum $md5sum" >> ../../../DownloadErrors_$date.logs
	          fi

	          cd ../

	      done < <(sed "s/\t/;/g" assembly_summary_$date.txt)

	      cd ../

    done < <(awk '{print $9}' ../Bacteria-NCBI-$date.listing | sed "s/\\r//g")

fi


## NCBI Bacterial Plasmids ##
if [[ $plasmid -eq 1 ]]
then
    mkdir Plasmids-NCBI-$date
    cd Plasmids-NCBI-$date
    echo "Downloading Bacterial Plasmids"

    wget -q ftp://ftp.ncbi.nlm.nih.gov/genomes/Plasmids/plasmids.all.gff.tar.gz
    tar xvf plasmids.all.gff.tar.gz
    mv am/ftp-genomes/Plasmids/gff/ ./
    rm -fr am

    wget -q ftp://ftp.ncbi.nlm.nih.gov/genomes/Plasmids/plasmids.all.faa.tar.gz
    tar xvf plasmids.all.faa.tar.gz
    mv am/ftp-genomes/Plasmids/faa/ ./
    rm -fr am

    wget -q ftp://ftp.ncbi.nlm.nih.gov/genomes/Plasmids/plasmids.all.fna.tar.gz
    tar xvf plasmids.all.fna.tar.gz
    mv am/ftp-genomes/Plasmids/fna/ ./
    rm -fr am

    cd ../
fi

## EBI Phage Genomes ##
if [[ $phages -eq 1 ]]
then
    mkdir Phages-EBI-$date
    cd Phages-EBI-$date
    echo "Downloading Bacterial Phages"
    wget -q -r -l1 --no-parent -nd \
	       --passive-ftp ftp://ftp.ebi.ac.uk/pub/databases/fastafiles/embl_genomes/genomes/Phage/
    echo "Extracting Phages and Renaming fasta files"
    for i in *.gz
    do
        gunzip $i
        f=${i%.fasta.gz}
        header=`head -1 $f.fasta`
        name=`echo $header | \
                sed "{
                s/>.* STD://g
                s/ complete .*//g
                s/ /_/g
                s/,//g
                s/\//-/g
                }"`
        c=0
        while read line
        do
            if [ ${line:0:1} == ">"  ]
            then
                c=$(($c+1))
            fi
            echo $line >> $name-$c.fasta
        done < $f.fasta
        rm $f.fasta
    done
    cd ../
fi


# ## Create Blast Database ##
# if [[ $blast -eq 1 ]]
# then
#     # export PATH=/is1/commonPrograms/ncbi-blast-2.2.29+-src/c++/BUILD/bin:$PATH
#     # Creating BlastDB
#     mkdir BlastDB
#     cd BlastDB
#     # find ../ -name *.fna.gz | parallel -j $core "zcat {} >> all-contigs.fna"
#     # makeblastdb -in uniprotKB_$date.fasta -input_type fasta -dbtype prot
#     # mv uniprotKB_$date.fasta.* BlastDB
# fi
