#!/bin/bash
# author : maxime déraspe
# email : maximilien1er@gmail.com
#

date=$(date +%Y-%m-%d)

DIR=$(dirname "$(readlink -f $0)")

cli_help=$(cat << EOF

     microbe-dbs.sh <database> <output>
        databases :

          [Bacteria]
            refseq_bacteria           RefSeq bacteria from NCBI ftp
            refseq_plasmid            RefSeq plasmid from NCBI ftp
            cog                       Clusters of Orthologous Groups (COG)
            pub_mlst                  Public databases for molecular typing and microbial genome diversity
            vfdb                      Virulence Factor Database
            mibig                     Minimum Information about a Biosynthetic Gene cluster
            ncbi_arg                  NCBI Bacterial Antimicrobial Resistance Reference Gene Database
            bacdive                   DSMZ Bacterial Diversity Metadatabase [Need User/Password]

          [Virus]
            refseq_viral              RefSeq viral from NCBI ftp

          [Other]
            ncbi_taxonomy             NCBI Taxonomy from ftp
            refseq_archaea            RefSeq archaea from NCBI ftp
            refseq_fungi              RefSeq fungi from NCBI ftp
            ebi_idmapping             ID mappings for proteins (Uniprot -> RefSeq, PDB, BioCyc)

     Example :
        $ ./microbe-dbs refseq_bacteria refseq_bacteria_output

EOF
        )

if [[ $1 == "" ]]; then
    echo "Need at least 1 databases !!"
    echo "$cli_help"
    exit 0
else
    if [ -z $2 ]; then
        output="$1"
    else
        output=$2
    fi
    case "$1" in
        refseq_bacteria)
            $DIR/db-scripts/refseq.sh $output "bacteria"
            ;;
        refseq_plasmid)
            $DIR/db-scripts/refseq.sh $output "plasmid"
            ;;
        refseq_viral)
            $DIR/db-scripts/refseq.sh $output "viral"
            ;;
        refseq_archaea)
            $DIR/db-scripts/refseq.sh $output "archaea"
            ;;
        refseq_fungi)
            $DIR/db-scripts/refseq.sh $output "fungi"
            ;;
        cog)
            $DIR/db-scripts/cog.sh $output
            ;;
        vfdb)
            $DIR/db-scripts/vfdb.sh $output
            ;;
        mibig)
            $DIR/db-scripts/mibig.sh $output
            ;;
		bacdive)
            $DIR/db-scripts/bacdive.sh $output
            ;;
        ncbi_arg)
            $DIR/db-scripts/ncbi_arg.sh $output
            ;;
        pub_mlst)
            $DIR/db-scripts/pub_mlst.sh $output
            ;;
        ebi_idmapping)
            $DIR/db-scripts/ebi_idmapping.sh $output
            ;;
        ncbi_taxonomy)
            $DIR/db-scripts/ncbi_taxonomy.sh $output
            ;;
        help)
            echo "$cli_help"
		        ;;
        \?)
		        echo "Invalid Database !"
		        exit 0
		        ;;
    esac
fi
