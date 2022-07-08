#!/bin/bash
# author : maxime déraspe
# email : maximilien1er@gmail.com
#
# required :
# - bash
# - python and adp3.py dependencies

# download adp
function download_adp() {

    out_dir=""
    if [ -z $1 ]
    then
        outdir=adp3
    else
        outdir=$1/adp3
    fi
    mkdir -p $outdir

    dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
    python $dir/adp3.py $outdir

}

download_adp $1
