#!/bin/bash
# author : maxime déraspe
# email : maximilien1er@gmail.com
#
# required :
# - bash
# - python and dbaasp.py dependencies

# download adp
function download_dbaasp() {

    out_dir=""
    if [ -z $1 ]
    then
        outdir=dbaasp
    else
        outdir=$1/dbaasp
    fi
    mkdir -p $outdir

    dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
    python $dir/dbaasp.py $outdir

}

download_dbaasp $1
