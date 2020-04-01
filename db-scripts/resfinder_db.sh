#!/bin/bash
# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maximilien1er@gmail.com


# download files from web site
function download_files() {

    if [ -z $(command -v git) ]
    then
        echo "Sorry install git for this one"
        exit 1
    fi

    git clone https://bitbucket.org/genomicepidemiology/resfinder_db.git _resfinder_db_tmp
    release=$(git --git-dir=./_resfinder_db_tmp/.git log -1 --format=%at | xargs -I{} date -d @{} +%Y-%m-%d)

    output=""
    if [ -z $1 ]
    then
        # mv _resfinder_db_tmp resfinder_db_$release
        output="resfinder_db_$release"
    else
        if [ ! -d $1 ]
        then
            mkdir $1
        fi
        # mv _resfinder_db_tmp $1/resfinder_db_$release
        output="$1/resfinder_db_$release"
    fi
    mv _resfinder_db_tmp $output
    cd $output
}

# organize files into sub-directories
function organize_files () {
    # dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
    # python $dir/isfinder.py IS.faa IS.csv > IS.tsv
    ls
}

download_files $1
organize_files
