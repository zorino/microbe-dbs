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

    if [ -z $1 ]
    then
        mv _resfinder_db_tmp resfinder_db_$release
    else
        if [ ! -d $1 ]
        then
            mkdir $1
        fi
        mv _resfinder_db_tmp $1/resfinder_db_$release
    fi

}

download_files $1
