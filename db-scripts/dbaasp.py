#!/usr/bin/env python
# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maximilien1er@gmail.com
# date:    	2022-07-13
# version: 	0.01

import sys
import os
import libxml2
import requests
from lxml import etree
from lxml.html import fromstring, tostring
import json
import urllib3

urllib3.disable_warnings()

URL = "https://dbaasp.org/peptides/"


def get_number_of_proteins():
    r = requests.get(URL)
    data = r.json()
    print("# %s proteins to download .. " % data['totalCount'])
    return int(data['totalCount'])


def get_peptides(nb_prot, outdir):

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    iterator = int(nb_prot / 1000)
    for i in range(0, iterator):
        _url = URL + ("?offset=%d&limit=1000" % (i * 1000))
        print(_url)
        r = requests.get(_url)
        data = r.json()
        ids = [d['id'] for d in data['data'] if 'id' in d]
        for id in ids:
            _prot_url = URL + str(id)
            _r = requests.get(_prot_url)
            _prot = _r.json()
            with open(outdir + '/' + str(id) + '.json', 'w') as outfile:
                json.dump(_prot, outfile)

        print(ids)
    return 0


# Main #
if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("dbaasp.py <output dir>")
        exit(1)

    outdir = sys.argv[1]

    nb_prot = get_number_of_proteins()
    get_peptides(nb_prot, outdir)
