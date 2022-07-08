#!/usr/bin/env python
# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maximilien1er@gmail.com
# date:    	2022-04-10
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

URL = "https://aps.unmc.edu"


def get_protein_list():

    print("# Downloading ADP3 AMP list..")
    peptites_by_family = {}
    r = requests.get(URL + "/home", verify=False)
    page = fromstring(r.text)
    for db in page.xpath('//table//tr//td//form'):
        params = db.xpath('input/@value')
        activity = params[1].replace(" ", "+")
        _r = requests.post(URL + "/database/anti",
                           data={
                               'activity': params[0],
                               'name': activity
                           },
                           verify=False)
        _page = fromstring(_r.text)
        peptites_by_family[activity] = []
        for prot in _page.xpath('//form'):
            peptites_by_family[activity].append({
                'url': URL + "/database/peptide",
                'data': {
                    'ID': prot.xpath('input/@value')[0]
                }
            })

    return peptites_by_family


def get_peptide_family(family, outdir):

    print("# Downloading ADP3 families..")
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    for fam in peptides_by_family.keys():
        print("   - %s" % fam)
        _outfamdir = outdir + "/" + fam
        if not os.path.exists(_outfamdir):
            os.mkdir(_outfamdir)

        prot_info = []

        for prot in peptides_by_family[fam]:

            _r = requests.post(prot['url'], data=prot['data'], verify=False)
            _p = fromstring(_r.text)
            prot = {}
            for info in _p.xpath("//tr"):
                key = str(
                    info.xpath('td')[0].text_content()).strip().strip(':')
                value = str(info.xpath('td')[1].text_content()).strip()
                prot[key] = value
                if key == "SwissProt ID":
                    value = value.replace("SwissProt ID: ",
                                          "").replace(" Go to SwissProt", "")
                elif key == "Additional Info":
                    value = info.xpath('td')[1].text_content()
                elif key == "Reference":
                    value = info.xpath('td')[1].xpath('a')

            with open(_outfamdir + '/' + prot['APD ID'] + '.json',
                      'w') as outfile:
                json.dump(prot, outfile)


# Main #
if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("adp3.py <output dir>")
        exit(1)

    outdir = sys.argv[1]

    peptides_by_family = get_protein_list()
    get_peptide_family(peptides_by_family, outdir)
