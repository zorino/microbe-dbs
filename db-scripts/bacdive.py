#!/usr/bin/env python
# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maximilien1er@gmail.com
# date:    	2017-11-16
# version: 	0.01

import sys
import requests
import json
import time
from random import shuffle

bacdive_url = "https://bacdive.dsmz.de/api/bacdive/bacdive_id/"

def crawl_info(u, p, ids, outdir):

    print("# Crawling Bacterial MetaData..")
    header = {'Accept': 'application/json'}

    outfile = outdir+"/zz-ids.txt"
    with open(outfile, 'w') as f:
        f.write("\n".join(ids))

    for _id in ids:
        worked = False
        tried = 1
        while not worked:
            try:
                time.sleep(2*tried)
                resp = requests.get(bacdive_url+"/%s/"%_id, auth=(u, p), headers=header)
                if(resp.ok):
                    data = resp.json()
                    outfile = outdir+"/%s.json"%_id
                    with open(outfile, 'w') as f:
                        f.write(json.dumps(data, indent=2, separators=(',', ': ')))
                    worked = True
            except:
                print("ID: %s  failed.. retrying" % _id)
                tried += 5

def crawl_ids(u,p):
    ids = []
    data = {}
    bacdive_current_url = bacdive_url
    header = {'Accept': 'application/json'}
    iterator = 0

    print("# Crawling IDs..")
    while bacdive_current_url != "null":
        resp = requests.get(bacdive_current_url, auth=(u,p), headers=header)
        if(resp.ok):
            data = resp.json()
            for r in data['results']:
                ids.append(r['url'].replace(bacdive_url, "")[0:-1])
            bacdive_current_url = data['next']
            time.sleep(3)
        else:
            print('Response not ok')
            bacdive_current_url = "null"

        if bacdive_current_url is None:
            bacdive_current_url = "null"
        elif bacdive_current_url[0:4] != "http":
            bacdive_current_url = "null"

    print(" ")

    return ids

# Main #
if __name__ == "__main__":
    user = sys.argv[1]
    password = sys.argv[2]
    output_dir = sys.argv[3]
    ids = crawl_ids(user, password)
    print(" Number of IDs: %i" % len(ids))
    crawl_info(user, password, ids, output_dir)
