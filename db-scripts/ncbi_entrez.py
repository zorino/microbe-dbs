#!/usr/bin/env python
# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maximilien1er@gmail.com
# date:    	2017-05-25
# version: 	0.01

import sys
from time import sleep
from Bio import Entrez


# Main #
if __name__ == "__main__":

    Entrez.email = "microbesdbs@example.com"

    project_id = sys.argv[1]
    output = sys.argv[2]

    handle = Entrez.elink(dbfrom="bioproject", id=project_id, linkname="bioproject_nuccore")
    records = Entrez.read(handle)
    # print(records[0]['LinkSetDb'][0]['Link'])
    for r in records[0]['LinkSetDb'][0]['Link']:
        sleep(1)
        print("".join(Entrez.efetch(db='nuccore', id=r['Id'], rettype=output).readlines()))
    handle.close()
