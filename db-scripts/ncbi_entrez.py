#!/usr/bin/env python
# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maximilien1er@gmail.com
# date:    	2017-05-25
# version: 	0.01

import sys
from time import sleep
from Bio import Entrez


def bioproject_nuccore(project_id, output):

    Entrez.email = "microbesdbs@example.com"

    handle = Entrez.elink(dbfrom="bioproject", id=project_id, linkname="bioproject_nuccore")
    records = Entrez.read(handle)
    # print(records[0]['LinkSetDb'][0]['Link'])
    for r in records[0]['LinkSetDb'][0]['Link']:
        sleep(1)
        print("".join(Entrez.efetch(db='nuccore', id=r['Id'], rettype=output).readlines()))
    handle.close()



# Main #
if __name__ == "__main__":


    usage = """
ncbi_entrez.py program [options]

program:

  bioproject_nuccore  <project_id> <output>

"""


    if len(sys.argv) < 2:
        print(usage)
        sys.exit(1)


    if sys.argv[1] == "bioproject_nucccore":
        if len(sys.argv) < 4:
            print(usage)
            sys.exit(1)

        project_id = sys.argv[2]
        output = sys.argv[3]

        bioproject_nuccore(project_id, output)
