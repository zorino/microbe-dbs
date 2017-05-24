#!/usr/bin/env python
# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maximilien1er@gmail.com
# date:    	2017-05-24
# version: 	0.01

import sys
import xml.etree.ElementTree as ET

# Main #
if __name__ == "__main__":

    tree = ET.parse('dbases.xml')
    root = tree.getroot()
    for species in root.findall('species'):
        species_name = species.text.strip()
        count = species.findall(".//profiles/count")[0].text
        url = species.findall(".//profiles/url")[0].text
        for locus in species.findall(".//loci/locus"):
            locus_name = locus.text.strip()
            locus_url = locus.findall(".//url")[0].text
            #print(locus.findall(".//url").text.strip())
            print("%s\t%s\t%s\t%s\t%s" % (species_name, count, url, locus_name, locus_url))
