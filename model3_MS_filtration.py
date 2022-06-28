#!/usr/bin/env python
# -*- coding: utf-8 -*-i


"""
Mass Spectrometry Data Filtering:
"""

import os
import sys

def handle01_create_db():
    cmd1='cat ./mut_result/mut_pro_pl.csv.fasta ./mut_result/mut_pro_nl.csv.fasta ./reference/ref_uniprot_crap.fasta > ./mut_result/mut_ref_db.fasta'
    os.system(cmd1)


def handle02_xml():
    cmd1='python ./biosoft/gen_mqpar.py ./biosoft/labelfree.xml /path/to/PGNneo/ms -o ./biosoft/mqpar.xml -t 6' 
    os.system(cmd1)


def handle03_maxquant():
    cmd1= "mono ./biosoft/MaxQuant/bin/MaxQuantCmd.exe ./biosoft/mqpar.xml"
    os.system(cmd1)


if __name__ == '__main__':
    handle01_create_db()
    handle02_xml()
    handle03_maxquant()
    print("Model3 Run Successful!")
    #python model3_MS_filtration.py
