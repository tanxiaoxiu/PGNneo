#!/usr/bin/env python
# -*- coding: utf-8 -*-i


"""
Neoantigen prediction and selection:
"""
import sys
import os
import pandas as pd
import re
import subprocess


def handle01_maxqpep(inputPath,outPath):
    inputFile = inputPath + "peptides.txt"
    outFile = outPath + "maxqpep.txt"
    readFile = open(inputFile, 'r')
    writeFile = open(outFile, 'w+')
    for line in readFile:
        ls1 = line.strip('\n').split('\t')
        if ls1[0] == 'Sequence':  
            writeFile.write(
                'Sequence' + '\t' + 'Gene' + '\t' + 'AA Change' + '\t' + 'Mut Position' + '\t' + 'Praw' + '\t' + 'Chain' + '\t' + 'Base' + '\t' + '\n')
        ls2 = re.split('[|,;]', ls1[34])  
        if ('uniprot' not in ls2):
            ls2 = ls1[34].split(';')
            for i in range(len(ls2)):
                ls3 = ls2[i].split('|')  
                if len(ls3) >= 8:  
                    ls4 = ls3[-2]
                    ls5 = ls3[-1].split('-')
                    ls6 = ls5[0]
                    if len(ls1[36]) != 0:
                        if float(ls1[36]) <= int(str(ls6)) and int(str(ls6)) <= float(ls1[37]):
                            ls7 = float(ls6) - float(ls1[36])
                            if float(ls7) < float(len(ls1[0])):
                                if ls1[0][int(ls7)] == ls4:
                                    writeFile.write(
                                        ls1[0] + '\t' + ls3[5] + '\t' + ls4 + '\t' + str(int(ls7) + 1) + '\t' + ls3[
                                            -2] + ls3[-1] + '\t' + ls3[4] + '\t' + ls3[0] + ls3[2] + ls3[1] + ls3[
                                            3] + '\n')



def handle02_cutpep(outPath):
    os.system("Rscript cutpep.R " + outPath)


def handle03_preneo(outPath):
    hla_allele1 = input("please input an HLA class I allele like 'HLA-A02:01' or multiple alleles like 'HLA-A02:01,HLA-B15:01,HLA-C01:02':")
    hla_allele1 = hla_allele1.replace(' ', '')
    for i in range(8,12):
        cmd1='netMHCpan -a '+hla_allele1+' -l '+str(i)+' '+outPath+'maxqpep_'+str(i)+'.fasta > '+outPath+'netM_'+str(i)+'.csv'
        os.system(cmd1)
    cmd2='cat '+outPath+'netM_8.csv'+' '+outPath+'netM_9.csv'+' ' +outPath+'netM_10.csv'+' ' +outPath+'netM_11.csv > '+outPath+'netM.csv'
    os.system(cmd2)  


def handle04_bindingneo(outPath):
    filename = outPath + "netM.csv"
    outFile = outPath + "neo1.csv"
    readFile = open(filename, 'r')
    writeFile = open(outFile, 'w+')

    for line in readFile:
        ls1 = line.strip(" ").split('\00')
        ls2 = ls1[0].strip(' ').split(' ')
        if "WB" in ls2[-1] or 'SB' in ls2[-1]:
            ls3 = '\t'.join(ls2[:])
            writeFile.write(ls3)

    readFile.close()
    writeFile.close()


def handle05_bindingneo(outPath):
    infile = outPath + "neo1.csv"
    outfile = outPath + "neo2.csv"

    data = pd.read_table(infile, sep="\s+", header=None)
    temp = data.iloc[:, [1, 2, -5, -3, -1]]
    temp.to_csv(outfile, sep="	", header=['HLA', 'Peptide', 'Gene', '%Rank', 'BindLevel'], index=None)


def handle06_bindingneo(outPath):
    filename = outPath + "neo2.csv"
    outname = outPath + "neo.csv"
    readFile = open(filename, 'r')
    writeFile = open(outname, 'w+')

    for line in readFile:
        ls = line.strip().split("\t")
        ls1 = list(ls[1])
        if 'X' not in ls1:
            ls2 = ''.join(ls1)
            writeFile.write(ls[0] + '\t' + ls2 + '\t' + ls[2] + '\t' + ls[3] + '\t' + ls[4] + '\n')


def handle07_getfasta(outPath):
    filename = outPath + "neo.csv"
    outFile = outPath + "neofa.fasta"
    readFile = open(filename, 'r')
    writeFile = open(outFile, 'w+')

    lines = readFile.readlines()[1:]
    for line in lines:
        ls = line.strip('\n').split('\t')
        writeFile.write(">" + ls[0] + "|" + ls[2] + '\n' + ls[1] + '\n')
    readFile.close()
    writeFile.close()

def handle08_filterneo(outPath,refdbPath):
	cmd1 = 'makeblastdb -in '+ refdbPath + "hcneo_db.fasta" +' -dbtype prot -title hcneo_db -parse_seqids -out '+refdbPath + 'hcneo_db'
	cmd2 = 'blastp -query '+outPath+ "neofa.fasta"+ ' -db '+ refdbPath + 'hcneo_db'+' -outfmt "6 qacc qseq sacc sseq evalue length pident" -evalue 100000000 -gapopen 11 -gapextend 1  > '+outPath + 'filterneo.txt'
	os.system(cmd1)
	os.system(cmd2) 
    

if __name__ == '__main__':
    inputPath = "./ms_resultmqpar/combined/txt/"
    outPath = "./preneo/"
    refdbPath = "./reference/hcneodb/"
    handle01_maxqpep(inputPath,outPath)
    handle02_cutpep(outPath)
    handle03_preneo(outPath)
    handle04_bindingneo(outPath)
    handle05_bindingneo(outPath)
    handle06_bindingneo(outPath)
    handle07_getfasta(outPath)
    handle08_filterneo(outPath,refdbPath)
    ##Delete intermediate files
    os.remove(outPath+"netM_8.csv")
    os.remove(outPath+"netM_9.csv")
    os.remove(outPath+"netM_10.csv")
    os.remove(outPath+"netM_11.csv")
    os.remove(outPath+"neo1.csv")
    os.remove(outPath+"neo2.csv")
    print("Model4 Run Successful!")
    #python model4_neoantigen_prediction_filtration.py
