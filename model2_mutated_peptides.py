#!/usr/bin/env python
# -*- coding: utf-8 -*-i


"""
Generation of mutated peptides:
"""

from os import remove
from functools import reduce
import os
import pandas as pd
import sys
import subprocess


def rev_comp(s):
    sc = ''
    for i in reversed(s):
        if i == 'A':
            sc += 'T'
        elif i == 'C':
            sc += 'G'
        elif i == 'G':
            sc += 'C'
        elif i == 'T':
            sc += 'A'
    return sc


def get_mutation(f, p_df):
    train_pos_data_path = f 
    alt_value = p_df 
    fh = open(train_pos_data_path, 'r')
    seq = {}
    for line in fh:
        if line.startswith('>'):
            name = line.replace('>', '').split()[0]
            seq[name] = ''
        else:
            seq[name] += line.replace('\n', '')
    fh.close()
    wild_list = []
    for key in seq:
        wild_list.append(seq[key])
    mut_list = list()
    mut_pos_list = alt_value['ALT'].tolist()
    for i in range(0, len(wild_list)):
        re_front = wild_list[i][0:100]
        re_back = wild_list[i][-100:]
        re = re_front + mut_pos_list[i] + re_back
        re = re.upper()
        mut_list.append(re)
    out_dict = {"ID": alt_value["ID"].tolist(), "POS": alt_value["POS"], "X.CHROM": alt_value["X.CHROM"].tolist(),
                "ALT": alt_value["ALT"].tolist(), "REF": alt_value["REF"].tolist(), "Gene": alt_value["Gene"].tolist(),
                "wild": wild_list, "mut": mut_list
                }
    out_df = pd.DataFrame(out_dict)
    return out_df


def sixFrame_translate_fun(strs):
    lst1 = strs.strip("\n")
    lst1 = lst1.upper()
    lst1 = lst1.replace('N', 'A')
    trans_cod1 = "".join([codon_table["".join(lst1[i:i + 3])] for i in range(0, len(lst1) - len(lst1) % 3, 3)])
    trans_cod2 = "".join(
        [codon_table["".join(lst1[i:i + 3])] for i in range(1, (len(lst1) - 1) - (len(lst1) - 1) % 3, 3)])
    trans_cod3 = "".join(
        [codon_table["".join(lst1[i:i + 3])] for i in range(2, (len(lst1) - 2) - (len(lst1) - 2) % 3, 3)])
    return ([trans_cod1, trans_cod2, trans_cod3])


def align_translate_func(srs, order_trans="pl"):
    if order_trans == "pl":
        mut_base_str = srs["mut"]
        wild_base_str = srs["wild"]
    elif order_trans == "nl":
        mut_base_str = rev_comp(srs["mut"])
        wild_base_str = rev_comp(srs["wild"])
    else:
        exit("Either pl or nl is permitted in align_translate_func")

    tr1_list = sixFrame_translate_fun(mut_base_str)
    tr2_list = sixFrame_translate_fun(wild_base_str)
    tmp_list = list(map(lambda x, y: x if (x != y) else "", tr1_list, tr2_list))
    tmp_seq = ["trans" + str(i) for i in range(1, len(tr1_list) + 1)]
    out_dict = {'ID': srs["ID"] + "|" + order_trans, "POS": srs["POS"], "X.CHROM": srs["X.CHROM"],
                "ALT": srs["ALT"], "REF": srs["REF"], "Gene": srs['Gene'],
                'mut_pro': tmp_list, 'mut_base': mut_base_str,
                'wild_pro': tr2_list, 'wild_base': wild_base_str,
                "trans": tmp_seq}
    out_df = pd.DataFrame(out_dict)
    return (out_df)


def mk_fastaStr(pdSeries):
    fastaStr_lst = []
    for strs in ["prot3", "prot2", "prot1"]:
        ant_prot = strs + "-" + 'ant'
        mut_prot = strs + "-" + 'mut'
        mut_prot_pos = mut_prot + "-" + "pos"
        if (pd.isna(pdSeries[ant_prot]) or len(pdSeries[ant_prot]) < 7):
            continue
        else:
            try:
                header = pdSeries["ID"] + "|" + pdSeries["gene"] + "|" + pdSeries[mut_prot] + "|" + pdSeries[
                    mut_prot_pos] + "\n"
            except:
                print("Mutated protein %s translaste error", pdSeries["ID"])
            header = header.replace(";", "|")
            fastaStr = header + pdSeries[ant_prot]
        fastaStr_lst.append(fastaStr)
    return (fastaStr_lst)


if __name__ == '__main__':
    con=sys.argv[1]
    case = sys.argv[2]
    con = con.split('_')[0]
    case = case.split('_')[0]
    con_case=con+'_'+case
    initial_path=os.getcwd()
    os.chdir(initial_path+ '/mut_result/')
    index = initial_path+"/reference/hg38/hg38.fa"
    codon_table = {"TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
                   "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
                   "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
                   "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
                   "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
                   "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
                   "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
                   "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
                   "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
                   "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
                   "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
                   "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
                   "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
                   "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
                   "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
                   "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"}

    filename = initial_path+"/mut_result/"+con_case+"_anno_out.hg38_multianno.txt" 
    seq_df = pd.read_table(filename, sep="\t", header=0)
    non_exonic_df = seq_df[~seq_df['Func.refGene'].isin(('exonic', 'exonic;splicing'))]
    non_exonic_df = non_exonic_df[['Otherinfo4', 'Otherinfo5', 'Otherinfo7', 'Otherinfo8', 'Gene.refGene']]
    non_exonic_df.columns = ('X.CHROM', 'POS', 'REF', 'ALT', 'Gene')
    non_exonic_df['begin'] = non_exonic_df['POS'] - 101
    non_exonic_df['tail'] = non_exonic_df['POS'] + non_exonic_df['REF'].astype('str').str.len() - 1 + 100
    tmp_cat_series = '>' + non_exonic_df["X.CHROM"].str.cat(
        [non_exonic_df['POS'].astype("str"), non_exonic_df['REF'].astype('str'), non_exonic_df['ALT'].astype("str")],
        sep="|")
    non_exonic_df.insert(0, 'ID', tmp_cat_series)
    mut_pos_bed = non_exonic_df[['X.CHROM', 'begin', 'tail']]
    mut_pos_bed.to_csv(path_or_buf="./pos.bed", sep="\t", header=False, index=False)
    run_status = subprocess.run(
        ['bedtools', 'getfasta', '-fi', index, '-bed', './pos.bed', '-fo', 'wild_base_tmp.fasta'])
    non_exonic_df.to_pickle("non_exonic_df")
    if os.path.exists("./wild_base_tmp.fasta"):
        mut_seq_df = get_mutation("./wild_base_tmp.fasta", non_exonic_df)
        prot_pl_seq_srs = mut_seq_df.apply(align_translate_func, axis=1, order_trans="pl")
        mut_wild_prot_pl_df = pd.concat(prot_pl_seq_srs.tolist())
        prot_nl_seq_srs = mut_seq_df.apply(align_translate_func, axis=1, order_trans="nl")
        mut_wild_prot_nl_df = pd.concat(prot_nl_seq_srs.tolist())
   
        mut_prot_pl_pivot = mut_wild_prot_pl_df.pivot(index="ID", columns="trans", values="mut_pro")
        mut_prot_pre_tmp = mut_wild_prot_pl_df[
            ['ID', 'X.CHROM', 'POS', 'ALT', 'REF', 'Gene', 'mut_base']].drop_duplicates()
        mut_prot_pl_df = pd.merge(mut_prot_pre_tmp, mut_prot_pl_pivot, on="ID", how="inner")
        mut_prot_pl_df.to_csv(path_or_buf="mut_pro_pl.csv", index=False)
        mut_prot_nl_pivot = mut_wild_prot_nl_df.pivot(index="ID", columns="trans", values="mut_pro")
        mut_prot_pre_tmp = mut_wild_prot_nl_df[
            ['ID', 'X.CHROM', 'POS', 'ALT', 'REF', 'Gene', 'mut_base']].drop_duplicates()
        mut_prot_nl_df = pd.merge(mut_prot_pre_tmp, mut_prot_nl_pivot, on="ID", how="inner")
        mut_prot_nl_df.to_csv(path_or_buf="mut_pro_nl.csv", index=False)

    i = 1
    for f in ["mut_pro_pl.csv", "mut_pro_nl.csv"]:
        outFile = 'tst' + str(i) + '.csv'
        print("Mutation position finding ...")
        run_status = subprocess.run(['Rscript', '../mutpeppos.R', f, outFile])
        print("done")
        tmp = pd.read_csv(outFile, index_col=0)
        tmp.insert(0, "ID", tmp.index)
        tmp.to_pickle("tst_df")
        strSeries = tmp.apply(mk_fastaStr, 1)
        strSeries2lst = list(strSeries)
        strSeries2lst = [item for a in strSeries for item in a]
        mergedStr = reduce(lambda x, y: x + "\n" + y, strSeries2lst)
        fo = open(f + ".fasta", "w")
        fo.write(mergedStr)
        i = i + 1

    ##Delete intermediate files
    os.remove(initial_path+"/mut_result/con_case_1.vcf.idx")
    os.remove(initial_path+"/mut_result/con_case_1.vcf.stats")
    os.remove(initial_path+"/mut_result/con_case.vcf.idx")
    os.remove(initial_path+"/mut_result/con_case.vcf.filteringStats.tsv")
    os.remove(initial_path+"/mut_result/con_case_anno_out.avinput")
    os.remove(initial_path+"/mut_result/tst_df")    
    os.remove(initial_path+"/mut_result/mut_pro_pl.csv")
    os.remove(initial_path+"/mut_result/mut_pro_nl.csv")
    print("Model2 Run Successful!")
    ##python model2_mutated_peptides.py con_R1.fastq.gz case_R1.fastq.gz
