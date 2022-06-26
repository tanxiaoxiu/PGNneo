#!/usr/bin/env python
# -*- coding: utf-8 -*-i


"""
Noncoding somatic variant calling and HLA typing:
"""

import os
import sys

def rnaseq_processing(input1,input2):
    con_case=input1.split('_')[0]
    cmd1='java -jar ./biosoft/trimmomatic-0.39-1/trimmomatic.jar PE -phred33 ./rnaseq/'+input1+' ./rnaseq/'+input2\
         +' ./rna_result/'+con_case+'_cut_R1.fastq.gz ./rna_result/'+con_case+'_cut_unpaired_R1.fastq.gz ./rna_result/'+con_case+'_cut_R2.fastq.gz ./rna_result/'+con_case+'_cut_unpaired_R2.fastq.gz ILLUMINACLIP:trimmomatic-0.39-1/adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:5:20 LEADING:5 TRAILING:5 HEADCROP:10 MINLEN:50'
    cmd2="bwa mem -t 16 -R '@RG\\tID:foo\\tSM:bar\\tLB:library1' ./reference/hg38/hg38.fa ./rna_result/"+con_case+"_cut_R1.fastq.gz ./rna_result/"+con_case+"_cut_R2.fastq.gz > ./rna_result/"+con_case+"_cut.sam"
    cmd3='samtools fixmate -O bam ./rna_result/'+con_case+'_cut.sam ./rna_result/'+con_case+'_cut_fixmate.bam'

    cmd4 = 'samtools sort -O bam -o ./rna_result/'+con_case+'_cut_sorted.bam -T ./rna_result/'+con_case+'_cut_temp ./rna_result/'+con_case+'_cut_fixmate.bam'

    cmd5 = 'gatk MarkDuplicates -I ./rna_result/'+con_case+'_cut_sorted.bam -O ./rna_result/'+con_case+'_cut_marked_duplicates.bam -M ./rna_result/'+con_case+'_cut_marked_dup_metrics.txt'

    cmd6 = 'gatk AddOrReplaceReadGroups -I ./rna_result/'+con_case+'_cut_marked_duplicates.bam -O ./rna_result/'+con_case+'_cut_marked_duplicates_1.bam -ID 4 -LB lib1 -PL illumina -PU unit1 -SM 20'
    cmd7 = 'gatk BaseRecalibrator -I ./rna_result/'+con_case+'_cut_marked_duplicates_1.bam -R ./reference/hg38/hg38.fa --known-sites ./reference/dbsnp_146.hg38.vcf.gz --known-sites ./reference/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz -O ./rna_result/'+con_case+'_cut_recal_data.table'
    cmd8 = 'gatk ApplyBQSR -R ./reference/hg38/hg38.fa -I ./rna_result/'+con_case+'_cut_marked_duplicates_1.bam --bqsr-recal-file ./rna_result/'+con_case+'_cut_recal_data.table -O ./rna_result/'+con_case+'_cut_recal.bam'

    cmd9 = 'picard AddOrReplaceReadGroups I=./rna_result/'+con_case+'_cut_recal.bam O=./rna_result/'+con_case+'_recal.bam RGID='+con_case+' RGLB=library1 RGPL=illumina RGPU=unit1 SORT_ORDER=coordinate RGSM='+con_case 
    cmd10 = 'samtools index ./rna_result/'+con_case+'_recal.bam' 
    cmd11='rm ./rna_result/'+con_case+'_cut_unpaired_R1.fastq.gz'
    cmd12 = 'rm ./rna_result/'+con_case+'_cut_unpaired_R2.fastq.gz'
    cmd13= 'rm ./rna_result/'+con_case+'_cut.sam'
    cmd14='rm ./rna_result/'+con_case+'_cut_fixmate.bam'
    cmd15 = 'rm ./rna_result/'+con_case+'_cut_sorted.bam'
    cmd16 = 'rm ./rna_result/'+con_case+'_cut_marked_duplicates.bam'
    cmd17 = 'rm ./rna_result/'+con_case+'_cut_marked_duplicates_1.bam'
    cmd18 = 'rm ./rna_result/'+con_case+'_cut_recal.bam'
    command=[cmd1,cmd2,cmd3,cmd4,cmd5,cmd6,cmd7,cmd8,cmd9,cmd10,cmd11,cmd12,cmd13,cmd14,cmd15,cmd16,cmd17,cmd18]
    for i in command:
        os.system(i)
def noncoding_mutation(input1,input2,input3,input4):
    con = input1.split('_')[0]
    case = input3.split('_')[0]
    combine=con+'_'+case
    cmd1='gatk Mutect2 -R ./reference/hg38/hg38.fa -I ./rna_result/'+case+'_recal.bam -tumor '+case+' -I ./rna_result/'+con+'_recal.bam -normal '+con+' -O ./mut_result/'+combine+'_1.vcf'
    cmd2 = 'gatk FilterMutectCalls -R ./reference/hg38/hg38.fa -V ./mut_result/'+combine+'_1.vcf -O ./mut_result/'+combine+'.vcf'
    cmd3 ="""cat ./mut_result/"""+combine+""".vcf | perl -alne '{if(/^#/){print}else{next unless $F[6] eq "PASS";next if $F[0] =~/_/;print } }' > ./mut_result/"""+combine+"""_filter_somatic.vcf"""
    cmd4 = 'perl ./biosoft/annovar/table_annovar.pl ./mut_result/'+combine+'_filter_somatic.vcf ./biosoft/annovar/humandb -buildver hg38 -out ./mut_result/'+combine+'_anno_out -remove -protocol refGene,cytoBand,exac03,avsnp150,dbnsfp30a,ALL.sites.2015_08,AFR.sites.2015_08,AMR.sites.2015_08,EUR.sites.2015_08,EAS.sites.2015_08 -operation g,r,f,f,f,f,f,f,f,f -nastring . -vcfinput'
    command = [cmd1, cmd2, cmd3, cmd4]
    for i in command:
        os.system(i)
def HLA_allele(input1,input2):
    case = input1.split('_')[0]
    cmd1='razers3 -i 95 -m 1 -dr 0 -o ./HLAtype/'+case+'_fished_1.bam ./biosoft/OptiType-1.3.5/data/hla_reference_rna.fasta ./rna_result/'+case+'_cut_R1.fastq.gz'
    cmd2='razers3 -i 95 -m 1 -dr 0 -o ./HLAtype/'+case+'_fished_2.bam ./biosoft/OptiType-1.3.5/data/hla_reference_rna.fasta ./rna_result/'+case+'_cut_R2.fastq.gz'
    cmd3='samtools fastq ./HLAtype/'+case+'_fished_1.bam > ./HLAtype/'+case+'_fished_1.fastq'
    cmd4='samtools fastq ./HLAtype/'+case+'_fished_2.bam > ./HLAtype/'+case+'_fished_2.fastq'
    cmd5='python ./biosoft/OptiType-1.3.5/OptiTypePipeline.py -i ./HLAtype/'+case+'_fished_1.fastq ./HLAtype/'+case+'_fished_2.fastq -r -o ./HLAtype/'+case+'_resultsDir -v'
    command = [cmd1, cmd2, cmd3, cmd4, cmd5]
    for i in command:
        os.system(i)
if __name__ == '__main__':
    con_input1 = sys.argv[1]
    con_input2 = sys.argv[2]     
    case_input1 = sys.argv[3]
    case_input2 =  sys.argv[4]     
    rnaseq_processing(con_input1,con_input2)
    rnaseq_processing(case_input1,case_input2)
    noncoding_mutation(con_input1,con_input2,case_input1,case_input2)
    HLA_allele(case_input1, case_input2)
    print('Model1 Run Successful!')
    ### python model1_rnaseq_mutation_hla.py con_R1.fastq.gz con_R2.fastq.gz case_R1.fastq.gz case_R2.fastq.gz


