#biosoft
wget -c http://119.3.70.71/PGNneo/data/biosoft.tar.gz
tar -zxvf biosoft.tar.gz

#reference
wget -c http://119.3.70.71/PGNneo/data/reference.tar.gz
wget -c http://118.31.70.55/ProGeo-neo/data/dbsnp_146.tar.gz
tar -zxvf reference.tar.gz
tar -zxvf dbsnp_146.tar.gz
mv dbsnp_146/dbsnp_146.hg38.vcf reference/
mv dbsnp_146/dbsnp_146.hg38.vcf.idx reference/

#test
wget -c http://119.3.70.71/PGNneo/data/test.tar.gz
tar -zxvf test.tar.gz
wget -c https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos3/sra-pub-run-20/SRR8382658/SRR8382658.2
wget -c https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos3/sra-pub-run-20/SRR8382623/SRR8382623.2
./biosoft/sratoolkit/bin/fastq-dump SRR8382658.2 --split-3 --gzip -O ./test/rnaseq
./biosoft/sratoolkit/bin/fastq-dump SRR8382623.2 --split-3 --gzip -O ./test/rnaseq
mv ./test/rnaseq/SRR8382658.2_1.fastq.gz ./test/rnaseq/case_R1.fastq.gz
mv ./test/rnaseq/SRR8382658.2_2.fastq.gz ./test/rnaseq/case_R2.fastq.gz
mv ./test/rnaseq/SRR8382623.2_1.fastq.gz ./test/rnaseq/con_R1.fastq.gz
mv ./test/rnaseq/SRR8382623.2_2.fastq.gz ./test/rnaseq/con_R2.fastq.gz
mv ./test/rnaseq ./
mv ./test/ms ./

rm -r test
rm -r dbsnp_146
rm -f *.2
rm -f *.gz
