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
RUN wget -c http://119.3.70.71/PGNneo/data/test.tar.gz && \
    tar -xzf test.tar.gz && \
    mv ./test/rnaseq ./ && \
    mv ./test/ms ./ && \

rm -r test
rm -r dbsnp_146
rm -f *.gz
