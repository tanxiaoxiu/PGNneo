FROM xiaoxiutan/pgnneo:v1

LABEL maintainer="tanxiaoxiu" email="tanxiaoxiu@sjtu.edu.cn"

ENV PATH /root/miniconda3/bin:$PATH

RUN yum install -y java-1.8.0-openjdk
RUN conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free/
RUN conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main/
RUN conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/
RUN conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/msys2/
RUN conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/pytorch/
RUN conda config --add channels https://mirrors.ustc.edu.cn/anaconda/cloud/bioconda/
RUN conda config --set show_channel_urls yes
RUN conda install -c r-base=4.0.3
RUN conda install -c bioconda bwa samtools gatk4 picard bedtools blast
RUN pip install --default-timeout=100 numpy pandas future pyomo matplotlib -i https://pypi.tuna.tsinghua.edu.cn/simple

#Download
RUN wget -c http://119.3.70.71/PGNneo/data/reference.tar.gz && \
    tar -xzf reference.tar.gz -C /home/PGNneo && \ 
    rm reference.tar.gz

RUN wget -c http://118.31.70.55/ProGeo-neo/data/dbsnp_146.tar.gz && \
    tar -xzf dbsnp_146.tar.gz && \
    mv ./dbsnp_146/dbsnp_146.hg38.vcf /home/PGNneo/reference && \
    mv ./dbsnp_146/dbsnp_146.hg38.vcf.idx /home/PGNneo/reference && \
    rm dbsnp_146.tar.gz

RUN wget -c http://119.3.70.71/PGNneo/data/test.tar.gz && \
    tar -xzf test.tar.gz && \
    mv ./test/rnaseq /home/PGNneo && \
    mv ./test/ms /home/PGNneo && \
    rm test.tar.gz


WORKDIR /home/PGNneo/

CMD [ "/bin/bash" ]



