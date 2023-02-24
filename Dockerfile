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
RUN conda install -c "conda-forge/label/gcc7" r-base
RUN conda install -c bioconda bwa samtools gatk4 picard bedtools blast
RUN pip install --default-timeout=100 numpy pandas future pyomo matplotlib -i https://pypi.tuna.tsinghua.edu.cn/simple

#Download
RUN wget -c http://119.3.70.71/PGNneo/data/reference.tar.gz && \
    tar -xzf reference.tar.gz && \
    rm reference.tar.gz

RUN wget -c http://118.31.70.55/ProGeo-neo/data/dbsnp_146.tar.gz && \
    tar -xzf dbsnp_146.tar.gz -C /home/PGNneo/reference && \
    rm dbsnp_146.tar.gz

RUN wget -c http://119.3.70.71/PGNneo/data/test.tar.gz && \
    tar -xzf test.tar.gz && \
    mv ./test/rnaseq ./ && \
    mv ./test/ms ./ && \
    rm test.tar.gz


WORKDIR /home/PGNneo/

CMD [ "/bin/bash" ]



