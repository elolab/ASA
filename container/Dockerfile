FROM ubuntu:20.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y apt-utils \
    && apt-get install -y dialog

RUN apt-get update \
    && apt-get -y install wget\
    && apt-get -y install git\
    && apt-get -y install emacs-nox\
    && apt-get -y install build-essential\
    && apt-get -y install clang\
    && apt-get -y install tmux\
    && apt-get -y install g++\
    && apt-get -y install zlib1g-dev\
    && apt-get -y install libghc-bzlib-dev\
    && apt-get -y install unzip\
    && apt-get -y install locales\
    && apt-get -y install sudo \
    && apt-get -y install mlocate \
    && apt-get -y install cmake \
    && apt-get -y install time \
    && apt-get -y install python \
    && apt-get -y install python3-biopython

RUN locale-gen en_US.UTF-8 en fi_FI.UTF-8

RUN echo hello
RUN git clone https://github.com/elolab/ASA.git /opt/asa

WORKDIR /opt/asa/
RUN git clone --branch 1.2 https://github.com/samtools/htslib.git
RUN git clone --branch 0.1.19 https://github.com/samtools/samtools.git samtools-0.1.19
RUN git clone --branch 1.2 https://github.com/samtools/samtools.git samtools-1.2


COPY bowtie2-2.2.3-with-n-award.diff /root/
COPY bowtie2-2.2.3-source.zip /root/

RUN unzip /root/bowtie2-2.2.3-source.zip
WORKDIR /opt/asa/bowtie2-2.2.3
RUN cat /root/bowtie2-2.2.3-with-n-award.diff | patch -p 1
RUN make -j 16

WORKDIR /opt/asa/htslib
RUN make

WORKDIR /opt/asa/samtools-0.1.19
RUN make

WORKDIR /opt/asa/samtools-1.2
RUN make

WORKDIR /opt/asa/
RUN ln -s bowtie2-2.2.3/bowtie2 bowtie2
RUN ln -s bowtie2-2.2.3/bowtie2-build bowtie2-build
RUN ln -s samtools-0.1.19/samtools samtools-git-0.1.19
RUN ln -s samtools-0.1.19/bcftools/bcftools bcftools-0.1.19
RUN ln -s samtools-0.1.19/bcftools/vcfutils.pl vcfutils-0.1.19.pl
RUN ln -s samtools-1.2/samtools samtools-git

RUN chmod u+x phylotool.py aligner-asa.py filter_database.py identify-v4-16S-region.py trees_from_families.sh 

#RUN ln -s /hugedata/src/asa/phylotool.py .
#RUN ln -s /hugedata/src/asa/aligner-asa.py .
#RUN ln -s /hugedata/src/asa/filter_database.py .
#RUN ln -s /hugedata/src/asa/identify-v4-16S-region.py .
#RUN ln -s /hugedata/src/asa/trees_from_families.sh .
#RUN ln -s /hugedata/src/asa/create_con_ref_tree.R .


