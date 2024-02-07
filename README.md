# Adaptive Sequence Alignment (ASA)

## Software
The analysis flow logic is implemented in aligner-tg program. The program requires Bowtie2 and Samtools utilities as well as RDP II database. The software is designed to run under Linux operating system.  

### Bowtie2

Download bowtie2 version 2.2.3: 

```
wget \
https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.3/bowtie2-2.2.3-source.zip
```
Unzip the package.

```
unzip bowtie2-2.2.3-source.zip
```

Next, patch Bowtie2 with following command:

```
cd bowtie2-2.2.3
cat bowtie2-2.2.3-with-n-award.diff | patch -p 1
cd ..

ln -s bowtie2-2.2.3/bowtie2 bowtie2
ln -s bowtie2-2.2.3/bowtie2-build bowtie2-build
```

Finally, compile Bowtie2:

```
make
```

### Samtools


Compile htslib 1.2

```
git clone https://github.com/samtools/htslib.git
cd htslib
make
git checkout 1.2
cd ..
```

Compile Samtools 0.1.19

```

git clone https://github.com/samtools/samtools.git samtools-0.1.19
cd samtools-0.1.19
git checkout 0.1.19
make
cd ..

ln -s samtools-0.1.19/samtools samtools-git-0.1.19
ln -s samtools-0.1.19/bcftools/bcftools bcftools-0.1.19
ln -s samtools-0.1.19/bcftools/vcfutils.pl vcfutils-0.1.19.pl
```

Compile Samtools 1.2

```
git clone https://github.com/samtools/samtools.git samtools-1.2
cd samtools-1.2
git checkout 1.2
make
cd ..

ln -s  samtools-1.2/samtools samtools-git

```

### RDP database filtering

Download and unzip Bacterial 16S sequences (as Genbank and FASTA) from RDP. The files are: current_Bacteria_unaligned.fa.gz and current_Bacteria_unaligned.gb.gz.

```
cd ..

wget \
https://rdp.cme.msu.edu/download/current_Bacteria_unaligned.gb.gz \
--no-check-certificate

wget \
https://rdp.cme.msu.edu/download/current_Bacteria_unaligned.fa.gz \
--no-check-certificate
```

Next, convert and filter the data to the format required by ASA:

```
gunzip current_Bacteria_unaligned.gb.gz current_Bacteria_unaligned.fa.gz

./filter_database.py \
--output_gb rdp_11_5_T_I_L1200.gb \
-f rdp_11_5_T_I_L1200.fasta \
--rdp_gb current_Bacteria_unaligned.gb \
--rdp_fasta current_Bacteria_unaligned.fa \
-t rdp_11_5_T_I_L1200_taxonomy.txt \
-s rdp_11_5_T_I_L1200_seqs.fa \
--filter_min_seqlen 1200 \
--filter_db_typestrain \
--filter_db_isolates
```
### MBARC26_public

The pre-processed MBARC-26 dataset is available in the NCBI SRA archive with the accession number SRX1836716. It can be downloaded using the SRA toolkit



### ASA tool example run

Create a directory for ASA to store temporary items in

```
mkdir ./tmp
mkdir ./out

./aligner-tg.py \
-i R1.fastq.gz,R2.fastq.gz \
-t rdp_11_5_T_I_L1200_taxonomy.txt \
--db_fasta rdp_11_5_T_I_L1200_seqs.fa \
--top_refs 2 \
--rounds 2 \
--best_match_only \
--candlist_bt2_score_min 0.99 \
--top_refs_file top_refs-example.txt \
--aligner_threads 6 \
--embed_in_n \
-o out/consensus.txt
```

## Output

First, the tool aligns the reads to the RDP database in order to discover the best aligning references. The result can be found in  `top_refs-example.txt`

Next, the tool chooses a number of the best matching references (specified by `--top_refs`) and starts iteratively building each sequence. The maximum amount of iterations can be set by the `--rounds` parameter

### Creating the RDP reference database

First, the RDP database was split into subsets on family level in order to decrease computational load. 

Next, one tree per family was generated using clustalw and fasttree. Finally, the minimum distances between genera were calculated using the `phylotool.py` script. It takes a directory of trees and creates a pickled object (`min_family_dists.pickle`) containing the cross-genera cutoff values.


### Consensus-reference mapping 

The consensus-reference tree is generated using the `create_con_ref_tree.R` script. Modify the paths in this file to point to the correct locations


`phylotool.py` processes the consensus-reference tree and reports the species mapping. Furthermore, a note is made if the mapping is within the cross-genera cutoff. If the cutoff is exceeded, then the mapping is not too confdent, since it is possibly mapped to the wrong genera.
