# Adaptive Sequence Alignment (ASA)

## Software
The analysis flow logic is implemented in aligner-asa program. The program requires Bowtie2 and Samtools utilities as well as RDP II database. The software is designed to run under Linux operating system.  

### MBARC26_public

The pre-processed MBARC-26 dataset is available in the NCBI SRA archive with the accession number SRX1836716. It can be downloaded using the SRA toolkit

### ASA tool installation

```
docker pull elolab/asa
```
Once download installation completes, installation can be checked:

```
docker run -it --rm elolab/asa /opt/asa/aligner-asa.py --help
```

### ASA tool example run

The ASA container can be run with following command. In this example the computer has /data -folder, which is visible also to the container. The output is written to /data/out/ folder, which is required to exist.

```

docker run -it --rm -v /data:/data elolab/asa /opt/asa/aligner-asa.py \
-i /data/R1.fastq.gz,R2.fastq.gz \
-t /data/pro_rdp_11_5_T_I_L1200_taxonomy.txt \
--db_fasta /data/pro_rdp_11_5_T_I_L1200_seqs.fa \
--top_refs 2 \
--rounds 2 \
--best_match_only \
--candlist_bt2_score_min 0.99 \
--top_refs_file top_refs-example.txt \
--aligner_threads 6 \
--embed_in_n \
-o /data/out/consensus.txt
```

## Output

First, the tool aligns the reads to the RDP database in order to discover the best aligning references. The result can be found in  `top_refs-example.txt`

Next, the tool chooses a number of the best matching references (specified by `--top_refs`) and starts iteratively building each sequence. The maximum amount of iterations can be set by the `--rounds` parameter

### Creating the RDP reference database

First, the RDP database was split into subsets on family level in order to decrease computational load. 

Next, one tree per family was generated using clustalw and fasttree. Finally, the minimum distances between genera were calculated using the `phylotool.py` script. It takes a directory of trees and creates a pickled object (`min_family_dists.pickle`) containing the cross-genera cutoff values.


### Consensus-reference mapping 

The consensus-reference tree is generated using the `/opt/asa/create_con_ref_tree.R` script. Modify the paths in this file to point to the correct locations


`/opt/asa/phylotool.py` processes the consensus-reference tree and reports the species mapping. Furthermore, a note is made if the mapping is within the cross-genera cutoff. If the cutoff is exceeded, then the mapping is not too confdent, since it is possibly mapped to the wrong genera.


## Appendix

#### Creating a new reference database files

ASA can use any fasta file as reference file. However, if there is taxonomy associated with fasta sequences, there are two files: 1) fasta file and 2) related taxonomy file.

Fasta file example:
```
>0
AGAGTTTGATCATGGCTCAGGACGA
>1
TATTTCGACGGGTTCCG
>2
TAAAGAGCTCG
```
Corresponding taxonomy file example:
```
0[TAB]r__Root;d__Bacteria;p__Actinobacteria;c__Actinobacteria;b__Acidimicrobidae;o__Acidimicrobiales;u__Acidimicrobineae;f__Acidimicrobiaceae;g__Acidimicrobium;s__Acidimicrobium_ferrooxidans_(T)_ICP
1[TAB]r__Root;d__Bacteria;p__Actinobacteria;c__Actinobacteria;b__Acidimicrobidae;o__Acidimicrobiales;u__Acidimicrobineae;f__Acidimicrobiaceae;g__Ferrimicrobium;s__Ferrimicrobium_acidiphilum_(T)_T23
2[TAB]r__Root;d__Bacteria;p__Actinobacteria;c__Actinobacteria;b__Acidimicrobidae;o__Acidimicrobiales;u__Acidimicrobineae;f__Acidimicrobiaceae;g__Ferrithrix;s__Ferrithrix_thermotolerans_(T)_Y005
```
Where [TAB] is the tabulator character.

### Example, RDP database filtering

Download and unzip Bacterial 16S sequences (as Genbank and FASTA) from RDP. The files are: current_Bacteria_unaligned.fa.gz and current_Bacteria_unaligned.gb.gz. Note. RDP II project website seems to be available currently.

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
