#!/usr/bin/env Rscript

library(Biostrings)


# Needed files and directories:
work.dir <- "dir"
fasta.out <- "con_ref.fasta"
tree.out <- "con_ref.tree"
ref.file.in <- "MBARC26_reflist.txt"
cand.file.in <- "consensus_sequences/pro_top150_ids.txt"
rdp.seqs.file.in <- "RDP/pro_rdp_11_5_T_I_L1200_seqs.fa"
rdp.tax.file.in <- "RDP/pro_rdp_11_5_T_I_L1200_taxonomy.txt"

# % of N allowed in consensus sequences
n.cutoff <- 0.03

# Auto-generated intermediate files
phylip.out <- gsub(fasta.out, pattern = ".fasta",fixed = T,
                    replacement = ".phy")

setwd(work.dir)

# Load in both RDP sequences and taxonomy strings together with the reference
# file. Required for generation of the tree and marking species from the
# reference sequences
rdp.seqs <- readDNAStringSet(rdp.seqs.file.in, format = "fasta")

rdp.tax <- read.delim(rdp.tax.file.in, sep = "\t", quote = "",
                      stringsAsFactors = F, header = F)

cand.ids <- read.table(cand.file.in)

# Change the column names to something more comprehensible
names(rdp.tax) <- c("id", "tax")
names(cand.ids) <- c("id")


# Get all the candidate ID sequences directly from the RDP database
# and write them to file
writeXStringSet(rdp.seqs[cand.ids$id], filepath = fasta.out, format = "fasta")

# Loop over the consensus_sequences directory and get the last iterations of
# each consensus
for (file in dir(path = "consensus_sequences", pattern = "_consensus.txt", full.names = T)) {

    # Get last iteration from consensus sequence
    con.seq <- tail(readDNAStringSet(file),1)
    
    # Discard all consensus sequences with "N%" over the cutoff value
    if ((alphabetFrequency(con.seq)[,"N"] / width(con.seq)) <= n.cutoff) {
        
        # Change the consensus sequence FASTA header to ID_C
        con.seq.id <- unlist(strsplit(names(con.seq), split = "_"))[1]
        names(con.seq) <- paste0(con.seq.id,"_C")

        # Write the entry to the same output FASTA as the candidates
        writeXStringSet(con.seq, filepath = fasta.out, append = T,
                        format = "fasta")
    }
}

# Create the big candidate-consensus tree (con_ref tree) from the
# new fasta file

# Alignment
system2(command = "clustalw",
        args = paste(fasta.out, "-output=PHYLIP"),
        wait = T)

# Tree generation
system2(command = "fasttree",
        args = paste("-gtr -gamma -cat 4 -nt", phylip.out , ">", tree.out),
        wait = T)