#! /usr/bin/env python3
from Bio import SeqIO
import argparse
from primerlib import find_primer
from primerlib import revcompPrimer

__author__ = "Sami Pietila"
__copyright__ = "Copyright 2016"
__credits__ = ["Sami Pietila"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Sami Pietila"
__email__ = "sampie@iki.fi"

#V3+V4:
#Forward CCTACGGGNGGCWGCAG
#Reverse GACTACHVGGGTATCTAATCC

#V4:
#forward: GTGCCAGCMGCCGCGGTAA
#reverse: GGACTACHVGGGTWTCTAAT

#V4+V5
#forward: GTGCCAGCMGCCGCGGTAA
#reverse: CCCGTCAATTCMTTTRAGT

primers = {#("v3v4","forward"): "CCTACGGGNGGCWGCAG",
            ("v4","forward"): "GTGCCAGCMGCCGCGGTAA",
            #("v4v5","forward"): "GTGCCAGCMGCCGCGGTAA",
            #("v3v4","revcomp"):"GACTACHVGGGTATCTAATCC",
            ("v4","revcomp"):"GGACTACHVGGGTWTCTAAT",
            #("v4v5","revcomp"):"CCCGTCAATTCMTTTRAGT"
            }


def find_primer_from_fasta(input_fasta, max_missbases):

    records = []

    counts = {}

    for region, ptype in primers:
        counts[region, ptype] = {"found":0, "missed":0}
    
    with open(input_fasta, "r") as input_fh:
        
        for record in SeqIO.parse(input_fh, "fasta"):

            seq = str(record.seq)

            forward_primer_found = False
            reverse_primer_found = False
            
            for region, ptype in primers:

                if ptype == "revcomp":
                    primer = revcompPrimer(primers[region, ptype])
                elif ptype == "forward":
                    primer = primers[region, ptype]
                else:
                    raise Exception("Unknown primer treatment")
            
                _,inds = find_primer(primer, record.seq, int(max_missbases))


                offset = 0
                
                for ind in inds:
                    from_ind = None
                    to_ind = None

                    from_ind = ind+offset
                    to_ind = ind+offset + len(primer)

                    if from_ind != None or to_ind != None:
                        if ptype == "forward":
                            forward_primer_found = True
                            seq = seq[:from_ind] + "[" + seq[from_ind:]
                            offset+=1
                        if ptype == "revcomp":
                            reverse_primer_found = True
                            seq = seq[:to_ind] + "]" + seq[to_ind:]
                            offset+=1
                            
            if forward_primer_found and reverse_primer_found:
                print (">" + record.id)
                print (seq)
            
 #    if records:
 #       SeqIO.write(records, output_fasta, "fasta")

                    
    return

if __name__=="__main__":
    
    parser = argparse.ArgumentParser(description='This is a script to find primer sequences from fasta')

    parser.add_argument('--input-fasta', 
                        action='store', 
                        dest='input_fasta', 
                        required=True, 
                        default=None,
                        help='Input fasta file')

    parser.add_argument('--max_missbases', 
                        action='store', 
                        dest='max_missbases', 
                        required=False, 
                        default="1",
                        help='Maxum amount of wrong bases [default 1]')


    args = parser.parse_args()

    find_primer_from_fasta(args.input_fasta, args.max_missbases)

