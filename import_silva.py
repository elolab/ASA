#! /usr/bin/env python3
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import argparse
import shlex
import subprocess
import os
import io
import sys

__author__ = "Sami Pietila"
__copyright__ = "Copyright 2024"
__credits__ = ["Sami Pietila"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Sami Pietila"
__email__ = "sampie@iki.fi"

__DEFAULT_SILVA_VERSION__ = "138.1"

class NonZeroReturnValueException(Exception):
    def __init__(self, returnvalue, msg):
        self.msg = msg
        self.returnvalue = returnvalue
        return

if __name__=="__main__":
    
    parser = argparse.ArgumentParser(description='This is a script to import SILVA database to ASA')

    parser.add_argument('--silva-version', 
                        action='store', 
                        dest='silva_version', 
                        required=False, 
                        default=__DEFAULT_SILVA_VERSION__,
                        help='Silva version to be imported. [Default: ' + __DEFAULT_SILVA_VERSION__ + ' ]')


    args = parser.parse_args()

    cwd = "./"

    input_fasta = "SILVA_" + args.silva_version + "_SSURef_NR99_tax_silva.fasta"

    silva_gz = "SILVA_" + args.silva_version + "_SSURef_NR99_tax_silva.fasta.gz"
    silva_fasta = "SILVA_" + args.silva_version + "_SSURef_NR99_tax_silva.fasta"

    if not (os.path.exists(silva_gz) or os.path.exists(silva_fasta)): 
        print ("Downloading SILVA ", args.silva_version, "....")

        cmd = shlex.split(
            "wget https://www.arb-silva.de/fileadmin/silva_databases/release_138_1/Exports/SILVA_" + args.silva_version + "_SSURef_NR99_tax_silva.fasta.gz")

        proc = subprocess.Popen(
            cmd, cwd=cwd)

        proc.wait()

        if proc.returncode != 0:
            raise NonZeroReturnValueException(proc.returncode, 'SILVA Download failed. Is wget installed?')

    if not os.path.exists(os.path.join(cwd, silva_fasta)): 

        print ("Unpacking SILVA...")

        cmd = shlex.split(
            "gzip -d " + silva_gz )
        
        proc = subprocess.Popen(
            cmd, cwd=cwd)
        
        proc.wait()

        if proc.returncode != 0:
            raise NonZeroReturnValueException(proc.returncode, 'SILVA uncompress failed. Is gzip installed?')

    print ("Converting " + input_fasta + " into ASA fasta and taxonomy.")

    out_records = []
    in_records = SeqIO.parse(os.path.join(cwd, input_fasta), "fasta")
    taxonomy_mapping = []

    for in_record in in_records:

        ID = in_record.id

        dna_sequence = in_record.seq.back_transcribe()
        descL = in_record.description.split(" ", 1)
        taxonomy_string = descL[1]
        assert (descL[0] == in_record.id)

        out_record = SeqRecord(Seq(dna_sequence), id=ID, description="")
        out_records.append(out_record)

        taxonomy_mapping.append((ID, taxonomy_string))


    ASA_taxonomy_file = "silva_" + args.silva_version + "_taxonomy.txt"

    print("Writing " + ASA_taxonomy_file)

    with open(os.path.join(cwd, ASA_taxonomy_file), "w") as fh:
        for entry in taxonomy_mapping:
            fh.write(entry[0] + "\t" + entry [1] + "\n")

    ASA_fasta_file = "silva_" + args.silva_version + ".fasta"

    print("Writing " + ASA_fasta_file)

    with open(os.path.join(cwd, ASA_fasta_file), "w") as fh:
        SeqIO.write(out_records, fh, "fasta")
