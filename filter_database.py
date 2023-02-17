#!/usr/bin/env python

import sys
import os.path
import Bio
from Bio import SeqIO
from Bio import GenBank
from Bio.Seq import Seq
from optparse import make_option
from cogent.util.misc import parse_command_line_parameters

__author__ = "Sami Pietila"
__copyright__ = "Copyright 2014"
__credits__ = ["Sami Pietila"]
__license__ = "GPL"
__version__ = "0.5-preview"
__maintainer__ = "Sami Pietila"
__email__ = "sampie@iki.fi"
__status__ = "Preview"

script_info={}
script_info['brief_description']="""Convert RDP unaligned fasta into ASA format"""
script_info['script_description']="""Convert RDP unaligned fasta into ASA format"""

script_info['script_usage']=[]

script_info['script_usage'].append(("""Standard Example usage:""","""""","""%prog\n\t--output_gb rdp_11_5_T_I_L1200.gb\n\t-f rdp_11_5_T_I_L1200.fasta\n\t--rdp_gb current_Bacteria_unaligned.gb\n\t--rdp_fasta current_Bacteria_unaligned.fa\n\t-t rdp_11_5_T_I_L1200_taxonomy.txt\n\t-s rdp_11_5_T_I_L1200_seqs.fa\n\t--filter_min_seqlen 1200\n\t--filter_db_typestrain\n\t--filter_db_isolates"""))

script_info['output_description']="""ASA taxonomy and sequence files."""


script_info['required_options']=[\
    # input fasta filepath(s)
    make_option(None, '--rdp_gb',
        help='RDP genbank file.'),

    make_option(None, '--rdp_fasta',
        help='Read taxonomic lineage from RDP FASTA file.'),

    make_option('-t', '--taxonomy',
        help='File name for taxonomy mapping to be written.'),

    make_option('-s', '--seqs',
         help='File name for sequences to be written.'), 
                                 
    ]
                    
script_info['optional_options']=[\
    make_option('-f', '--fasta',
        help='Name for FASTA file with taxonomy and sequences to be written.'),

    make_option(None, '--output_gb',
        help='Name for Genbank file to be written.'),

    make_option(None, '--lineage_level', default="species", 
        help='Write rank until given lineage level [domain,phylum,class,subclass,order,suborder,family,genus,species].'),

    make_option(None, '--filter_by_accession',
        help='A file containing list of accession numbers that are used in filtering.'),

    make_option(None, '--allow_only_first_match',action="store_true", default=False,
        help='When filtering by accession or species names, only the first match is allowed per accession/species'),

    make_option(None, '--filter_lowqual',action="store_true", default=False,
        help='Filter out low quality sequences that are marked in RDP with lowqual.'),

    make_option(None, '--filter_min_seqlen', default=0,
        help='Filter out RDP sequences below given threshold. 1200 is considered as a good sequence length'),

    make_option(None, '--filter_db_typestrain',action="store_true", default=False,
        help='Accept only type strains.'),

    make_option(None, '--filter_db_isolates',action="store_true", default=False,
        help='Accept only isolates.'),

     make_option(None, '--write_only_last_rank',action="store_true", default=False,
         help='Write only last rank instead of the whole lineage.'),
                                 
    make_option(None, '--separate_fasta',action="store_true", default=False,
        help='Split RDP entries into separate FASTA files, each entry in its own file.'),

    make_option(None, '--filter_species',
        help='Accept only listed species (comma separated, no spaces).')
    ]
    

script_info['version'] = __version__

# domain (d)
# phylum (p)
# class (c)
# sub class (b)
# order (o)
# sub order (u)
# family (f)
# genus (g)
# species (s)

taxonomy_prefix_list = "r__", "d__", "k__", "p__", "c__", "b__", "o__", "u__","f__", "g__", "s__" 

rankc={}
rankc['rootrank']='r__'
rankc['domain']='d__'
rankc['phylum']='p__'
rankc['class']='c__'
rankc['subclass']='b__'
rankc['order']='o__'
rankc['suborder']='u__'
rankc['family']='f__'
rankc['genus']='g__'
rankc['species']='s__'


def read_lineage(faFileName):

    taxa_ranks_by_id = {}
    records = SeqIO.parse(faFileName, "fasta")

    for record in records:

        ls = record.description.split("\t")[0].rstrip(";")[len(record.id)+1:].split(";")
        ls = ls[:2]
        for i in range(len(ls)):
            ls[i] = ls[i].strip()
            ls[i] = ls[i].replace(" ", "_")
        species = "_".join(ls)    
        
        
        taxl = []
        lineage = record.description.split("\t")[1].split("=")[1].strip().rstrip(";").split(";")
        for i in range(0, len(lineage)-1, 2):
            rank_label = lineage[i+1]
            rank = lineage[i].strip().strip('"')
            taxl.append(rankc[rank_label]+rank)
        taxl.append("s__"+species)
        taxa_ranks_by_id[record.id] = taxl
        
    return taxa_ranks_by_id
    

def filter_gb_by_typestrain(records):

    filtered_records = []

    for record in records:
        if "(T)" in record.description:
            filtered_records.append(record)
    
    return filtered_records


def filter_gb_by_accession_codes(records, accession_codes, allow_only_first_match):

    filtered_records = []
    acs = []
    
    for record in records:
        accession_code = record.annotations["accessions"][0]
        if accession_code in accession_codes:
            if accession_code not in acs:
                filtered_records.append(record)
                acs.append(accession_code)
            elif not allow_only_first_match:
                filtered_records.append(record)

    return filtered_records, acs


def filter_gb_by_isolate(records):
    filtered_records = []

    for record in records:

        if not "uncultured" in record.annotations["source"]:
            filtered_records.append(record)
    
    return filtered_records
    

def filter_gb_by_seq_length(records, min_seq_length):

    filtered_records = []
    
    for record in records:

        seq_length = None

        for f in record.features:

            if f.type == "rRNA":
                if seq_length is None and f.location != False:
                    seq_length = len(f.location)
                else:
                    print record.name
                    raise Exception("Value error", "GB records contains multiple sequence locations")
        
        if int(seq_length) >= int(min_seq_length):
            filtered_records.append(record)

    return filtered_records


def filter_gb_by_qual(records):

    filtered_records = []
    
    for record in records:

        is_highqual = True
        if record.annotations.has_key("comment") and "islowqual: true" in record.annotations["comment"]:
            is_highqual = False

        if is_highqual is True:
            filtered_records.append(record)

    return filtered_records


def filter_gb_by_species(records, species_names, allow_only_first_match):

    found_species = []
    filtered_records = []
    
    for record in records:
        
        assert lineages_by_id.has_key(record.name)
        species = lineages_by_id[record.name][-1].strip()
        
        for a_species in species_names:
            if "s__"+a_species in species:
                if a_species not in found_species:
                    filtered_records.append(record)
                    found_species.append(a_species)
                elif not allow_only_first_match:
                    filtered_records.append(record)
                    
    return records, found_species


def make_db(records, lineages_by_id):
    db = []
    
    for record in records:
   
        strain = None
        for f in record.features:

            if f.type == "source":
                if f.qualifiers.has_key("strain"):
                    if strain is None:
                        strain = f.qualifiers["strain"][0]
                    else:
                        raise Exception("Value error", "GB records contains multiple strain definitions for a record")

        assert lineages_by_id.has_key(record.name)
        lineageL = lineages_by_id[record.name]
        lineage = ";".join(lineageL)
        
        db.append((lineage, record.seq, strain))

    return db


def print_db(db):
    for t in db:
        print t[0], t[6]


def filter_rank(db, rank_level, only_last):
    f_db = []
    rank_prefix = rankc[rank_level]
    for t in db:
        taxonomy = t[0].split(";")
        seq = t[1]
        for i, r in enumerate(taxonomy):
            if rank_prefix == r.strip()[:3]:
                if only_last:
                    f_db.append((r.strip(), seq))
                else:
                    f_db.append((";".join(taxonomy[:i+1]), seq))
    return f_db

# def strip_taxonomy(db):
#     f_db = []
#     for t in db:
#         taxonomy = t[0]
#         seq = t[1]
#         a_species_name = taxonomy.split(";")[-1].strip()
#         f_db.append((a_species_name, seq))
#     return f_db


def write_db(db, seqsFN, taxaFN):
    seqsFH = open(seqsFN, "w")
    taxaFH = open(taxaFN, "w")
    
    for i, entry in enumerate(db):
        taxonomy = entry[0]
        seq = entry[1]
        seqsFH.write(">"+str(i) + "\n")
        seqsFH.write(str(seq)+"\n")
        taxaFH.write(str(i)+"\t"+taxonomy+"\n")
    
    seqsFH.close()
    taxaFH.close()


def write_fasta(db, fastaFN):
    fastaFH = open(fastaFN, "w")
    
    for i, entry in enumerate(db):
        taxonomy = entry[0]
        seq = entry[1]
        fastaFH.write(">"+ taxonomy + "\n")
        fastaFH.write(str(seq)+"\n")
    
    fastaFH.close()

def write_gb(records, output_gb):
    fh = open(output_gb, "w")
    SeqIO.write(records, fh, "gb")
    fh.close()
    return

def write_separate_fasta(db):
    uniqs = {}
    
    for entry in db:
        taxonomy = entry[0]
        seq = entry[1]
        
        if uniqs.has_key(taxonomy) == True:
            uniqs[taxonomy] += 1
            uq = uniqs[taxonomy]
            fname = taxonomy +"_"+str(uq) + ".fasta"
        else:
            uniqs[taxonomy] = 0
            fname = taxonomy+".fasta"
        
        fastaFH = open(fname, "w")
        fastaFH.write(">"+ taxonomy + "\n")
        fastaFH.write(str(seq)+"\n")
        fastaFH.close()

    
if __name__=="__main__":

    option_parser, opts, args = parse_command_line_parameters(**script_info)

    lineages_by_id = None
    if opts.rdp_fasta != None:
        lineages_by_id = read_lineage(opts.rdp_fasta)
    
    records = SeqIO.parse(opts.rdp_gb, "gb")

    if opts.filter_db_typestrain == True:
        records = filter_gb_by_typestrain(records)

    if opts.filter_db_isolates == True:
        records = filter_gb_by_isolate(records)

    if int(opts.filter_min_seqlen) > 0:
        records = filter_gb_by_seq_length(records, opts.filter_min_seqlen)

    if opts.filter_lowqual == True:
        records = filter_gb_by_qual(records)

    if opts.filter_by_accession is not None:
        accessionfh = open(opts.filter_by_accession, "r")
        accessionlist = accessionfh.read().split("\n")
        accessionfh.close()
        query_acs = [x.strip() for x in accessionlist]
        records,acs = filter_gb_by_accession_codes(records, query_acs, opts.allow_only_first_match)
        missing_acs = set(query_acs).difference(set(acs))
        if len(missing_acs)>0:
            print "Warning: Following accession codes were not found from database: " + ",".join(missing_acs)
    
    if opts.filter_species != None:
        query_species = opts.filter_species.split(",")
        records, found_species = filter_gb_by_species(records, query_species, opts.allow_only_first_match)
        missing_species = set(query_species).difference(set(found_species))
        if len(missing_species)>0:
            print "Warning: Following species were not found from database: " + ",".join(missing_species)

    db = make_db(records, lineages_by_id)
    db = filter_rank(db, opts.lineage_level, opts.write_only_last_rank)
        
    write_db(db, opts.seqs, opts.taxonomy)
    
    if opts.output_gb != None:
        write_gb(records, opts.output_gb)
    
    if opts.fasta != None:
        write_fasta(db, opts.fasta)
        
    if opts.separate_fasta == True:
        write_separate_fasta(db)
        
        
