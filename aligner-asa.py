#! /usr/bin/env python3

import sys
import os
import os.path
import shutil
import subprocess
import tempfile
import uuid
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import errno
import argparse

__author__ = "Medical Bioinformatics Centre"
__copyright__ = "Copyright 2024"
__credits__ = ["Sami Pietila", "Niklas Paulin"]
__license__ = "GPL"
__version__ = "3.0"
__maintainer__ = "Medical Bioinformatics Centre"
__email__ = "mbc@utu.fi"
__status__ = ""

install_path = "/opt/asa/"

def read_taxonomy(fname):
    taxFH = open(fname, "r")
    taxData = taxFH.read()
    taxFH.close()
    lines = taxData.split("\n")
    taxonomy = {}
    for l in lines:
        if l == "":
            continue
        ll = l.split("\t")
        taxonomy[ll[0]]=ll[1]
    return taxonomy

def read_taxonomy_from_fasta(fname):
    with open(fname, "r") as fh:
        lines = fh.read().split("\n")
    taxonomy = {}
    for l in lines:
        if l.strip().startswith(">"):
            taxon = l.strip()[1:].strip()
            assert taxon not in taxonomy
            taxonomy[taxon]=taxon
    return taxonomy
    

def writeFasta(fhandle, records):
    for record in records:
        fhandle.write(">"+str(record.id)+"\n")
        fhandle.write(str(record.seq)+"\n")

def readFasta_list(fastaFN):
    seqs = []
    handle = open(fastaFN, "r")
    for record in SeqIO.parse(handle, "fasta") :
        newkey = record.description.split(" ")[0]
        seqs.append(newkey, str(record.seq))
    handle.close()    
    return seqs

def readFasta_dict(fastaFN):
    seqs = {}
    handle = open(fastaFN, "r")
    for record in SeqIO.parse(handle, "fasta"):
        newkey = record.description.split(" ")[0]
        seqs[newkey]=str(record.seq)
    handle.close()    
    return seqs

def readFastaFH_dict(handle):
    seqs = {}
    for record in SeqIO.parse(handle, "fasta"):
        newkey = record.description.split(" ")[0]
        seqs[newkey]=str(record.seq)
    return seqs

def readFastq_dict(fastqFN):
    seqs = {}
    handle = open(fastqFN, "r")
    for record in SeqIO.parse(handle, "fastq"):
        newkey = record.description.split(" ")[0]
        seqs[newkey]=str(record.seq)
    handle.close()
    
    for key in list(seqs.keys()):
        seq = seqs.pop(key)
        newkey = key.split(" ")[0]
        seqs[newkey] = seq
    return seqs


def readFastqFH_dict(handle):
    seqs = {}
    for record in SeqIO.parse(handle, "fastq"):
        newkey = record.description.split(" ")[0]
        seqs[newkey]=str(record.seq)
    
    for key in list(seqs.keys()):
        seq = seqs.pop(key)
        newkey = key.split(" ")[0]
        seqs[newkey] = seq
    return seqs


def aggregate(match_fh, paired_end, taxonomy, bt2_score):
    ref_matches = {}

    if paired_end:
        
        prev_match = None
        for line in match_fh:
            match = line.strip().split("\t")

            # Skip header lines                
            if match[0] == "@SQ" or match[0] == "@PG" or match[0] == "@HD":
                continue
    
            flags = int(match[1])
            if (flags & 0b10000011 == 0b10000011): # second pair of a properly matching pair

                queryID = prev_match[0]
                assert queryID == match[0]
                refID = prev_match[2]
                assert refID == match[2]
            
                if refID in ref_matches:
                    ref_matches[refID]+=1
                else:
                    ref_matches[refID]=1
              
            prev_match = match
                      
    else:

        for line in match_fh:
            match = line.strip().split("\t")        
            # Skip header lines
            if match[0] == "@SQ" or match[0] == "@PG" or match[0] == "@HD":
                continue
            queryID = match[0]
            refID = match[2]
            if refID not in taxonomy:
                continue

            if refID in ref_matches:
                ref_matches[refID]+=1
            else:
                ref_matches[refID]=1
            
            
    all_refs = sorted(list(ref_matches.items()), key=lambda field: field[1], reverse=True)
    return all_refs


def build_conseq(sam_filename, ref_fasta):
    
    #tmp_dir = '/elo/technodrome/genomics/B17013_Orion/tmp/'

    # try:
    #     os.mkdir(tmp_dir)
    # except OSError as exc:
    #     if exc.errno != errno.EEXIST:
    #         raise
    #     pass


    print("CONSENSUS")
    print("SAM:", sam_filename) 
    print("REF FASTA:", ref_fasta)
    # Modified 14.8.19 delete = True (was False for all below)
    mpileup_temp = tempfile.NamedTemporaryFile(mode="w+",delete=False, dir=tmp_dir)
    bcftools_temp = tempfile.NamedTemporaryFile(mode="w+",delete=False, dir=tmp_dir)
    conseq_temp = tempfile.NamedTemporaryFile(mode="w+",delete=False, dir=tmp_dir)
    bam_temp = tempfile.NamedTemporaryFile(mode="w+",delete=False, dir=tmp_dir)
    sorted_bam_temp = tempfile.NamedTemporaryFile(mode="w+",delete=False, dir=tmp_dir)

    #print mpileup_temp.name, bcftools_temp.name, conseq_temp.name, bam_temp.name, sorted_bam_temp.name

    cmd = [os.path.join(install_path, "samtools-git"), "view", "-bS", sam_filename]
    print("Running command: ", " ".join(cmd))
    proc1 = subprocess.Popen(cmd, stdout=bam_temp)
    proc1.wait()
    bam_temp.seek(0)
    bam_temp.flush()
    
    cmd = [os.path.join(install_path, "samtools-git"), "sort", bam_temp.name, "-f", sorted_bam_temp.name]
    print("Running command: ", " ".join(cmd))
    proc2 = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    proc2.wait()

    #print "SORTED BAM temp:", sorted_bam_temp.name

    cmd = [os.path.join(install_path, "samtools-git-0.1.19"), "mpileup", "-Auf", ref_fasta, sorted_bam_temp.name]
    print("Running command: ", " ".join(cmd))
    proc3 = subprocess.Popen(cmd, stdout=mpileup_temp)
    proc3.wait()
    mpileup_temp.flush()
    mpileup_temp.seek(0)

    #print "MPILEUP temp:", mpileup_temp.name
    
    cmd = [os.path.join(install_path, "bcftools-0.1.19"), "view", "-cg", "-"]
    print("Running command: ", " ".join(cmd))
    proc4 = subprocess.Popen(cmd, stdin=mpileup_temp, stdout=bcftools_temp)
    proc4.wait()
    bcftools_temp.flush()
    bcftools_temp.seek(0)

    #print "BCFTOOLS temp:", bcftools_temp.name
    
    cmd = [os.path.join(install_path, "vcfutils-0.1.19.pl"), "vcf2fq"]
    print("Running command: ", " ".join(cmd))
    proc5 = subprocess.Popen(cmd, stdin=bcftools_temp, stdout=conseq_temp)
    proc5.wait()
    conseq_temp.flush()
    conseq_temp.seek(0)
    print("VCFutils done.")
    #print "CONSEQ temp:", conseq_temp.name

    conseq_dict = readFastqFH_dict(conseq_temp)

    #print conseq_dict

    conseq_temp.close()
    print("1")
    os.unlink(conseq_temp.name)

    bcftools_temp.close()
    os.unlink(bcftools_temp.name)
    print("2")
    mpileup_temp.close()
    os.unlink(mpileup_temp.name)
    print("3")
    sorted_bam_temp.close()
    os.unlink(sorted_bam_temp.name)
    print("4")
    bam_temp.close()
    os.unlink(bam_temp.name)
    
    return conseq_dict

    
def create_bt2_index(fasta_filename, index_path, index_name):

    cmd = [os.path.join(install_path, "bowtie2-build"), fasta_filename, index_path+"/"+index_name]
    print("Running command: ", " ".join(cmd))
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    proc.wait()



def create_conseq_from_fastq(ref_filename, ref_id, conseq_fh, R1, R2, conseq_nceil, conseq_bt2_score, aligner_threads, sensitivity):

    cmd = ["bash", os.path.join(install_path, "aligner-conseq.sh"),
           ref_filename,
           ref_id,
           conseq_fh.name,
           R1,
           R2,
           conseq_nceil,
           conseq_bt2_score,
           aligner_threads,
           "--"+sensitivity
          ]
    
    print("Running: ", cmd)
    
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    proc.wait()
    return

def bt2_align(tempfilename, ref, r1fname, r2fname, bt2_score, n_treatment, n_ceil, sensitivity, best_match_only, aligner_threads):

    params = []
    
    if r2fname is None:

        inputtype2params = {'.fastq': " -q ", '.fasta': " -f "}
        fileExtension = os.path.splitext(r1fname.lower())[1]
        f_type = inputtype2params[fileExtension]
        params.extend([f_type, r1fname])
        
    else:
        params.extend([" -q ", " -1 " , r1fname, " -2 ", r2fname])
        
    if int(n_treatment) > 0:
        params.extend([" --na ", n_treatment])
    else:
        params.extend([" --np ", str(abs(int(n_treatment)))])
        
    if not best_match_only:
        params.extend([" -a "])

        
    cmd = [os.path.join(install_path, "bowtie2"),
    " --threads ", aligner_threads,
    " --n-ceil ", " L,0,"+n_ceil+" ",
    " --"+sensitivity+" ",
    " --quiet ",
    " --no-unal ",
    " --ma 1 ",
    " --mp 0,0 ",
    " --rdg 1,1 ",
    " --rfg 0,0 ",
    " --score-min L,0,"+ bt2_score + " ",
    " -x ", ref]
    cmd.extend(params)
    cmd.extend([" -S ",tempfilename])

    print("Running align command: ", " ".join(cmd))
    
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    proc.wait()

    print("Bowtie2 stdout:", proc.stdout.read())

    return


def BT2opts(match):
    opts = {}
    for opt in match[11:]:
        entry = opt.split(":")
        opts[":".join(entry[0:2])] = entry[2]
    return opts


def filter_fasta(fasta_in_fh, ref_id, fasta_out_fh):
    all_refs = readFastaFH_dict(fasta_in_fh)
    records = [SeqRecord(Seq(all_refs[ref_id]), id=ref_id, description="")]
    SeqIO.write(records, fasta_out_fh, "fasta")
    return


def sam_process_template(all_lines, template, ref_id):
    matches = []
    sam_content = ""
    lines=[]

    prevmate = 0
    for line in all_lines:
        fields = line.strip().split("\t")        
        flags = int(fields[1].strip())
        align_ref = fields[2].strip()
        if align_ref != "*" and (ref_id is None or ref_id == align_ref):
            sam_content += line

    return sam_content


def filter_sam(sam_in_fh, ref_id, sam_out_fh):
    tag_content = ""
    match_content = ""
    lines = []
    
    template = None
    
    for line in sam_in_fh:
        fields = line.strip().split("\t")
         
        if fields[0][0] == "@":
                tag_content += line
        else:
            if template == fields[0].strip():
                lines.append(line)
            else:
                if template is not None:
                    match_content += sam_process_template(lines, template, ref_id)
                lines = [line]
                template = fields[0].strip()
            
    match_content += sam_process_template(lines, template, ref_id)

    sam_out_fh.write(tag_content + match_content)
    sam_out_fh.flush()


def fastq_ok(fastq_fh):
    lines = fastq_fh.read().split("\n")
    return len(lines) >= 4 # we have at least 4 lines
    
    
if __name__=="__main__":
    
    #option_parser, opts, args = parse_command_line_parameters(**script_info)
    parser = argparse.ArgumentParser(description='ASA')

    parser.add_argument('-i',"--input",
                        action='store',
                        dest='input',
                        required=True,
                        default=None,
                        help='FASTA/FASTQ files with sequences to be aligned against database".')

    parser.add_argument('-o',"--output",
                        action='store',
                        dest='output',
                        default=None,
                        help='File to be written containing db matches for each query sequence within filter range')

    parser.add_argument('-t', '--taxonomy',
                    action='store',
                    dest='taxonomy',
                    default=None,
                    help='A text file with seq ID to taxonomy".')

    parser.add_argument('--rounds', dest = 'rounds', action='store', default = "10",
         help='How many rounds to grow consensus [default: 10]"')

    parser.add_argument('--only_top_refs_ids', dest = 'only_top_refs_ids', action='store', default = None,
         help='Produce only top references id list and do not proceed to building consensus sequences')

    parser.add_argument('--top_refs', action='store', dest = 'top_refs',
         help='Cut-off for top refs"')

    parser.add_argument('--top_refs_file', dest = 'top_refs_file', action='store', default = "top_refs.txt",
         help='File name for top ref list"')

    parser.add_argument('--db_fasta', dest = 'db_fasta',
        help='A FASTA file containing taxonomy sequences.')

    parser.add_argument('--n_treatment', action='store', dest = 'n_treatment', default = "1",
         help='Award / penalty for n. For example --n_treatment -1 for penalty and --n_treatment 1 for award"')

    parser.add_argument('--aligner_threads', action='store', dest = 'aligner_threads', default = "16",
         help='Max threads the aligner can use. Default 16"')
 
    parser.add_argument('--best_match_only', default=False, action="store_true", dest = 'best_match_only',
         help='Align a query only to a best match. Default is to allow query to all refs within given threshold.')

    parser.add_argument('--candlist_bt2_score_min', default="0.97", dest = 'candlist_bt2_score_min',
         help='For candidate list: Score threshold to give to bowtie2 [default 0.97]')

    parser.add_argument('--conseq_bt2_score', default="0.97", dest = 'conseq_bt2_score',
         help='For conseq building: Score threshold to give to bowtie2 [default 0.97]')

    parser.add_argument('--conseq_round0_bt2_score_min', dest = 'conseq_round0_bt2_score_min', default=None,
         help='For conseq building: Score threshold to give to bowtie2 at round 0')
                                 
    parser.add_argument('--conseq_nceil', default="0.15", dest = 'conseq_nceil',
         help='For conseq building: n-ceil value [default 0.15]')

    parser.add_argument('--conseq_round0_nceil', default=None, dest = 'conseq_round0_nceil',
         help='For conseq building: n-ceil value at round 0')

    parser.add_argument('--candlist_nceil', default="0.30", dest = 'candlist_nceil',
         help='For candlist building: n-ceil value [default 0.30]')

    parser.add_argument('--conseq_bt2_sensitivity', default="sensitive-local", dest = 'conseq_bt2_sensitivity',
         help='For conseq building: bt2 sensitivity {sensitive-local,very-sensitive-local} [default sensitive-local]')

    parser.add_argument('--candlist_bt2_sensitivity', default="very-sensitive-local", dest = 'candlist_bt2_sensitivity',
         help='For candlist building: bt2 sensitivity {sensitive-local,very-sensitive-local} [default very-sensitive-local]')

    parser.add_argument('--candidate_refs', default=None, dest = 'candidate_refs',
         help='Do not determine candidate list, but set it. Ref is a comma separated list, for example 9663,3256')

    parser.add_argument('--taxonomy_from_fasta', default=False, action="store_true", dest = 'taxonomy_from_fasta',
         help='Take taxonomy from fasta when taxonomy file is not available.')

    parser.add_argument('--embed_in_n', default=False, action="store_true", dest = 'embed_in_n',
         help='Make sure ref sequence has Ns before and after the sequence.')

    opts = parser.parse_args()

    tmp_dir = "/tmp/"
    #tmp_dir = '/elo/technodrome/genomics/B17013_Orion/tmp/'

    fnames = opts.input.split(",")
    
    assert len(fnames) > 0
    assert len(fnames) <= 2

    if len(fnames) == 1:
        paired_end = False
        r1name = fnames[0]
        r2name = None
        r1N, fileExtension = os.path.splitext(r1name)
        
    else:
        paired_end = True
        r1name = fnames[0]
        r2name = fnames[1]

    if opts.candidate_refs is not None:
        
        run_refs = opts.candidate_refs.split(",")
        print("Using refs "+opts.candidate_refs) 

    else:
        clist_sam_temp = tempfile.NamedTemporaryFile(mode="w+",delete=False, dir = tmp_dir)
        dbindex_temp_path = tempfile.mkdtemp(dir = tmp_dir)
        create_bt2_index(opts.db_fasta, dbindex_temp_path, "refs")
        bt2_align(clist_sam_temp.name, os.path.join(dbindex_temp_path, "refs"), r1name, r2name, opts.candlist_bt2_score_min, opts.n_treatment, opts.candlist_nceil, opts.candlist_bt2_sensitivity, opts.best_match_only, opts.aligner_threads)
        shutil.rmtree(dbindex_temp_path)
        
        if opts.taxonomy_from_fasta:
             taxonomy = read_taxonomy_from_fasta(opts.db_fasta)
        else:
             taxonomy = read_taxonomy(opts.taxonomy)

        clist_sam_temp.seek(0)
        all_refs = aggregate(clist_sam_temp, paired_end, taxonomy, opts.candlist_bt2_score_min)
        clist_sam_temp.close()
        os.unlink(clist_sam_temp.name)
    
        top_refs_count = None
        if opts.top_refs:
            top_refs_count = int(opts.top_refs)
        
        all_refs_fh = open(opts.top_refs_file, "w")
        if top_refs_count:
            all_refs_fh.write("Using top "+str(top_refs_count)+" refs. \n\n")
        all_refs_fh.write("Match Count\tREF_ID\tLineage\n")
        for t in all_refs:
            ref, match_count = t
            lineage = taxonomy[ref]
            all_refs_fh.write(str(match_count)+"\t"+str(ref)+"\t"+lineage+"\n")
        all_refs_fh.close()
            
        top_refs = []
        for t in all_refs:
            ref, match_count = t
            top_refs.append(ref)
        if top_refs_count:
            top_refs = top_refs[0:top_refs_count]   

        if opts.only_top_refs_ids is not None:
            top_ref_ids_fh = open(opts.only_top_refs_ids, "w")
            for ref in top_refs:
                top_ref_ids_fh.write(str(ref)+"\n")
            top_ref_ids_fh.close()            
            print("Candidate list ids written to a file " + opts.only_top_refs_ids + ". Not proceeding to building consensus as requested.")
            sys.exit(0)
    
        run_refs = top_refs
    
    refs_temp = tempfile.NamedTemporaryFile(mode="w+",delete=False, dir = tmp_dir)
    refs_dict = readFasta_dict(opts.db_fasta)
    records = []
    for ref_id in run_refs:
        records.append(SeqRecord(Seq(refs_dict[ref_id]), id=ref_id, description=""))
    SeqIO.write(records, refs_temp, "fasta")
    refs_temp.flush()

    sam_temp = tempfile.NamedTemporaryFile(mode="w+",delete=False, dir = tmp_dir)
    
    for j, ref_id in enumerate(run_refs):

        print("NEW REF: " + ref_id + " " +str(j)+ "/"+str(len(run_refs)))

        single_fastaref_temp = tempfile.NamedTemporaryFile(mode="w+",delete=False, dir = tmp_dir)

        refs_temp.seek(0)
        print("Filtering FASTA", refs_temp.name)
        filter_fasta(refs_temp, ref_id, single_fastaref_temp)
        print("Filtering FASTA DONE", single_fastaref_temp.name)
        single_fastaref_temp.flush()

        output_path = os.path.split(opts.output)
        output_fh = open(os.path.join(output_path[0], "ref_"+ref_id+"_"+output_path[1]), "w")
    
        prevseq = None
        for i in range(int(opts.rounds)): 
    
            print("ROUND ", i, ", REF", ref_id)
            index_temp_path = tempfile.mkdtemp(dir = tmp_dir)
            create_bt2_index(single_fastaref_temp.name, index_temp_path, "refs")
    
            conseq_bt2_score = opts.conseq_bt2_score
            conseq_nceil = opts.conseq_nceil
            if i == 0:
                if opts.conseq_round0_bt2_score_min is not None: 
                    conseq_bt2_score = opts.conseq_round0_bt2_score_min
                if opts.conseq_round0_nceil is not None:
                    conseq_nceil = opts.conseq_round0_nceil
    
            bt2_align(sam_temp.name, index_temp_path+"/refs", r1name, r2name, conseq_bt2_score, opts.n_treatment, conseq_nceil, opts.conseq_bt2_sensitivity, False, opts.aligner_threads)
            
            print("Created a SAM",sam_temp.name) 
            shutil.rmtree(index_temp_path)

            try:
                fseqs = build_conseq(sam_temp.name, single_fastaref_temp.name)
                print("build consensus done")
                if fseqs[ref_id] == prevseq:
                    print("Consensus not changing, discontinuing consensus loop. Ref:", ref_id)
                    break
                prevseq = fseqs[ref_id]
                print("after prevseq")
                assert len(list(fseqs.keys())) == 1
                single_fastaref_temp.seek(0)
                single_fastaref_temp.truncate()



                fseq = fseqs[ref_id]
                print("after fseq") 
                if opts.embed_in_n:
                    print("within if statement adding N")
                    if not fseq.upper().startswith("N"*100):
                        fseq = "N"*100 + fseq
                        print("added N to START")
                    # if not fseq.upper().endswith("N"*100):
                    #     fseq = fseq + "N"*100
                    #     print "added N to END"
                #print "Passed N if statement"
                conseq_ref = [SeqRecord(Seq(fseq), id=ref_id, description="")]
                #print "after conseq_ref"
                writeFasta(single_fastaref_temp, conseq_ref)
                #print "after writefasta"
                single_fastaref_temp.flush()
                conseq_records = [SeqRecord(Seq(fseqs[ref_id]), id=ref_id+"_round_"+str(i), description="")]
                #print "after conseq_records"
                SeqIO.write(conseq_records, output_fh, "fasta")
                
            except Exception as e:
                print(e)
                print("Dropping consensus sequence "+ref_id)
                break

        output_fh.close()
                
        single_fastaref_temp.close()
        os.unlink(single_fastaref_temp.name)


    sam_temp.close()
    os.unlink(sam_temp.name)

    refs_temp.close()
    os.unlink(refs_temp.name)
