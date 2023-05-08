from Bio import SeqIO
import gzip
import re
import os
from progress.bar import Bar
a = '''
----> Into module of multifasta_to_fasta
'''


def description():
    '''
    This has follwoing 3 functions
        1. ensemble_write_seq(source, out_dir,txid)
        Source: source file, can be .gz or txt
        outdir: where results shoudl be stored, it sghould have txid appended in the nd, esle new directyory will be created with that appended
        txid: human or other, eg 9606
        Ensemble file if retrieved from biomart, can have multiple paareters in header line
        2. swissprot_write_seq(source, out_dir, txid):
        3. def refseqP_write_seq(source, out_dir, txid):
        Highlight: uses Biopython module to split sequneces
    '''


print a


def ensemble_write_seq(source, out_dir, splitHeader=False):
    print "In a ensemble_write_seq Funct()"
    print source, ": is a source File"
    print out_dir, ": is a output Directory"
    empty=0
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    temp = os.getcwd()
    source1 = os.path.join(temp,"ens_data.tmp")
    if re.search(r"\.gz$", source):
        print ("yes")
        with gzip.open(source, "rb") as fin:
            dat = fin.read()
        with open(source1, "w") as fin:
            fin.write(dat)
        
    if len(os.listdir(out_dir)) < 10:
        print "Writing individual fasta files from the source ensemble into %s" % out_dir
        with open(source1) as fin:
            temp_dat = fin.read()
        bar = Bar('Processing', max=temp_dat.count(">"))
        del temp_dat
        for record in SeqIO.parse(source1, "fasta"):
            header = record.id
            record.seq = record.seq.strip("*")
            # print header, "hh"\
            # print bool(re.search(r"ENSP\w+", header))
            protein_id = re.search(r"ENSP\w+", header).group() if re.search(r"ENSP\w+", header) else record.id
            if splitHeader:
                # print (header, splitHeader)
                protein_id = header.split('|')[splitHeader]
            if protein_id:
                #record.id = ">"+protein_id if record.id[0]!='>' else ""
                with open(os.path.join(out_dir,"%s" % (protein_id)), "w") as output_handle:
                    SeqIO.write(record, output_handle, "fasta")
                bar.next()
            else: 
                empty+=1
        bar.finish()
        if os.path.isfile(source1):
            os.remove(source1)
    else:
        print "sequences are already Prased"
    print ("empty sequences were %s"%empty)
    """
    if confdition because of following 
        >ENSG00000008735.14|ENST00000008876.7||ENSE00003773674;ENSE00003772801;ENSE00003768317;ENSE00003765641;ENSE00003769486;ENSE00003763622;ENSE00003731955;ENSE00003773228;ENSE00003772161;ENSE00003608148|50605368;50603841;50603626;50605562;50606921;50605825;50610707;50606658;50610212;50603133|50605443;50605064;50603719;50605734;50606991;50605934;50613981;50606765;50610310;50603498|||4;3;2;5;8;6;10;7;9;1||1
    Sequence unavailable

    """
# ensemble_write_seq()


def swissprot_write_seq(source, out_dir):
    print "In a SP_write_seq Funct()"
    print source, ": is a source File"
    print out_dir, ": is a output Directory"
    temp = "/home/paras/mysql/temp/"
    if re.search(r"\.gz$", source):
        # file is gzipped
        # unzip it
        # print "in"
        with gzip.open(source, "rb") as fin:
            dat = fin.read()
        with open(temp+"swiss_data.tmp", "w") as fin:
            fin.write(dat)
        source = temp+"swiss_data.tmp"
    if len(os.listdir(out_dir)) < 10:
        print "Writing individual fasta files from the source SP into %s" % out_dir
        with open(source) as fin:
            temp_dat = fin.read()
        bar = Bar('Processing', max=temp_dat.count(">"))
        del temp_dat
        for record in SeqIO.parse(source, "fasta"):
            header = record.id
            record.seq = record.seq.strip("*")
            protein_id = header.split("|")[1] if header.split("|") and len(header.split("|")) > 1 else header
            with open(out_dir+"%s" % (protein_id), "w") as output_handle:
                SeqIO.write(record, output_handle, "fasta")
            bar.next()
        bar.finish()
        if os.path.isfile(temp+"swiss_data.tmp"):
            os.remove(temp+"swiss_data.tmp")
    else:
        print "SP sequences are already Prased"


def refseqP_write_seq(source, out_dir):
    print "In a RefseqP_write_seq Funct()"
    print source, ": is a source File"
    print out_dir, ": is a output Directory"
    temp = "/home/paras/mysql/temp/"
    if re.search(r"\.gz$", source):
        # file is gzipped
        # unzip it
        # print "in"
        with gzip.open(source, "rb") as fin:
            dat = fin.read()
        with open(temp+"ref_data.tmp", "w") as fin:
            fin.write(dat)
        source = temp+"ref_data.tmp"
    if len(os.listdir(out_dir)) < 10:
        print "Writing individual fasta files from the source Refseq_file into %s" % out_dir
        with open(source) as fin:
            temp_dat = fin.read()
        bar = Bar('Processing', max=temp_dat.count(">"))
        del temp_dat
        for record in SeqIO.parse(source, "fasta"):
            header = record.id
            record.seq = record.seq.strip("*")
            with open(out_dir+"%s" % (header), "w") as output_handle:
                SeqIO.write(record, output_handle, "fasta")
            bar.next()
        bar.finish()
        if os.path.isfile(temp+"ref_data.tmp"):
            os.remove(temp+"ref_data.tmp")
    else:
        print "Refseq sequences are already Parsed"
