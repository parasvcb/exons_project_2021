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


def ensemble_write_seq(source, out_dir):
    print "In a ensemble_write_seq Funct()"
    print source, ": is a source File"
    print out_dir, ": is a output Directory"
    temp = "/home/paras/mysql/temp/"
    if re.search(r"\.gz$", source):
        # file is gzipped
        # unzip it
        # print "in"
        with gzip.open(source, "rb") as fin:
            dat = fin.read()
        with open(temp+"ens_data.tmp", "w") as fin:
            fin.write(dat)
        source = temp+"ens_data.tmp"
    if len(os.listdir(out_dir)) < 10:
        print "Writing individual fasta files from the source ensemble into %s" % out_dir
        with open(source) as fin:
            temp_dat = fin.read()
        bar = Bar('Processing', max=temp_dat.count(">"))
        del temp_dat
        for record in SeqIO.parse(source, "fasta"):
            header = record.id
            record.seq = record.seq.strip("*")
            # print header, "hh"\
            # print bool(re.search(r"ENSP\w+", header))
            protein_id = re.search(r"ENSP\w+", header).group() if re.search(r"ENSP\w+", header) else record.id
            with open(out_dir+"%s" % (protein_id), "w") as output_handle:
                SeqIO.write(record, output_handle, "fasta")
            bar.next()
        bar.finish()
        if os.path.isfile(temp+"ens_data.tmp"):
            os.remove(temp+"ens_data.tmp")
    else:
        print "sequences are already Prased"

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
