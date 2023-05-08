import constructing_data.multifasta_to_fasta as parse
import sys,os
if len (sys.argv)!=3:
    print ('Please enter correct args 1. ensemblMultifasta (biomart paras version, it has prot id at pos 3, hardcoded), 2. where to store the single sequences')
    sys.exit()
prog,inp,out=sys.argv


parse.ensemble_write_seq(inp, out, splitHeader=2)