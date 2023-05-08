import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import Counter
import re,os

if len(sys.argv)!=4:
    print ("Please netre correct args outputMmseqs inpscaffold outdir")
    sys.exit()
prog, mmseqs, scaffoldF, out = sys.argv

total = 0
flag=0

mmseqhas={}
with open (mmseqs) as fin:
    dat=[i for i in fin.read().split('\n') if len(i)>0]
for i in dat:
    # print (i)
    key, value = i.split('\t')
    if key not in mmseqhas:
        mmseqhas[key]=[]
    mmseqhas[key]+=[value]
print ('Total were %s NR were %s'%(len(dat),len(mmseqhas)) )

# NP_000005.3     NP_000005.3
# NP_000005.3     NP_001334352.2
# NP_000005.3     ENSP00000323929.8
# NP_000006.2     NP_000006.2
# NP_000006.2     XP_016868427.1

def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]
    # use xrnage if python2

list5000=list(chunks(list(mmseqhas.keys()),5000))

if scaffoldF.split('.')[-1]=='gz':
    os.system('cp %s testEnsSeq.gz'%(scaffoldF))
    os.system('gzip -d -f testEnsSeq.gz')
    scaffoldF='testEnsSeq'
    flag=1

accepted_sequences=[]
# >ENSG00000185864.18|ENST00000165086.8||ENSE00002315477;ENSE00003545592;ENSE00002287434;ENSE00002258591;ENSE00003659406|21843386;21839105;21837474;21839940;21839250|21843412;21839140;21837744;21840008;21839310|||1;4;5;2;3||-1
# Sequence unavailable
record_dict = SeqIO.index(scaffoldF, "fasta")

for ind,val in enumerate(list5000):
    accepted_seq=[]
    cluster=ind
    for header in val:
        accepted_seq+=[record_dict[header]]
        # break
    dire = os.path.join(out,'clsuter_%s'%(cluster))
    if not os.path.isdir(dire):
        os.makedirs(dire)
    # print (accepted_seq)
    SeqIO.write(accepted_seq, os.path.join(out,'clsuter_%s.faa'%(cluster)), "fasta")
    # print (dire)
    for i in accepted_seq:
        SeqIO.write(i, os.path.join(dire,'%s'%i.id), "fasta")
    # break
   
if flag:
    os.remove(scaffoldF)