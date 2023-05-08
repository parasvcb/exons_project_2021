import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import Counter
import re,os

if len(sys.argv)!=3:
    print ("Please netre correct args inpscaffold, out")
    sys.exit()
prog, scaffoldF, out = sys.argv

total = 0
flag=0
if scaffoldF.split('.')[-1]=='gz':
    os.system('cp %s testEnsSeq.gz'%(scaffoldF))
    os.system('gzip -d -f testEnsSeq.gz')
    scaffoldF='testEnsSeq'
    flag=1

accepted_sequences=[]
# >ENSG00000185864.18|ENST00000165086.8||ENSE00002315477;ENSE00003545592;ENSE00002287434;ENSE00002258591;ENSE00003659406|21843386;21839105;21837474;21839940;21839250|21843412;21839140;21837744;21840008;21839310|||1;4;5;2;3||-1
# Sequence unavailable
for record in SeqIO.parse(scaffoldF, "fasta"):
    # print (record.seq)
    # print (record.seq.find('Sequence unavailable'))
    if record.seq.find('Sequenceunavailable')<0:
        total+=1
        newId = record.id.split('|')[2]        
        record.seq = record.seq.strip('*')
        record.id = newId
        accepted_sequences +=[record] 
        record.name = ""
        record.description = ""

if flag:
    os.remove(scaffoldF)
print (total)
SeqIO.write(accepted_sequences, out, "fasta")

#### Writing blcok below top get rid of asterisk in betwee the seuqnces (yes ensembl adds them to caution about stop gain variants)
datnew=''
with open (out) as fin:
    dat=fin.read().split('\n')
    for i in dat:
        if len(i)>0 and i[0]!='>':
            nline=re.sub(r'\*','X',i)
            datnew+='%s\n'%nline
        else:
            datnew+='%s\n'%i
with open (out,'w') as fin:
    fin.write('%s'%datnew)




####



