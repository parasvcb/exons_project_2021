from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.SubsMat import MatrixInfo as matlist
import sys
matrix = matlist.blosum62
seqA= "GIMDIEAYLERIGYKKSRNKLDLETLTDILQHQIRAVPFENLNIHCGDAMDLGLEAIFDQVVRRNRGGWCLQVNHLLYWALTTIGFETTMLGGYVYSTPAKKYSTGMIHLLLQVTIDGRNYIVDAGFGRSYQMWQPLELISGKDQPQVPCVFRLTEENGFWYLDQIRREQYIPNEEFLHSDLLEDSKYRKIYSFTLKPRTIEDFESMNTYLQTSPSSVFTSKSFCSLQTPDGVHCLVGFTLTHRRFNYKDNTDLIEFKTLSEEEIEKVLKNIFNISLQRKLVPKHGDRFFTI"
seqB= "MDIEAYLERIGYKKSRNKLDLETLTDILQHQIRAVPFENLNIHCGDAMDLGLEAIFDQVVRRNRGGWCLQVNHLLYWALTTIGFETTMLGGYVYSTPAKKYSTGMIHLLLQVTIDGRNYIVDAGFGRSYQMWQPLELISGKDQPQVPCVFRLTEENGFWYLDQIRREQYIPNEEFLHSDLLEDSKYRKIYSFTLKPRTIEDFESMNTYLQTSPSSVFTSKSFCSLQTPDGVHCLVGFTLTHRRFNYKDNTDLIEFKTLSEEEIEKVLKNIFNISLQRKLVPKHGDRFFTI"

for a in pairwise2.align.globaldx(seqA, seqB, matrix):
    a1= format_alignment(a, full_sequences=True)
    print(a1)
    sys.exit()


def stats (aln):
    # identity, cov1, cov2
    seq1 = aln[0]
    seq2 = aln[1]
    match = 0
    nongaps = 0 
    for ind, val in enumerate (seq1):
        if seq2[ind] == val:
            match +=1
        if seq2[ind]!='-' and val!='-':
            nongaps +=1
    
    seq1Length = len(seq1) - seq1.count('-')
    seq2Length = len(seq2) - seq2.count('-')

    identity1 = round(float(match)/seq1Length,3)
    identity2 = round(float(match)/seq2Length,3)
    
    cov1 = round(float(nongaps)/seq1Length,3)
    cov2 = round(float(nongaps)/seq2Length,3)
    
    return "sln1:%s, sln2: %s, idn1: %s, idn2: %s, cov1:%s, cov2:%s" %(seq1Length, seq2Length, identity1, identity2, cov1, cov2)

a = pairwise2.align.globaldx(seqA, seqB, matrix)

print (len (a))
print (a[0])

print (stats(a[0]))