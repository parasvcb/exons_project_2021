import cPickle as pickle
with open ("../outdir/latest_5_org_2022_27_nov/objectsToDB/objectsave_10090__noStructure.pick") as fin:
	dat=pickle.load(fin)
# with open ("../outdir/latest_5_org_2022_27_nov/objectsToDB/objectsave_9606__noStructure.pick") as fin:
# 	dat=pickle.load(fin)

def UTRVar(ob):
	resdat=[]
	for gene in ob:
		has={}
		for trans in ob[gene].transcripts:
			seq=''
			for exons in trans.exons:
				seq+=exons.seq
			if seq not in has:
				has[seq]=[trans.ID]
			else:
				has[seq]+=[trans.ID]
		for i in has:
			if len(has[i])>4:
				npcount=",".join(has[i]).count("NP")
				resdat+=[(len(has[i]),str(gene), ob[gene].detail, npcount,len(has[i])-npcount, has[i])]
	resdat.sort(reverse=True)
	print ("totalTransListed\tGeneId\tGeneName\tNPListed\tOthersListed\tTransIDS")
	for i in resdat:
		print '%s\t%s\t%s\t%s\t%s\t%s'%(i[0],i[1],i[2],i[3],i[4],i[5])
UTRVar(dat)
