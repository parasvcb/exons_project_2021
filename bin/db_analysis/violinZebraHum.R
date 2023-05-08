library(ggplot2)

df1 = read.csv("tableViolin.tsv",sep='\t')
df1$type='X'
# GeneId	humanIsf	FishIsf	diff	frac
# 28990	4	2	2	0.5
# 23440	1	4	-3	-3.0
# 55760	5	2	3	0.6
# 151556	7	4	3	0.429
# 348	5	3	2	0.4
# 962	5	4	1	0.2

print (head(df1))

p1<-ggplot(data=df1, aes(x=type,y=diff)) + geom_violin() + labs(y="isf count subtraction", x="487 genes") + geom_boxplot(width=0.1)
ggsave(p1,filename="Violin_zebra_diff.pdf")

p2<-ggplot(data=df1, aes(x=type, y=frac)) + geom_violin() + labs(y="isf count subtraction normalized by human isf", x="487 genes") + geom_boxplot(width=0.1)
ggsave(p2,filename="Violin_zebra_frac.pdf")
