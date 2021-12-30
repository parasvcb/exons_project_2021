args = commandArgs(trailingOnly=TRUE)
if (length(args)<2) {
  stop("Please supply two arguements 1. inpfile and 2.outputFile", call.=FALSE)
} else if (length(args)==2) {
  # default output file
  systemFile=args[1]
  out=args[2]
}

library(ggplot2)

df4<-read.csv(systemFile,sep='\t')
# Range	Total	Freq	CumTotal	CumFreq
# 0.0_0.01	6482	1.0	6482	1.0
# 0.01_0.1	0	0.0	6482	1.0
# 0.101_0.2	0	0.0	6482	1.0

df4$Range <- factor(df4$Range, levels = unique(df4$Range))


sumVal=sum(df4$Total)
gg = ggplot(data=df4, aes(x=Range, y=Freq)) + geom_bar(stat="identity") 
#gg <- gg + facet_wrap(~Type,nrow=1,labeller = as_labeller(hospital_names))
gg <- gg + guides(fill=FALSE)+  geom_text(aes(label=Freq), vjust=-0.3, color="black", size=3.5) 
gg <- gg + theme (axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) 
gg <- gg + labs(y="Alt_exon_freq", x="Inclusion_freq", subtitle = paste0("ChangeInExonFrationCount, Total=", sumVal)) 
ggsave(filename = out)
