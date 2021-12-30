args = commandArgs(trailingOnly=TRUE)
if (length(args)<2) {
  stop("Please supply two arguements 1. inpfile and 2.outFile", call.=FALSE)
} else if (length(args)==2) {
  # default output file
  systemFile=args[1]
  out=args[2]
}
df4<-read.csv(systemFile,sep='\t')
library(ggplot2)

gg = ggplot(data=df4, aes(x=FractionCovered, y=Freq)) + geom_bar(stat="identity") 
gg <- gg + facet_wrap(~Tag,nrow=2)
gg <- gg + guides(fill=FALSE)+  geom_text(aes(label=Freq), vjust=-0.3, color="black", size=3.5) 
gg <- gg + theme (axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) 
#gg <- gg + labs(y="Alt_exon_freq", x="Inclusion_freq", subtitle = "Inclusion freq of alt_exons") 
ggsave(filename = paste0(out,'ATITCore_freq.pdf'))

gg = ggplot(data=df4, aes(x=FractionCovered, y=CumFreq)) + geom_bar(stat="identity") 
gg <- gg + facet_wrap(~Tag,nrow=2)
gg <- gg + guides(fill=FALSE)+  geom_text(aes(label=CumFreq), vjust=-0.3, color="black", size=3.5) 
gg <- gg + theme (axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) 
#gg <- gg + labs(y="Alt_exon_CumFreq", x="Inclusion_CumFreq", subtitle = "Inclusion CumFreq of alt_exons") 
ggsave(filename = paste0(out,'ATITCore_CumFreq.pdf'))
