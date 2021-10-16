args = commandArgs(trailingOnly=TRUE)
if (length(args)<2) {
  stop("Please supply two arguements 1. inpPanelsString and 2.outputFile", call.=FALSE)
} else if (length(args)==2) {
  # default output file
  systemFile=args[1]
  out=args[2]
}

library(ggplot2)

#df4<-read.csv("results/fig1/Alternate_exons_WEF.csv")
df4<-read.csv(systemFile)
#Inclusion_Range,Type,Total,Frequency
# 0-0.2,AltOnly,2518,0.28
# 0.201-0.4,AltOnly,2101,0.234
df4$Inclusion_Range <- factor(df4$Inclusion_Range, levels = unique(df4$Inclusion_Range))


hospital_names <- c(
                    `AltOnly` = paste("AltOnly ",sum(df4[df4$Type=="AltOnly",]$Total)),
                    `AltWithSS` = paste("AltWithSS ",sum(df4[df4$Type=="AltWithSS",]$Total)))
                    

gg = ggplot(data=df4, aes(x=Inclusion_Range, y=Frequency, fill=Inclusion_Range)) + geom_bar(stat="identity") 
gg <- gg + facet_wrap(~Type,nrow=1,labeller = as_labeller(hospital_names))
gg <- gg + guides(fill=FALSE)+  geom_text(aes(label=Frequency), vjust=-0.3, color="black", size=3.5) 
gg <- gg + theme (axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) 
gg <- gg + labs(y="Alt_exon_freq", x="Inclusion_freq", subtitle = "Inclusion freq of alt_exons") 
ggsave(filename = file.path(out,'Fig1_panelD_withoutExtras.pdf'))
