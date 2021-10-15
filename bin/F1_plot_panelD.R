library(ggplot2)

df4<-read.csv("results/fig1/Alternate_exons_WEF.csv")
df4$Inclusion_Range <- factor(df4$Inclusion_Range, levels = unique(df4$Inclusion_Range))

gg = ggplot(data=df4, aes(x=Inclusion_Range, y=Frequency, fill=Inclusion_Range)) + geom_bar(stat="identity") 
gg <- gg + facet_wrap(~Type,nrow=1)
gg <- gg + guides(fill=FALSE)+  geom_text(aes(label=Frequency), vjust=-0.3, color="black", size=3.5) 
gg <- gg + theme (axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) 
gg <- gg + labs(y="Alt_exon_freq", x="Inclusion_freq", subtitle = "Inclusion freq of alt_exons") 
ggsave(filename = 'results/Fig1_panelD_withoutExtras.pdf')
