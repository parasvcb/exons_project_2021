library(ggplot2)

df5<-read.csv("results/fig1/panelE_Length_distribution_exonsRAW.csv")
df5$Category <- factor(df5$Category, levels = unique(df5$Category))


print (max(df5$Length))
gg = ggplot(df5, aes(x=Category,y=Length, fill=Category)) + guides(fill=FALSE) 
gg <- gg + geom_point(aes(fill=Category,color=Category), size=0.5, shape=21,
                position=position_jitter(width=0.2, height=0.1))
gg <- gg + geom_boxplot(lwd=0.1,alpha=0.5,notch=TRUE, width =0.3) + stat_summary(fun=mean, geom="point", shape=23, size=1)
gg <- gg + scale_y_continuous(breaks = c(seq(min(df5$Length), 100, by = 10),seq(100,max(df5$Length),by=500)),trans='log2')
gg <- gg + theme (axis.text.x = element_text(angle=45,hjust = 1, size = 8), axis.text.y = element_text(hjust = 1, size = 5))
gg <- gg + labs(y="Length_of_exons",x="Category",subtitle="Length_of_constitutive_and_alternate_exons")
ggsave(gg,filename = 'results/Fig1_panelE_Length_distribution.pdf')
