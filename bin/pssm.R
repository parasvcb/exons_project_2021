library(ggplot2)

df1 = read.csv("results/fig2/Individual_alljunct_window3.csv")
df1$Type <- "1_alljunct"
df2 = read.csv("results/fig2/Individual_background_window3.csv")
df2$Type <- "2_background"
df3 = read.csv("results/fig2/Individual_AltAlt_window3.csv")
df3$Type <- "3_AltAlt"
df4 = read.csv("results/fig2/Individual_ConsAlt_window3.csv")
df4$Type <- "4_ConsAlt"
df5 = read.csv("results/fig2/Individual_ConsCons_window3.csv")
df5$Type <- "5_ConsCons"
df6 = read.csv("results/fig2/Individual_ConsCons_window3.csv")
df6$Type <- "6_ConsCons"

df7=rbind(df1,df2,df3,df4,df5,df6)

df7$SSType <- factor(df7$SSType, levels = unique(df7$SSType))

p1<-ggplot(data=df7, aes(x=SSType, y=Freq, fill=SSType)) + geom_bar(stat="identity") +
guides(fill=FALSE)+  geom_text(aes(label=Freq), vjust=-0.3, color="black", size=3.5) + facet_wrap(~Type,nrow=3) +
labs(y="Frequency", x="SSType", subtitle="Junction_Bw_Const_const_exons") + theme (axis.text.x = element_text(angle = 45, hjust = 1, size = 10))
ggsave(p1,filename="results/fig2_ssjunctions_PI.pdf")