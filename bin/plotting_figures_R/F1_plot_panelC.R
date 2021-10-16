args = commandArgs(trailingOnly=TRUE)
if (length(args)<2) {
  stop("Please supply two arguements 1. inpFile and 2.outFolder", call.=FALSE)
} else if (length(args)==2) {
  # default output file
  systemFile=args[1]
  out=args[2]
}

library('ggplot2')
df3<-read.csv(systemFile)
#df3<-read.csv("results/fig1/Unique_Coods.csv")
print (head (df3))

#df3$Category <- factor(df3$Category, levels = unique(df3$Category))

gg=ggplot(df3, aes(x=Category,y=Exons, fill=Category)) + guides(fill=FALSE)
gg <- gg + facet_wrap(~Type,nrow=1)
gg <- gg + geom_boxplot(notch=FALSE) + stat_summary(fun=mean, geom="point", shape=23, size=2)
gg <- gg + scale_y_continuous(breaks = seq(min(df3$Exons), max(df3$Exons), by = 2),trans='pseudo_log')
gg <- gg + theme (axis.text.x = element_text(angle = 45, hjust = 1, size = 5),axis.text.y = element_text(size = 5))
gg <- gg + ylab("Exon_Count/Gene")
gg <- gg + xlab("Categories")
#p3 <- p3 + coord_cartesian(ylim=c(0, 50))
ggsave(gg,filename = file.path(out,'Fig1_panelC.pdf'))

df3_1 = df3[!df3$Type=='Extras',]
df3_2 = df3_1[!df3_1$Exons<1,]

gg1=ggplot(data=df3_1,aes(x=Category,y=Exons, fill=Category)) + guides(fill=FALSE)
gg1 <- gg1 + facet_wrap(~Type,nrow=1)
gg1 <- gg1 + geom_point(aes(fill=Category,color=Category), size=0.5, shape=21,
                position=position_jitter(width=0.2, height=0.1))
gg1 <- gg1 + geom_boxplot(lwd=0.1,notch=FALSE,alpha=0.5) + stat_summary(fun=mean, geom="point", shape=23, size=1)
gg1 <- gg1 + scale_y_continuous(breaks = seq(min(df3$Exons), max(df3$Exons), by = 2),trans='pseudo_log')
gg1 <- gg1 + theme (axis.text.x = element_text(angle = 45, hjust = 1, size = 5),axis.text.y = element_text(size = 5))
gg1 <- gg1 + ylab("Exon_Count/Gene")
gg1 <- gg1 + xlab("Categories")
#p3 <- p3 + coord_cartesian(ylim=c(0, 50))
ggsave(gg1,filename = file.path(out,'Fig1_panelC_withoutExtras.pdf'))

gg2=ggplot(data=df3_2,aes(x=Category,y=Exons, fill=Category)) + guides(fill=FALSE)
gg2 <- gg2 + facet_wrap(~Type,nrow=1)
gg2 <- gg2 + geom_boxplot(notch=FALSE) + stat_summary(fun=mean, geom="point", shape=23, size=1)
gg2 <- gg2 + scale_y_continuous(breaks = seq(min(df3$Exons), max(df3$Exons), by = 2),trans='pseudo_log')
gg2 <- gg2 + theme (axis.text.x = element_text(angle = 45, hjust = 1, size = 5),axis.text.y = element_text(size = 5))
gg2 <- gg2 + ylab("Exon_Count/Gene")
gg2 <- gg2 + xlab("Categories")
#p3 <- p3 + coord_cartesian(ylim=c(0, 50))
ggsave(gg2,filename = file.path(out,'Fig1_panelC_withoutExtras_NonZeros.pdf'))


# gg <- ggplot(data=dfgen_primary,aes(x=polnew,y=mean,fill=hbtype))
# gg <- gg + geom_bar(position="stack", stat="identity") + geom_col(position = position_stack(reverse = TRUE))
# gg <- gg + scale_fill_manual(values=c("gray","green","red"))
# #gg <- gg + geom_text(aes(label=mean), position = position_stack(vjust = 0.5)) #vjust=1, color="black", size=1)
# gg <- gg + geom_errorbar(aes(x=polnew,y=vertmean,ymin=vertmean-std, ymax=vertmean+std,color = hbtype, width=0.3),position="identity") #+ geom_col(position = position_stack(reverse = TRUE))
# gg <- gg + scale_color_manual(values=c("brown","black","blue"))
# gg <- gg + ylab("Hb-count")
# gg <- gg + xlab("Polymer")
# gg <- gg + facet_wrap(~hbfrom,nrow=1)
# gg <- gg + theme (axis.text.x = element_text( hjust = 1, size = 7, angle = 45),axis.text.y = element_text( hjust = 1, size = 7), legend.position="top", panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey95",size = 0.25))
# ggsave(filename = file.path(out,pputfile,'_','Fig5_raw_plot_facet.pdf'))