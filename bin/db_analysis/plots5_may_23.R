args = commandArgs(trailingOnly=TRUE)
if (length(args)<2) {
  stop("Please supply two arguements 1. inpDir and 2.outDir", call.=FALSE)
} else if (length(args)==2) {
  # default output file
  systemDir=args[1]
  outDir=args[2]
}

#red, #blue, #green, #black, #purple

#col_pallette = c("#d1b4cc", "#0ce861", "#f5646b", "#f48b36", "#082b95")
col_pallette = c("pink", "cadetblue", "red",  "green", "brown")
col_pallette2 = c("white","pink","white",  "cadetblue","white",  "red", "white", "green", "white","brown")

text_pallette = c("magenta","blue", "darkred",  "darkgreen", "black")


library(ggplot2)
library(ggExtra)
library(stringr)

fname1 = file.path(systemDir, "stringA0_transData_supp.tsv")
# Category	Organism	Count	ValueMean	ValueStd
# 0_Total_ISF	3_Mouse	11494	6.528	6.912
df1 <- read.csv(fname1,sep='\t')
gg  = ggplot(data=df1, aes(x=Category, y=ValueMean, fill=Organism)) + scale_fill_manual(values=col_pallette)
gg <- gg + geom_bar(stat="identity", position=position_dodge())
gg <- gg + geom_text(aes(label=round(ValueMean,2), color=Organism), vjust=-0.3, size=2.5, position=position_dodge(width = .9)) + scale_color_manual(values=text_pallette)
gg <- gg + scale_y_continuous(breaks = seq(min(df1$ValueMean), max(df1$ValueMean), by = 1))
gg <- gg + geom_errorbar(aes(ymin=ValueMean, ymax=ValueMean+ValueStd), color = "#6e8b3d", width=.25, size=0.2, position=position_dodge(.9))
gg <- gg + theme_light() + theme (axis.text.x = element_text(angle = 45, hjust = 1, size = 6), legend.position = "bottom") + ggtitle("transData_supp")
# gg <- gg + labs(y="Alt_exon_freq", x="Inclusion_freq", subtitle = paste0("ChangeInExonFrationCount, Total=", sumVal)) 
ggsave(filename =file.path(outDir,"stringA0_transData_supp.pdf"), height = 5, width =7)


fname2 = file.path(systemDir, "stringA1_GT_exonData_supp.tsv")
# Category	Organism	Count	Value
# 001	3_Mouse	3717	0.169
df2 <-read.csv(fname2,sep='\t')
gg  = ggplot(data=df2, aes(x=Category, y=Value, fill=Organism)) + scale_fill_manual(values=col_pallette)
gg  <- gg + geom_bar(stat="identity", position=position_dodge())
gg <- gg + geom_text(aes(label=round(Value,2), color=Organism), vjust=-0.3, size=2, position=position_dodge(width = .9)) + scale_color_manual(values=text_pallette)
gg <- gg + theme_light() + theme (axis.text.x = element_text(angle = 45, hjust = 1, size = 10), legend.position = "bottom") + ggtitle("GT_exonData_supp")
# gg <- gg + labs(y="Alt_exon_freq", x="Inclusion_freq", subtitle = paste0("ChangeInExonFrationCount, Total=", sumVal)) 
ggsave(filename =file.path(outDir,"stringA1_GT_exonData_supp.pdf"), height = 6, width =10)


fname3 = file.path(systemDir, "stringA2_GT_transData_supp.tsv")
# Category	Organism	Value
# 2	Mouse	0.409
df3<-read.csv(fname3,sep='\t')
gg = ggplot(data=df3, aes(x=Category, y=Value, fill=Organism)) + scale_fill_manual(values=col_pallette)
gg <- gg + geom_bar(stat="identity", position=position_dodge())
gg <- gg + geom_text(aes(label=round(Value,2), color=Organism), vjust=-0.3, size=2, position=position_dodge(width = .9)) + scale_color_manual(values=text_pallette)
gg <- gg +theme_light() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10), legend.position = "bottom") + ggtitle("GT_transData_supp")
# gg <- gg + labs(y="Alt_exon_freq", x="Inclusion_freq", subtitle = paste0("ChangeInExonFrationCount, Total=", sumVal)) 
ggsave(filename =file.path(outDir,"stringA2_GT_transData_supp.pdf"), height = 6, width =10)

fname31 = file.path(systemDir, "stringA3_GT_exonData.tsv")
# Category	Tag	Organism	Count	Value
# 001	Total	3_Mouse	2583	0.118
df31<-read.csv(fname31,sep='\t', stringsAsFactors=FALSE)
df31 <- df31[df31['Tag'] != "Total", ]
# print (head(df31))
# df31$Tag[df31$Tag == 'PseudoUTR'] <- '0PseudoUTR'
df31[df31 == 'PseudoUTR'] <- '0PseudoUTR'
df31['filCode'] <- paste(df31$Organism,df31$Tag)
print ('100')
# print (head(df31))
gg = ggplot(data=df31, aes(x=Organism, y=Value, fill=filCode)) + scale_fill_manual(values=col_pallette2)
gg <- gg + geom_bar(stat="identity", color='black', size=0.25)
gg <- gg + facet_wrap(~Category, nrow=1)
gg <- gg + geom_text(aes(label=round(Value,2), color=Organism), vjust=1, size=1.5) + scale_color_manual(values=text_pallette)
gg <- gg +theme_light() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7), legend.position = "bottom") + ggtitle("GT_exonData")
# gg <- gg + labs(y="Alt_exon_freq", x="Inclusion_freq", subtitle = paste0("ChangeInExonFrationCount, Total=", sumVal)) 
ggsave(filename =file.path(outDir,"stringA3_GT_exonData.pdf"), height = 6, width =12)

fname32 = file.path(systemDir, "stringA4_GT_transData.tsv")
# Category	Tag	Organism	Count	Value
# 001	Total	3_Mouse	2583	0.118
df32<-read.csv(fname32,sep='\t', stringsAsFactors=FALSE)
df32 <- df32[df32['Tag'] != "Total", ]
# print (head(df32))
# df32$Tag[df32$Tag == 'PseudoUTR'] <- '0PseudoUTR'
df32[df32 == 'PseudoUTR'] <- '0PseudoUTR'
df32['filCode'] <- paste(df32$Organism,df32$Tag)
print ('100')
# print (head(df32))
gg = ggplot(data=df32, aes(x=Organism, y=Value, fill=filCode)) + scale_fill_manual(values=col_pallette2)
gg <- gg + geom_bar(stat="identity", color='black', size=0.25)
gg <- gg + facet_wrap(~Category, nrow=1)
gg <- gg + geom_text(aes(label=round(Value,2), color=Organism), vjust=1, size=1.5) + scale_color_manual(values=text_pallette)
gg <- gg +theme_light() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7), legend.position = "bottom") + ggtitle("GT_transData")
# gg <- gg + labs(y="Alt_exon_freq", x="Inclusion_freq", subtitle = paste0("ChangeInExonFrationCount, Total=", sumVal)) 
ggsave(filename =file.path(outDir,"stringA4_GT_transData.pdf"), height = 6, width =12)

fname33 = file.path(systemDir, "stringA5_GT_ContourExontransDataUniqueProteinCDS.tsv")
df33<-read.csv(fname33,sep='\t')
# print (sapply(df33, class))
# df33["CodingExonCount"] <- sapply(df33["CodingExonCount"],as.numeric)
# df33["NRIsfCount"] <- sapply(df33["NRIsfCount"],as.numeric)
# print (sapply(df33, class))
# print (head(df33))
# print (unique(df33$CodingExonCount))
# print (unique(df33$NRIsfCount))
# print (unique(df33$Organism))
# Gene	CodingExonCount	NRIsfCount	Organism
# 21846	23	1	3_Mouse
# 546143	3	3	3_Mouse

print (max(df33$'CodingExonCount'))
print (max(df33$'NRIsfCount'))
print (min(df33$'CodingExonCount'))
print (min(df33$'NRIsfCount'))
gg <- ggplot(data=df33, aes(x=NRIsfCount,y=CodingExonCount, color=Organism))
# scale_colour_manual(values=c("red","lightblue","green","yellow","black"))
# gg <- gg + geom_point(aes(color=Organism, shape=Organism), alpha=0.7, size=2)  + scale_colour_manual(values=col_pallette)
gg <- gg + geom_point(aes(color=Organism), alpha=0.3, size=2) + scale_colour_manual(values=col_pallette)
gg <- gg + scale_y_continuous(breaks = seq(0, 200, by = 10))
gg <- gg + scale_x_continuous(breaks = seq(0, 200, by = 10))
gg <- gg + geom_vline(xintercept = 2, color = "black", size=0.5)
gg <- gg + geom_vline(xintercept = 4, linetype="dashed", color = "blue", size=0.5)
gg <- gg + geom_hline(yintercept = 2, color = "black", size=0.5)
gg <- gg + geom_hline(yintercept = 4, linetype="dashed", color = "blue", size=0.5)
gg <- gg + theme (axis.text.x = element_text( hjust = 1, size = 7, angle = 45), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20), axis.text.y = element_text( hjust = 1, size = 7), legend.position="bottom", panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey95",size = 0.25))
ag <- ggMarginal(gg,type = "histogram", groupColour = TRUE, groupFill = TRUE)
pdf(file.path(outDir,"stringA5_GT_ContourExontransDataUniqueProteinCDS_noLog_1.pdf"), width=7, height=7)
print (ag)
dev.off()

gg <- gg + coord_trans(x="log2", y="log2")
ag1 <- ggMarginal(gg,type = "histogram", groupColour = TRUE, groupFill = TRUE)
pdf(file.path(outDir,"stringA5_GT_ContourExontransDataUniqueProteinCDS_log_1.pdf"), width=7, height=7)
print (ag1)
dev.off()
# gg <- gg + geom_point(aes(shape =Organism), alpha=0.7, size=2)
#  coord_trans(x="log2", y="log2")
# gg <- gg + geom_density_2d_filled(contour_var="count", alpha=0.5)
# gg <- gg + geom_density_2d_filled(aes(color=Organism))  + scale_fill_manual(values = col_pallette) + scale_colour_manual(values=col_pallette)
# gg <- gg + stat_density_2d(geom = "polygon", aes(x=NRIsfCount,y=CodingExonCount,alpha = ..level.., fill = Organism, colour=Organism))  + scale_fill_manual(values = col_pallette) + scale_colour_manual(values=col_pallette)
# gg <- gg + geom_point(data=df33, aes(x=NRIsfCount,y=CodingExonCount), alpha=1, size=5)
# gg <- gg + theme (axis.text.x = element_text( hjust = 1, size = 20, angle = 45), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20), axis.text.y = element_text( hjust = 1, size = 20), legend.position="top", panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey95",size = 0.25))
# ggsave(filename =file.path(outDir,"stringA5_GT_ContourExontransDataUniqueProteinCDS.tiff"), gg, height =10, width=12)

###################################################################################################################################################################################
print ("0")
fname4 = file.path(systemDir, "stringB1_exonData_exonGeneFrac.tsv")
# Category	Facet	Organism	ValueFractionExons	ValueFractionGenes
# Total_exons	0_extras	Mouse	1.0	1.0
df4<-read.csv(fname4,sep='\t')
# gg = ggplot(data=df4, aes(x=Category, y=ValueFractionExons, fill=Organism)) + scale_fill_manual(values=col_pallette)
# gg <- gg + geom_bar(stat="identity", position=position_dodge())
# gg <- gg + scale_y_continuous(breaks = seq(min(df4$ValueFractionExons), max(df4$ValueFractionExons), by = 0.10))
# gg <- gg + geom_text(aes(label=round(ValueFractionExons,2), color=Organism), vjust=-1, size=1.3, position=position_dodge(width = .9)) + scale_color_manual(values=text_pallette)
# gg <- gg + facet_wrap(~Facet, nrow=1)
# gg <- gg +theme_light() + theme (axis.text.x = element_text(angle = 45, hjust = 1, size = 8), legend.position = "bottom") + ggtitle("exonData_exFrac")
# # gg <- gg + labs(y="Alt_exon_freq", x="Inclusion_freq", subtitle = paste0("ChangeInExonFrationCount, Total=", sumVal)) 
# ggsave(filename =file.path(outDir,"stringB1_exonData_exonFrac.jpg"), width = 10, height = 5)

print (max(df4$ValueFractionGenes))
gg = ggplot(data=df4, aes(x=Category, y=ValueFractionGenes, fill=Organism)) + scale_fill_manual(values=col_pallette)
gg <- gg + geom_col(stat="identity", position=position_dodge()) + coord_flip()
gg <- gg + geom_col(data=df4, aes(x=Category, y=-1*(ValueFractionExons), fill=Organism), stat="identity", position=position_dodge()) + coord_flip()
gg <- gg + geom_text(aes(label=round(ValueFractionGenes,2), color=Organism),  vjust=-1.1, size=1.7, position=position_dodge(width = .9)) + scale_color_manual(values=text_pallette)
gg <- gg + geom_text(aes(y=-1*ValueFractionExons, label=round(ValueFractionExons,2), vjust = 1.1, color=Organism), size=1.7, position=position_dodge(width = .9)) + scale_color_manual(values=text_pallette)
gg <- gg + scale_y_continuous(breaks = seq(min(-1*df4$ValueFractionExons), max(df4$ValueFractionGenes), by = 0.10))
gg <- gg + geom_hline(yintercept=0, size=0.25)
gg <- gg + geom_hline(yintercept=0.25, size=0.1, color="blue")
gg <- gg + geom_hline(yintercept=0.50, size=0.1, color="green")
gg <- gg + geom_hline(yintercept=0.75, size=0.1, color="red")
gg <- gg + geom_hline(yintercept=-0.25, size=0.1, color="blue")
gg <- gg + geom_hline(yintercept=-0.50, size=0.1, color="green")
gg <- gg + geom_hline(yintercept=-0.75, size=0.1, color="red")
gg <- gg + facet_wrap(~Facet, nrow=1)
gg <- gg + theme_light() + theme (axis.text.x = element_text(angle = 45, hjust = 1, size = 8), axis.text.y = element_text(angle = 45, hjust = 1, size = 8), legend.position = "right")  + ggtitle("exonData_geneFracRight_exonFracLeft")
# gg <- gg + labs(y="Alt_exon_freq", x="Inclusion_freq", subtitle = paste0("ChangeInExonFrationCount, Total=", sumVal)) 
ggsave(filename =file.path(outDir,"stringB1_exonData_geneFracRight_exonFracLeft.jpg"), width = 9, height = 6)

fname41 = file.path(systemDir, "stringB2_exonData_exonGeneFrac_additional.tsv")
# Category	Facet	Organism	ValueFractionExons	ValueFractionGenes
# Total_exons	0_extras	Mouse	1.0	1.0
df41<-read.csv(fname41,sep='\t')
gg = ggplot(data=df41, aes(x=Category, y=ValueFractionGenes, fill=Organism)) + scale_fill_manual(values=col_pallette)
gg <- gg + geom_col(stat="identity", position=position_dodge()) + coord_flip()
gg <- gg + geom_col(data=df41, aes(x=Category, y=-1*(ValueFractionExons), fill=Organism), stat="identity", position=position_dodge()) + coord_flip()
gg <- gg + geom_text(aes(label=round(ValueFractionGenes,2), color=Organism),  vjust=-1.1, size=1.3, position=position_dodge(width = .9)) + scale_color_manual(values=text_pallette)
gg <- gg + geom_text(aes(y=-1*ValueFractionExons, label=round(ValueFractionExons,2), vjust = 1.1, color=Organism), size=1.3, position=position_dodge(width = .9)) + scale_color_manual(values=text_pallette)
gg <- gg + scale_y_continuous(breaks = seq(min(-1*df41$ValueFractionExons), max(df41$ValueFractionGenes), by = 0.10))
gg <- gg + geom_hline(yintercept=0, size=0.25)
gg <- gg + geom_hline(yintercept=0.25, size=0.1, color="blue")
gg <- gg + geom_hline(yintercept=0.50, size=0.1, color="green")
gg <- gg + geom_hline(yintercept=0.75, size=0.1, color="red")
gg <- gg + geom_hline(yintercept=-0.25, size=0.1, color="blue")
# gg <- gg + facet_wrap(~Facet, nrow=1)
gg <- gg + theme_light() + theme (axis.text.x = element_text(angle = 45, hjust = 1, size = 8), axis.text.y = element_text(angle = 45, hjust = 1, size = 8), legend.position = "right")  + ggtitle("exonData_additional_geneFracRight_exonFracLeft")
# gg <- gg + labs(y="Alt_exon_freq", x="Inclusion_freq", subtitle = paste0("ChangeInExonFrationCount, Total=", sumVal)) 
ggsave(filename =file.path(outDir,"stringB2_exonData_geneFracRight_exonFracLeft_additional.jpg"), width = 8, height = 12)

###################################################################################################################################################################################
print ("2")
fname5 = file.path(systemDir, "stringB1_exonData_meanStdCounts.tsv")
# Category	Facet	Organism	ValueMean	ValueStd
# Total_exons	0_extras	Mouse	14.671	9.919
df5<-read.csv(fname5,sep='\t')
print ("2.1")
gg = ggplot(data=df5, aes(x=Category, y=ValueMean, fill=Organism)) + scale_fill_manual(values=col_pallette)
gg <- gg + geom_bar(stat="identity", position=position_dodge())
gg <- gg + geom_text(aes(label=round(ValueMean,2)), color="black", vjust=-1.1, size=1.3, position=position_dodge(width = .9)) + scale_color_manual(values=text_pallette)
gg <- gg + scale_y_continuous(breaks = seq(min(df5$ValueMean), max(df5$ValueMean), by = 2))
gg <- gg + geom_errorbar(aes(ymin=ValueMean, ymax=ValueMean+ValueStd), color = "black", size = 0.1, width=.25, position=position_dodge(.9))
gg <- gg + facet_wrap(~Facet, nrow=1)
gg <- gg + theme_light() + theme (axis.text.x = element_text(angle = 45, hjust = 1, size = 7), legend.position = "bottom") + ggtitle("exonData_meanstd")
# gg <- gg + labs(y="Alt_exon_freq", x="Inclusion_freq", subtitle = paste0("ChangeInExonFrationCount, Total=", sumVal)) 
ggsave(filename =file.path(outDir,"stringB1_exonData_meanStdCounts.jpg"), width = 9, height = 5)

fname5 = file.path(systemDir, "stringB2_exonData_meanStdCounts_additional.tsv")
# Category	Facet	Organism	ValueMean	ValueStd
# Total_exons	0_extras	Mouse	14.671	9.919
df51<-read.csv(fname5,sep='\t')
print ("2.1")
gg = ggplot(data=df51, aes(x=Category, y=ValueMean, fill=Organism)) + scale_fill_manual(values=col_pallette)
gg <- gg + geom_bar(stat="identity", position=position_dodge())
gg <- gg + geom_text(aes(label=round(ValueMean,2)), color="black", vjust=-1.1, size=1.8, position=position_dodge(width = .9)) + scale_color_manual(values=text_pallette)
gg <- gg + scale_y_continuous(breaks = seq(min(df51$ValueMean), max(df51$ValueMean), by = 0.5))
gg <- gg + geom_errorbar(aes(ymin=ValueMean, ymax=ValueMean+ValueStd), color = "black", size = 0.1, width=.25, position=position_dodge(.9))
# gg <- gg + facet_wrap(~Facet, nrow=1)
gg <- gg + theme_light() + theme (axis.text.x = element_text(angle = 45, hjust = 1, size = 7), legend.position = "bottom") + ggtitle("exonData_meanstd")
# gg <- gg + labs(y="Alt_exon_freq", x="Inclusion_freq", subtitle = paste0("ChangeInExonFrationCount, Total=", sumVal)) 
ggsave(filename =file.path(outDir,"stringB2_exonData_meanStdCounts_additional.jpg"), width = 12, height = 4)

print ("3")
fname6 = file.path(systemDir, "stringC_exonLength.tsv")
# Category	Organism	ValueMean	ValueMedian	ValueStd
# CodingStrict	Mouse	50.653	39.0	68.749
df6<-read.csv(fname6,sep='\t')
gg = ggplot(data=df6, aes(x=Category, y=ValueMean, fill=Organism)) + scale_fill_manual(values=col_pallette)
gg <- gg + geom_bar(stat="identity", position=position_dodge())
gg <- gg + geom_errorbar(aes(ymin=ValueMean, ymax=ValueMean+ValueStd), width=.25, size =0.2, position=position_dodge(.9))
gg <- gg + geom_text(aes(label=round(ValueMean,1), color=Organism), vjust=-0.3, size=2.5, position=position_dodge(width = .9)) + scale_color_manual(values=text_pallette)
gg <- gg + theme_light() + theme (axis.text.x = element_text(angle = 45, hjust = 1, size = 10), legend.position = "bottom") + ggtitle("exonData_LengthMean")
# gg <- gg + labs(y="Alt_exon_freq", x="Inclusion_freq", subtitle = paste0("ChangeInExonFrationCount, Total=", sumVal)) 
ggsave(filename =file.path(outDir,"stringC_exonLength_mean.jpg"), width = 12, height = 5)

print ("4")
gg = ggplot(data=df6, aes(x=Category, y=ValueMedian, fill=Organism)) + scale_fill_manual(values=col_pallette)
gg <- gg + geom_bar(stat="identity", position=position_dodge()) + scale_fill_manual(values=col_pallette)
gg <- gg + geom_text(aes(label=ValueMedian, color=Organism), vjust=-0.3, size=2.5, position=position_dodge(width = .9)) + scale_color_manual(values=text_pallette)
gg <- gg + geom_errorbar(aes(ymin=ValueMedian, ymax=ValueMedian+ValueStd), width=.25, size =0.2, position=position_dodge(.9))
gg <- gg + theme_light() + theme (axis.text.x = element_text(angle = 45, hjust = 1, size = 10), legend.position = "bottom") + ggtitle("exonData_LengthMedian")
# gg <- gg + labs(y="Alt_exon_freq", x="Inclusion_freq", subtitle = paste0("ChangeInExonFrationCount, Total=", sumVal)) 
ggsave(filename =file.path(outDir,"stringC_exonLength_median.jpg"), width = 10, height = 5)


fname7 = file.path(systemDir, "stringD_exonWEF.tsv")
# Category	Organism    Range	Frequency	Total
# AltOnly	Mouse	0-0.2	0.104	3000
# AltOnly	Mouse	0.201-0.4	0.138	3962
df7<-read.csv(fname7,sep='\t')
gg = ggplot(data=df7, aes(x=Category, y=Frequency, fill=Organism)) + scale_fill_manual(values=col_pallette)
gg <- gg + geom_bar(stat="identity", position=position_dodge())
gg <- gg + facet_wrap(~Range, nrow=1)
gg <- gg + geom_text(aes(label=round(Frequency,2), color=Organism), vjust=-0.3, size=2, position=position_dodge(width = .9)) + scale_color_manual(values=text_pallette)
gg <- gg + theme_light() + theme (axis.text.x = element_text(angle = 45, hjust = 1, size = 10), legend.position = "bottom") + ggtitle("exonData_WEFFreq")
# gg <- gg + labs(y="Alt_exon_freq", x="Inclusion_freq", subtitle = paste0("ChangeInExonFrationCount, Total=", sumVal)) 
ggsave(filename =file.path(outDir,"stringD_exonWEF.jpg"), width = 10, height = 5)

fname8 = file.path(systemDir, "stringE_exonWEF.tsv")
# Category        Organism        Value
# AltOnly 3_Mouse 0.334
# AltOnly 3_Mouse 0.667
# AltOnly 3_Mouse 0.13

df7[c('First_Val', 'Last_Val')] <- str_split_fixed(df7$Range, '-', 2)
# print (unique(df7$Category))

df8 <- read.csv(fname8,sep='\t')
# print (unique(df8$Category))
# print (head (df7))
# print (head (df8))

# print (unique(df7$Category))

gg = ggplot(data=df8, aes(x=Value, y=..scaled.., color=Organism, fill=Organism), size = 2) + scale_color_manual(values=col_pallette) +scale_fill_manual(values=col_pallette)
print (1)
# gg <- gg + geom_histogram(aes(y=..density..), position=position_dodge())
# gg <- gg + geom_density(data=df8, aes(x=Value, y=..scaled..), alpha=0, lwd=1)
gg <- gg + geom_density(aes(y=..scaled..), alpha=0, lwd=1)
print (2)
gg <- gg + geom_bar(data=df7, aes(x=as.double(Last_Val)-0.1, y=Frequency), alpha= 0.5, stat="identity", position=position_dodge())
print (3)
# gg <- gg + geom_text(data=df7, aes(label=round(Frequency,2), color=Organism), vjust=-0.3, size=2, position=position_dodge(width = .9)) + scale_color_manual(values=text_pallette)
# gg <- gg + geom_text(data=df7, inherit.aes = F, aes(x=as.double(Last_Val)-0.1, y=Frequency, label=round(Frequency,2), color=Organism), vjust=-0.3, size=2, position=position_dodge(width = .9)) + scale_color_manual(values=text_pallette)
print (4)
gg <- gg + scale_x_continuous(breaks = seq(0, 1, by = 0.2))
gg <- gg + scale_y_continuous(breaks = seq(0, 1, by = 0.1))
print (5)
gg <- gg + facet_wrap(~Category, nrow=2)
print (6)
gg <- gg + theme_light() + theme (axis.text.x = element_text(angle = 45, hjust = 1, size = 10), legend.position = "right") + ggtitle("exonData_WEFFreq, histograms spans range 0-0.2, ... 0.8-1.0")
# gg <- gg + labs(y="Alt_exon_freq", x="Inclusion_freq", subtitle = paste0("ChangeInExonFrationCount, Total=", sumVal)) 
ggsave(filename =file.path(outDir,"stringE_exonWEF.jpg"), width = 8, height = 6)
