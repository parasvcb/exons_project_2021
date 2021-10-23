library(ggplot2)
args = commandArgs(trailingOnly=TRUE)
if (length(args)<2) {
  stop("Please supply two arguements 1. Dir and 2. inpPanelsString", call.=FALSE)
} else if (length(args)==2) {
  # default output file
  systemFile=args[1]
  string_panels=args[2]
  out=basename(systemFile)
}



df1 = read.csv(paste0(file.path(systemFile, "F2_Individual"),string_panels,"alljunct_window3.csv"))
df1$Type <- "1_alljunct"
df2 = read.csv(paste0(file.path(systemFile, "F2_Individual"),string_panels,"background_window3.csv"))
df2$Type <- "2_background"
df3 = read.csv(paste0(file.path(systemFile, "F2_Individual"),string_panels,"AltAlt_window3.csv"))
df3$Type <- "3_AltAlt"
df4 = read.csv(paste0(file.path(systemFile, "F2_Individual"),string_panels,"ConsAlt_window3.csv"))
df4$Type <- "4_ConsAlt"
df5 = read.csv(paste0(file.path(systemFile, "F2_Individual"),string_panels,"ConsCons_window3.csv"))
df5$Type <- "5_ConsCons"
df6 = read.csv(paste0(file.path(systemFile, "F2_Individual"),string_panels,"ConsCons_window3.csv"))
df6$Type <- "6_ConsCons"



df7=rbind(df1,df2,df3,df4,df5,df6)

count_df <- c(
                    `1_alljunct` = paste("1_alljunct",sum(df7[df7$Type=="1_alljunct",]$Total)),
                    `2_background` = paste("2_background",sum(df7[df7$Type=="2_background",]$Total)),
                    `3_AltAlt` = paste("3_AltAlt",sum(df7[df7$Type=="3_AltAlt",]$Total)),
                    `4_ConsAlt` = paste("4_ConsAlt",sum(df7[df7$Type=="4_ConsAlt",]$Total)),
                    `5_ConsCons` = paste("5_ConsCons",sum(df7[df7$Type=="5_ConsCons",]$Total)),
                    `6_ConsCons` = paste("6_ConsCons",sum(df7[df7$Type=="6_ConsCons",]$Total))
                    )
df7$SSType <- factor(df7$SSType, levels = unique(df7$SSType))

p1<-ggplot(data=df7, aes(x=SSType, y=Freq, fill=SSType)) + geom_bar(stat="identity") +
guides(fill=FALSE)+  geom_text(aes(label=Freq), vjust=-0.3, color="black", size=1.8) + facet_wrap(~Type,nrow=3,labeller=as_labeller(count_df)) +
labs(y="Frequency", x="SSType") + theme (axis.text.x = element_text(angle = 45, hjust = 1, size = 10))
ggsave(p1,filename=file.path(systemFile,"junctions_PI.pdf"))