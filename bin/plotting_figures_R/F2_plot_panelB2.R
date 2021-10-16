args = commandArgs(trailingOnly=TRUE)
if (length(args)<2) {
  stop("Please supply two arguements 1. inpPanelsString and 2.outputFile", call.=FALSE)
} else if (length(args)==2) {
  # default output file
  systemFile=args[1]
  out=args[2]
}


library('ggplot2')
df3<-read.csv(paste0("results/ExonsJunction",systemFile,"_and_majorSS_exonsJunction_conservation.csv"),sep='\t')

#df3<-read.csv("results/fig2/ExonsJunctionWEF_and_majorSS_exonsJunction_conservation.csv",sep='\t')
df3$wrap <- paste(df3$ExonType,df3$WEFRange)
print (head (df3))
print (unique(df3$wrap))
print (paste("cons_cons 0.25_0.5",sum(df3[df3$wrap=="cons_cons 0.25_0.5",]$Count)))

# ExonType	WEFRange	Junction	OccurenceRange	OccureneceTotal	OccureneceFreq
# cons_cons	0.75_1.01	CC	0.601_0.7	0	0.0
# cons_cons	0.75_1.01	CC	0_0.1	377	0.394

for (i in unique(df3$ExonType)) {
    print (i)
    dfUnique=df3[df3$ExonType==i,]
    dfUnique$wefJunction <- paste(dfUnique$WEFRange,dfUnique$Junction)
    #print (head(dfUnique))
    juncNames <- c(
                    `0.75_1.01 HH` = paste("0.75_1.01 HH",sum(dfUnique[dfUnique$wefJunction=="0.75_1.01 HH",]$OccureneceTotal)),
                    `0.75_1.01 CC` = paste("0.75_1.01 CC",sum(dfUnique[dfUnique$wefJunction=="0.75_1.01 CC",]$OccureneceTotal)),
                    `0.75_1.01 EE` = paste("0.75_1.01 EE",sum(dfUnique[dfUnique$wefJunction=="0.75_1.01 EE",]$OccureneceTotal)),
                    
                    `0_0.25 HH` = paste("0_0.25 HH",sum(dfUnique[dfUnique$wefJunction=="0_0.25 HH",]$OccureneceTotal)),
                    `0_0.25 CC` = paste("0_0.25 CC",sum(dfUnique[dfUnique$wefJunction=="0_0.25 CC",]$OccureneceTotal)),
                    `0_0.25 EE` = paste("0_0.25 EE",sum(dfUnique[dfUnique$wefJunction=="0_0.25 EE",]$OccureneceTotal)),
                    
                    `0.25_0.5 HH` = paste("0.25_0.5 HH",sum(dfUnique[dfUnique$wefJunction=="0.25_0.5 HH",]$OccureneceTotal)),
                    `0.25_0.5 CC` = paste("0.25_0.5 CC",sum(dfUnique[dfUnique$wefJunction=="0.25_0.5 CC",]$OccureneceTotal)),
                    `0.25_0.5 EE` = paste("0.25_0.5 EE",sum(dfUnique[dfUnique$wefJunction=="0.25_0.5 EE",]$OccureneceTotal)),
                    
                    `0.5_0.75 HH` = paste("0.5_0.75 HH",sum(dfUnique[dfUnique$wefJunction=="0.5_0.75 HH",]$OccureneceTotal)),
                    `0.5_0.75 CC` = paste("0.5_0.75 CC",sum(dfUnique[dfUnique$wefJunction=="0.5_0.75 CC",]$OccureneceTotal)),
                    `0.5_0.75 EE` = paste("0.5_0.75 EE",sum(dfUnique[dfUnique$wefJunction=="0.5_0.75 EE",]$OccureneceTotal))
    )
    print (unique(dfUnique$wefJunction))
    print (unique(dfUnique$WEFRange))
    #gg=ggplot(dfUnique, aes(x=OccurenceRange,y=OccureneceFreq)) + geom_bar(stat="identity") 
    gg=ggplot(dfUnique, aes(x=OccurenceRange,y=cumulativeFreq)) + geom_bar(stat="identity") 

    gg <- gg + facet_wrap(~wefJunction,nrow=4,labeller = as_labeller(juncNames))
    gg <- gg + guides(fill=FALSE)+  geom_text(aes(label=cumulativeFreq), vjust=-0.3, color="black", size=1.6) 
    gg <- gg + theme (axis.text.x = element_text(angle = 45, hjust = 1, size = 6)) 
    gg <- gg + labs(y="OccurenceRange", x="Freq") 
    ggsave(filename = paste0('results/Fig2_occurenceRange',systemFile,'_',i,'.pdf'))
    }
