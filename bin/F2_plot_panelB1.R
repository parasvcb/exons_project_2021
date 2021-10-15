args = commandArgs(trailingOnly=TRUE)
if (length(args)<2) {
  stop("Please supply two arguements 1. inpPanelsString and 2.outputFile", call.=FALSE)
} else if (length(args)==2) {
  # default output file
  systemFile=args[1]
  out=args[2]
}
#ExonsJunction_BuriedExposedWin1WEF_and_exonsJunction.csv

library('ggplot2')

df3<-read.csv(paste0("results/ExonsJunction",systemFile,"_and_exonsJunction.csv"),sep='\t')
df3$wrap <- paste(df3$ExonType,df3$WEFRange)
print (head (df3))
print (unique(df3$wrap))
print (paste("cons_cons 0.25_0.5",sum(df3[df3$wrap=="cons_cons 0.25_0.5",]$Count)))

print (sum(df3[df3$wrap=='cons_cons 0_0.25',]$Count))
hospital_names <- c(
                    `cons_cons 0_0.25` = paste("cons_cons 0_0.25",sum(df3[df3$wrap=="cons_cons 0_0.25",]$Count)),
                    `cons_cons 0.25_0.5` = paste("cons_cons 0.25_0.5",sum(df3[df3$wrap=="cons_cons 0.25_0.5",]$Count)),
                    `cons_cons 0.5_0.75` = paste("cons_cons 0.5_0.75",sum(df3[df3$wrap=="cons_cons 0.5_0.75",]$Count)),
                    `cons_cons 0.75_1.01` = paste("cons_cons 0.75_1.01",sum(df3[df3$wrap=="cons_cons 0.75_1.01",]$Count)),
                    
                    `alt_alt 0_0.25` = paste("alt_alt 0_0.25",sum(df3[df3$wrap=="alt_alt 0_0.25",]$Count)),
                    `alt_alt 0.25_0.5` = paste("alt_alt 0.25_0.5",sum(df3[df3$wrap=="alt_alt 0.25_0.5",]$Count)),
                    `alt_alt 0.5_0.75` = paste("alt_alt 0.5_0.75",sum(df3[df3$wrap=="alt_alt 0.5_0.75",]$Count)),
                    `alt_alt 0.75_1.01` = paste("alt_alt 0.75_1.01",sum(df3[df3$wrap=="alt_alt 0.75_1.01",]$Count)),
                    
                    `cons_alt 0_0.25` = paste("cons_alt 0_0.25",sum(df3[df3$wrap=="cons_alt 0_0.25",]$Count)),
                    `cons_alt 0.25_0.5` = paste("cons_alt 0.25_0.5",sum(df3[df3$wrap=="cons_alt 0.25_0.5",]$Count)),
                    `cons_alt 0.5_0.75` = paste("cons_alt 0.5_0.75",sum(df3[df3$wrap=="cons_alt 0.5_0.75",]$Count)),
                    `cons_alt 0.75_1.01` = paste("cons_alt 0.75_1.01",sum(df3[df3$wrap=="cons_alt 0.75_1.01",]$Count))
                    )
df3_1 = df3[!df3$Type=='Extras',]
df3_2 = df3_1[!df3_1$Exons<1,]

#ExonType	WEFRange	Junction	Count	Freq
#cons_cons	0_0.25	HX	35	0.012
gg=ggplot(df3, aes(x=Junction,y=Freq, fill=Freq)) + geom_bar(stat="identity") 
gg <- gg + facet_wrap(~wrap,nrow=3,labeller = as_labeller(hospital_names))
gg <- gg + guides(fill=FALSE)+  geom_text(aes(label=Freq), vjust=-0.3, color="black", size=1.6) 
gg <- gg + theme (axis.text.x = element_text(angle = 45, hjust = 1, size = 6)) 
gg <- gg + labs(y="SS Type", x="Freq") 

ggsave(filename = paste0('results/Fig2_WEF_andSS',systemFile,'.pdf'), width=10, height= 10)


