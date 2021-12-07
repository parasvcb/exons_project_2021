args = commandArgs(trailingOnly=TRUE)
if (length(args)<1) {
  stop("Please supply two arguements SourceDomDir", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  source_dir=args[1]
  out_dir=source_dir
}

library(ggplot2)
library(cowplot)
library(reshape2)
#ChangeInLength,TotalCases,Frequency,SSeqLength,Css,Hss,Ess,StrideLength,Cst,Hst,Est,Nter,Middle,Cter
df1n = read.csv(paste0(source_dir,"Fig2_1_BN_extchange_histogram.csv"))
df1n$ChangeInLength <- factor(df1n$ChangeInLength, levels = df1n$ChangeInLength)
su<-sum(df1n$Total)
p1n<-ggplot(data=df1n, aes(x=ChangeInLength, y=Frequency, fill=ChangeInLength)) + geom_bar(stat="identity") +
guides(fill=FALSE)+  geom_text(aes(label=Frequency), vjust=-0.3, color="black", size=2.5) +
labs(y="Frequency", x="Length_range", subtitle="aaChange_in_BNcases") + theme (axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

df1c = read.csv(paste0(source_dir,"Fig2_1_BC_extchange_histogram.csv"))
df1c$ChangeInLength <- factor(df1c$ChangeInLength, levels = df1c$ChangeInLength)
su<-sum(df1c$Total)
p1c<-ggplot(data=df1c, aes(x=ChangeInLength, y=Frequency, fill=ChangeInLength)) + geom_bar(stat="identity") +
guides(fill=FALSE)+  geom_text(aes(label=Frequency), vjust=-0.3, color="black", size=2.5) +
labs(y="Frequency", x="Length_range", subtitle="aaChange_in_BCcases") + theme (axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
#========================================
keeps_pos <- c("ChangeInLength", "Nter","Middle","Cter")
df2n=df1n[keeps_pos]
df2n <- melt(df2n, id="ChangeInLength")
p2n<-ggplot(data=df2n, aes(x=ChangeInLength, y=value, fill=variable)) + geom_bar(stat="identity") + geom_text(aes(label=value),size = 2, position = position_stack(vjust = 0.5)) +
labs(y="Frequency", x="ChangeinLength", subtitle="Position_affected") + theme (axis.text.x = element_text(angle = 45, hjust = 1, size = 8), legend.position="bottom")

df2c=df1c[keeps_pos]
df2c <- melt(df2c, id="ChangeInLength")
p2c<-ggplot(data=df2c, aes(x=ChangeInLength, y=value, fill=variable)) + geom_bar(stat="identity") + geom_text(aes(label=value),size = 2, position = position_stack(vjust = 0.5)) +
labs(y="Frequency", x="ChangeinLength", subtitle="Position_affected") + theme (axis.text.x = element_text(angle = 45, hjust = 1, size = 8), legend.position="bottom")
#=========================================

keeps_ss <- c("ChangeInLength", "Css","Hss","Ess")
df3n=df1n[keeps_ss]
df3n <- melt(df3n, id="ChangeInLength")
p3n<-ggplot(data=df3n, aes(x=ChangeInLength, y=value, fill=variable)) + geom_bar(stat="identity") + geom_text(aes(label=value),size = 2, position = position_stack(vjust = 0.5)) +
labs(y="Frequency", x="ChangeinLength", subtitle="Prediction_SSchange") + theme (axis.text.x = element_text(angle = 45, hjust = 1, size = 8), legend.position="bottom")

keeps_ss <- c("ChangeInLength", "Css","Hss","Ess")
df3c=df1c[keeps_ss]
df3c <- melt(df3c, id="ChangeInLength")
p3c<-ggplot(data=df3c, aes(x=ChangeInLength, y=value, fill=variable)) + geom_bar(stat="identity") + geom_text(aes(label=value),size = 2, position = position_stack(vjust = 0.5)) +
labs(y="Frequency", x="ChangeinLength", subtitle="Prediction_SSchange") + theme (axis.text.x = element_text(angle = 45, hjust = 1, size = 8), legend.position="bottom")
 #=====================================================

if (FALSE) {
keeps_st <- c("ChangeInLength", "Cst","Hst","Est")
df4c=df1c[keeps_st]
df4c <- melt(df4c, id="ChangeInLength")
p4c<-ggplot(data=df4c, aes(x=ChangeInLength, y=value, fill=variable)) + geom_bar(stat="identity") + geom_text(aes(label=value),size = 2, position = position_stack(vjust = 0.5)) +
labs(y="Frequency", x="ChangeinLength", subtitle="Structre_SS") + theme (axis.text.x = element_text(angle = 45, hjust = 1, size = 8), legend.position="bottom")

keeps_st <- c("ChangeInLength", "Cst","Hst","Est")
df4n=df1n[keeps_st]
df4n <- melt(df4n, id="ChangeInLength")
p4n<-ggplot(data=df4n, aes(x=ChangeInLength, y=value, fill=variable)) + geom_bar(stat="identity") + geom_text(aes(label=value),size = 2, position = position_stack(vjust = 0.5)) +
labs(y="Frequency", x="ChangeinLength", subtitle="Structre_SS") + theme (axis.text.x = element_text(angle = 45, hjust = 1, size = 8), legend.position="bottom")
}
#=============================================

df5 = read.csv(paste0(source_dir,"Fig2_1_BN_extconserved_fraction.csv"))
df5$Conserved_sequence_range <- factor(df5$Conserved_sequence_range, levels = df5$Conserved_sequence_range)
p5<-ggplot(data=df5, aes(x=Conserved_sequence_range, y=Frequency, fill=Conserved_sequence_range)) + geom_bar(stat="identity") +
guides(fill=FALSE)+  geom_text(aes(label=Frequency), vjust=-0.3, color="black", size=3.5) +
labs(y="Frequency", subtitle="fraction_conserved_in_region_overlapping") + theme (axis.text.x = element_text(angle = 45, hjust = 1, size = 8))



plot2by2 <- plot_grid(p1n,p2n, p3n, p5, labels = "AUTO", ncol = 2, nrow=2)

save_plot(paste0(out_dir,"Fig2_1_B_N.png"), plot2by2,
          ncol = 2, # we're saving a grid plot of 2 columns
          nrow = 2, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 1.5
          )

plot2by2 <- plot_grid(p1c,p2c,p3c, p5, labels = "AUTO", ncol = 2, nrow=2)
save_plot(paste0(out_dir,"Fig2_1_B_C.png"), plot2by2,
          ncol = 2, # we're saving a grid plot of 2 columns
          nrow = 2, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 1.5
          )

#first histogras
