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
df1 = read.csv(paste0(source_dir,"Fig2_1_C_extchange_histogram.csv"))
df1$ChangeInLength <- factor(df1$ChangeInLength, levels = df1$ChangeInLength)
su<-sum(df1$Total)
p1<-ggplot(data=df1, aes(x=ChangeInLength, y=Frequency, fill=ChangeInLength)) + geom_bar(stat="identity") +
guides(fill=FALSE)+  geom_text(aes(label=Frequency), vjust=-0.3, color="black", size=3.5) +
labs(y="Frequency", x="Length_range", subtitle="Change_in_aa_length_in_Ccases") + theme (axis.text.x = element_text(angle = 45, hjust = 1, size = 10))


keeps_pos <- c("ChangeInLength", "Nter","Middle","Cter")
df2=df1[keeps_pos]
df2 <- melt(df2, id="ChangeInLength")
p2<-ggplot(data=df2, aes(x=ChangeInLength, y=value, fill=variable)) + geom_bar(stat="identity") + geom_text(aes(label=value),size = 3, position = position_stack(vjust = 0.5)) +
labs(y="Frequency", x="ChangeinLength", subtitle="Position_affected_wrt_changeinLength") + theme (axis.text.x = element_text(angle = 45, hjust = 1, size = 10), legend.position="bottom")

keeps_ss <- c("ChangeInLength", "Css","Hss","Ess")
df3=df1[keeps_ss]
df3 <- melt(df3, id="ChangeInLength")
p3<-ggplot(data=df3, aes(x=ChangeInLength, y=value, fill=variable)) + geom_bar(stat="identity") + geom_text(aes(label=value),size = 3, position = position_stack(vjust = 0.5)) +
labs(y="Frequency", x="ChangeinLength", subtitle="Prediction_SSchange_wrt_changeinLength") + theme (axis.text.x = element_text(angle = 45, hjust = 1, size = 10), legend.position="bottom")

keeps_st <- c("ChangeInLength", "Cst","Hst","Est")
df4=df1[keeps_st]
df4 <- melt(df4, id="ChangeInLength")
p4<-ggplot(data=df4, aes(x=ChangeInLength, y=value, fill=variable)) + geom_bar(stat="identity") + geom_text(aes(label=value),size = 3, position = position_stack(vjust = 0.5)) +
labs(y="Frequency", x="ChangeinLength", subtitle="Structre_SS_wrt_changeinLength") + theme (axis.text.x = element_text(angle = 45, hjust = 1, size = 10), legend.position="bottom")

df5 = read.csv(paste0(source_dir,"Fig2_1_C_extconserved_fraction.csv"))
df5$Conserved_sequence_range <- factor(df5$Conserved_sequence_range, levels = df5$Conserved_sequence_range)
p5<-ggplot(data=df5, aes(x=Conserved_sequence_range, y=Frequency, fill=Conserved_sequence_range)) + geom_bar(stat="identity") +
guides(fill=FALSE)+  geom_text(aes(label=Frequency), vjust=-0.3, color="black", size=3.5) +
labs(y="Frequency", subtitle="fraction_conserved_in_region_overlapping") + theme (axis.text.x = element_text(angle = 45, hjust = 1, size = 10))


plot2by2 <- plot_grid(p1, p2, p3, p4, p5, labels = "AUTO", ncol = 2, nrow=3)

save_plot(paste0(out_dir,"Fig2_1_C.png"), plot2by2,
          ncol = 2, # we're saving a grid plot of 2 columns
          nrow = 3, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 1.3
          )
#first histogras
