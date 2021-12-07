# library(ggplot2)
args = commandArgs(trailingOnly=TRUE)
if (length(args)<2) {
  stop("Please supply two arguements SourceDomDir and PlotsOutDir", call.=FALSE)
} else if (length(args)==2) {
  # default output file
  source_dir=args[1]
  out_dir=args[2]
}
library(cowplot)
print ("!")
#ggplot(data=df,aes(x=Position_affected,y=Freq, group=1))+geom_line()+geom_point()
df1 = read.csv(paste0(source_dir,"N_extdensityvals_CDF.csv")) 
p1<-ggplot(data=df1,aes(x=Position_affected,y=Freq, group=1))+geom_line()+geom_point() +
labs(x="Protein_position", subtitle="N_cases_affecting_princiapl_isoform")
print ("2")
df2 = read.csv(paste0(source_dir,"C_extdensityvals_CDF.csv")) 
p2<-ggplot(data=df2,aes(x=Position_affected,y=Freq, group=1))+geom_line()+geom_point() +
labs(x="Protein_position", subtitle="C_cases_affecting_princiapl_isoform")
print ("3")
df3 = read.csv(paste0(source_dir,"BC_extdensityvals_CDF.csv"))
p3<-ggplot(data=df3,aes(x=Position_affected,y=Freq, group=1))+geom_line()+geom_point() +
labs(x="Protein_position", subtitle="B_cases_affecting_princiapl_isoform")
print ("4")

plot2by2 <- plot_grid(p1, p2, p3, labels = "AUTO", ncol = 3, nrow=1)

save_plot(paste0(out_dir,"fig2_density_CDF.png"), plot2by2,
          ncol = 3, # we're saving a grid plot of 2 columns
          nrow = 1, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 1.3
          )
