args = commandArgs(trailingOnly=TRUE)
if (length(args)<1) {
  stop("Please supply two arguements SourceDomDir", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  source_dir=args[1]
  out_dir=source_dir
}

library(cowplot)
library(ggplot2)

df1 = read.csv(paste0(source_dir,"Fig2_1_N_extdensityvals.csv"))
p1<-ggplot(data=df1, aes(x=position)) + geom_density(color="darkblue", fill="lightblue") + 
scale_y_continuous(breaks=seq(0,100,10)) +
labs(x="Protein_position", subtitle="N_cases_affecting_princiapl_isoform")

df2 = read.csv(paste0(source_dir,"Fig2_1_C_extdensityvals.csv"))
p2<-ggplot(data=df2, aes(x=position)) + geom_density(color="darkblue", fill="lightblue") + 
scale_y_continuous(breaks=seq(0,100,10)) +
labs(x="Protein_position", subtitle="C_cases_affecting_princiapl_isoform")

df3 = read.csv(paste0(source_dir,"Fig2_1_BC_extdensityvals.csv"))
p3<-ggplot(data=df3, aes(x=position)) + geom_density(color="darkblue", fill="lightblue") + 
scale_y_continuous(breaks=seq(0,100,10)) +
labs(x="Protein_position", subtitle="B_cases_affecting_princiapl_isoform")


plot2by2 <- plot_grid(p1, p2, p3, labels = "AUTO", ncol = 3, nrow=1)

save_plot(paste0(out_dir,"Fig2_1_densityNCB.png"), plot2by2,
          ncol = 3, # we're saving a grid plot of 2 columns
          nrow = 1, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 1.3
          )
