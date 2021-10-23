# Data generation
By and large the whole dataset and process was written mainly in python's 2.7 version, i would be constantly updating and revisiting such that this process can be made resuable

# The figures and data
###### This panel overall will talk about the the components of the figure, which modules were mused form which file and othe overall layout of the system


The condition genes were set in the Fig1 subpanel B 


every core figure program will be associated with the list of common modules that it may use and folder having codes for the subpanels
This will help improve the overall readability of the code and maintennace
 

### Figure 1
#### Codes and generation:
The program is bin/creating_figures_python/Fig1_general_stats.py and this program's documnettaion is enclosed within and will be briefly discussed in every panel also, for every panel,
The program divides every sub image to differnt module, (bin/modules/figPanels/supanels_fig1) listed in subpanels_fig1
 
This figure gives a brief layout of the extent of the AS in humans atleast, 
We would give a detailed know how to the 
a) length of possible prteins, 
b) the transcripts types 
c) types of exons
d) WEF inclusion freqincy of the exons
e) length of exons and their types
f) the trancript types and their change in lnegth

Among these panels above i would like to keep only the panel C), D), and E) in main Figure 1

## Panel A) if we consider the average length of huamn protein and choose trancipt that is overa and large 

## panel C)







# The Meeting dated 23 oct 
For the record, i am supposed to write the fugure legends and their description till the 26th october for the description 1st and 2nd figure, and the NCB layouts, 
lets start the planning,
Figure 1 is complete and will finish the layout tonight only, 
For the fig 2, We were inteersted in the aspect  of looking the junction so secondary structres at the ends, 
We ficussed only on the junction types that were laways coding in all thetranscripts of the proteins, 


## Fig 2 Secondary structre and junctions
### Fig2.1 Preponderant junction types in the Alt-Alt Const Const and Const Alt types [All junction and Background wil lalso be represneted]
### Fig 2.2 Will make inteher classification and check the fequency of every such type being prep






# The general programming and the logs being uploaded 


python bin/creating_figures_python/Fig1_general_stats.py output/derived_data/results/objectsave_9606_0.6.pick results/

Rscript bin/plotting_figures_R/F1_plot_panelE.R results/Condition_3000Pilength_4Isf_4ExCount/panelE_Length_distribution_exonsRAW.csv results/Condition_3000Pilength_4Isf_4ExCount/
Rscript bin/plotting_figures_R/F1_plot_panelE.R results/Condition_3000Pilength_2Isf_2ExCount/panelE_Length_distribution_exonsRAW.csv results/Condition_3000Pilength_2Isf_2ExCount/

Rscript bin/plotting_figures_R/F1_plot_panelD.R results/Condition_3000Pilength_4Isf_4ExCount/Alternate_exons_WEF.csv results/Condition_3000Pilength_4Isf_4ExCount/
Rscript bin/plotting_figures_R/F1_plot_panelD.R results/Condition_3000Pilength_2Isf_2ExCount/Alternate_exons_WEF.csv results/Condition_3000Pilength_2Isf_2ExCount/

Rscript bin/plotting_figures_R/F1_plot_panelC.R results/Condition_3000Pilength_2Isf_2ExCount/Unique_Coods.csv results/Condition_3000Pilength_2Isf_2ExCount/
Rscript bin/plotting_figures_R/F1_plot_panelC.R results/Condition_3000Pilength_4Isf_4ExCount/Unique_Coods.csv results/Condition_3000Pilength_4Isf_4ExCount/


#figure2,   ss and plot
python bin/creating_figures_python/Fig2_ss_junctions.py output/derived_data/results/objectsave_9606_0.6.pick results/Condition_3000Pilength_2Isf_2ExCount/condition_genes.pick results/Condition_3000Pilength_2Isf_2ExCount/

python bin/creating_figures_python/Fig2_ss_junctions.py output/derived_data/results/objectsave_9606_0.6.pick results/Condition_3000Pilength_4Isf_4ExCount/condition_genes.pick results/Condition_3000Pilength_4Isf_4ExCount/


Rscript bin/plotting_figures_R/F2_plot_panelB0.R results/Condition_3000Pilength_4Isf_4ExCount/ _
Rscript bin/plotting_figures_R/F2_plot_panelB0.R results/Condition_3000Pilength_2Isf_2ExCount/ _
#secondary strucutre default

Rscript bin/plotting_figures_R/F2_plot_panelB1.R results/Condition_3000Pilength_2Isf_2ExCount/ WEF
Rscript bin/plotting_figures_R/F2_plot_panelB1.R results/Condition_3000Pilength_4Isf_4ExCount/ WEF
#secondary structure and WEF of junctions

Rscript bin/plotting_figures_R/F2_plot_panelB2.R results/Condition_3000Pilength_2Isf_2ExCount/ WEF
Rscript bin/plotting_figures_R/F2_plot_panelB2.R results/Condition_3000Pilength_4Isf_4ExCount/ WEF


#unrefined
Rscript bin/F2_plot_panelB1.R _BuriedExposedWin1WEF results/fig2/
Rscript bin/F2_plot_panelB1.R _BuriedExposedWin3WEF results/fig2/
Rscript bin/F2_plot_panelB1.R StrideAF_WEF results/fig2/
Rscript bin/F2_plot_panelB1.R WEF results/fig2/

Rscript bin/F2_plot_panelB2.R _BuriedExposedWin1WEF results/fig2/
Rscript bin/F2_plot_panelB2.R _BuriedExposedWin3WEF results/fig2/
Rscript bin/F2_plot_panelB2.R StrideAF_WEF results/fig2/
Rscript bin/F2_plot_panelB2.R WEF results/fig2/

#Rscript bin/pssm.R _BuriedExposedWin1 results/
Rscript bin/pssm.R _BuriedExposedWin3 results/
Rscript bin/pssm.R StrideAF_ results/