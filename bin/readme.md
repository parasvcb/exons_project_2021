v: 27.09.22

Two new programs has been added,
prep.bash (which takes into account making preparations and standardizing directory structure)
    improvements needed:
    A) taking stockpiling for the current genome and the tax-ids, on broader level , it should ask for the modes, NCBI/ENSEMBL and the list of tax-Ids instead of working for the one, 
    B) every set should be given the separate directories, with record stamp of the date on the parent folder.
    C) all prelim file structure should not be created here but in the python program, as storing those addresses again will have the complications.
    D) use it to keep record of scripts used to derived data, PFAM, I-tasser (PSIPRED), PDB sift file, cross reference files, alphafold databases, colabold DB's

object builder
the most recent file is the object_builder_latest.py
changes:
    structure information is made optional, the current DB isn't harboring any of those releases.
    the nomenclature file has been changes (IR event nomenclature are recorded in middle for start of IR event from former exon but not from the last one [yes it was for last one])
    F event cases has been improved and I think in future releases I will remove the usage of the redefine F cases module
    Extensive documentation has been done for the many inside modules but still a detailed chart out is pending
    after some releases, the structure and cath information modules should be discarded.
    

STATUS:
F CASES DONE, IR EVENT CASE IMPROVED, EVENTS ARE BIT LIGHT ON COMPUTATION, DISORDER IS BEING MEASURED.

PENDING
FOREVER IMPROVEMENTS
Documentation of every module there and here is a must, 
ENSEMBL ready build

checkup of the latest to build colabfold files





# Data generation
By and large the whole dataset and process was written mainly in python's 2.7 version, i would be constantly updating and revisiting such that this process can be made resuable
(the biggest hurdle is pickling objects)


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

The main figure is decides, need to writ etheir methodlogy and where to invoke trhe modified detauils of the 

## Fig 2 Secondary structre and junctions
This whole section is being scrutinized again, 
The aspect 1, only T cases, 
aspect 2, T cases between first Coding exons and constitutive exons
aspect 3, T cases exlcuding
### Fig2.1 Preponderant junction types in the Alt-Alt Const Const and Const Alt types [All junction and Background wil lalso be represneted]
### Fig 2.2 Will make WEF based junction classification and check the fequency occurence of this junction to be maintained often
#### Fig 2.2.A 
   '''
    Screen all transcripts (by default they will be protein coding)
    fill the has (repr)
        iterate its exons in pairs, 
        if their length is grt 9 (both):
            get their ss seqs
            if they are T cases:
                get their linear sequqnce tags or close to them ()
                now get their junction sequqnces, both amino acid and SS (4tup)[It can be empty also or undefined, However this helps getting WEF count correct]
                make a key based on the linear tag and add a list [PI,coding_exons_count,texo[i].ID,texo[i+1].ID,SS_junction]

    itrerate the repr: key is pair of linear tag and values are the jujctions with different proerties here
    for pair in repr:
        get count of different times this pair has been called out for with defined sequqnce [not NA or doscard unsolved]
        if len(representative[exonSequencePair])>2 and numberofssNotfalse >2:
            get junction representative: sort the repr has and screen till we get a repr,
            if juncrepr:
                get its WEF (from all the junctions (solved unsolved))
                get other unique junctions also (wef unsolved and solved)
                this is in form of 
                has_exon_categories_combined[flagPairKey][gene,exonSequencePair,WEF]=(representative[exonSequencePair],reprSS,reprSS_count,reprSS_Freq, temphas)
                where temphas temphas[tempjunction]=round(float(temphas[tempjunction])/len(representative[exonSequencePair]),3)

    
    now this has_exon_categories_combined is send for being sorted screening

    in short (), these junctions wef will be segregated into 4 sets, ragings in steps of 0.25 from 0 to 1, and their represnetativeis shown in panel 2.1
    and all soncsetrvations from other pissble SS it may adoopt willbe showin in 2.2

    in detail()

    f1=open(fname_junctions,'w')
    f1.write('ExonType\tWEFRange\tJunction\tCount\tFreq\n')
    conservation_has={exonType:{} for exonType in has} 
    # this will be used to get data for the cvonservation/occurence of the Junctions in the isoforms, segregated first into exonTypes, later into their occurenec range (HH is observed in 23 times out of 100 when junction type is ALTALt) and then in to HH, how oftne this was presnerved in next junctions
    for exonType in has:
        #exonType is conscons altalt altcons
        WEFRangeHas={(0,0.25):{},(0.25,0.50):{},(0.50,0.75):{},(0.75,1.01):{}}
        #seg 1 begins
        for geneExpairWefKey in has[exonType]:
            #geneExpairWefKey is tuple of (gene,Expair,Wef)
            wef=geneExpairWefKey[2]
            exonsListdetailed,junction,junctionCount,JunctionFreq,temphas=has[exonType][geneExpairWefKey]
            #print (exonsListdetailed)
            #print (junction,junctionCount,JunctionFreq,temphas)
            for psiExonRange in WEFRangeHas:
                #print (psiExonRange,wef)
                if psiExonRange[0]<wef<=psiExonRange[1]:
                    #print ('yes')
                    if junction not in WEFRangeHas[psiExonRange]:
                        WEFRangeHas[psiExonRange][junction]=[0,[]]
                    WEFRangeHas[psiExonRange][junction][0]+=1
                    WEFRangeHas[psiExonRange][junction][1]+=[JunctionFreq]
                    break
            # seg1, WEFRangeHas updated
            if psiExonRange not in conservation_has[exonType]:
                conservation_has[exonType][psiExonRange]={}
            # seg2, initiated, we are still in loop 2 (iterating major exjunc types) then individual pairs in them, 
            # where we are considering adding PSI value of junction to anotherhas which will calculate the occurence of major HHTYPES in them 
            for junctionAll in temphas: 
                #print (junctionAll)
                #print (temphas[junctionAll])   
                if junctionAll in default_sstypes:
                    #seg 2 now being pro[perly resumed]
                    if junctionAll not in conservation_has[exonType][psiExonRange]:
                        conservation_has[exonType][psiExonRange][junctionAll]={(0,0.10):0,(0.101,0.20):0,(0.201,0.30):0,(0.301,0.40):0,(0.401,0.50):0,(0.501,0.60):0,(0.601,0.70):0,(0.701,0.80):0,(0.801,0.90):0,(0.901,1.00):0}
                    for vals2 in conservation_has[exonType][psiExonRange][junctionAll]:
                        if vals2[0]<=temphas[junctionAll]<=vals2[1]:
                            conservation_has[exonType][psiExonRange][junctionAll][vals2]+=1
                            break
        # seg2, is now updated
        # resuming seg1 below
        list_WEF=list(WEFRangeHas.keys())
        list_WEF.sort()
        for WEFrange in list_WEF:
            #f1.write('ExonType\tWEFRange\tJunction\tCount\tFreq\n') written already
            total_junctions=sum([WEFRangeHas[WEFrange][junc][0] for junc in WEFRangeHas[WEFrange]])
            for junctions in WEFRangeHas[WEFrange]:
                count=WEFRangeHas[WEFrange][junctions][0]
                f1.write('%s\t%s\t%s\t%s\t%s\n'%(exonType,"_".join(map(str,list(WEFrange))),junctions,count,gm.div_fact(count,total_junctions)))
    #seg 1 done
    #seg 2, now lets do the second aspect of calculating the fraction
    f2=open(fname_junc_conservations,'w')
    f2.write('ExonType\tWEFRange\tJunction\tOccurenceRange\tOccureneceTotal\tOccureneceFreq\tcumulativeCount\tcumulativeFreq\n')
    a=conservation_has.keys()
    # for key in a:
    #     print (key,conservation_has[key].keys())
    
    #print (conservation_has)
    for exonType in conservation_has:
        for WEFRangeJunction in conservation_has[exonType]:
            wefrangestr="_".join(map(str,list(WEFRangeJunction)))
            for junction in conservation_has[exonType][WEFRangeJunction]:
                totalOcc=sum([conservation_has[exonType][WEFRangeJunction][junction][ti] for ti in conservation_has[exonType][WEFRangeJunction][junction]])
                occRangeList=list(conservation_has[exonType][WEFRangeJunction][junction].keys())
                occRangeList.sort()
                previousValue=0
                for occRange in occRangeList:
                    occrangestr="_".join(map(str,list(occRange)))
                    count_occrange=conservation_has[exonType][WEFRangeJunction][junction][occRange]
                    cumulative_count=previousValue+count_occrange
                    previousValue+=count_occrange
                    #print (1,exonType,WEFRangeJunction)
                    #print (2,exonType,wefrangestr,junction,occrangestr,count_occrange,gm.div_fact(count_occrange,totalOcc), cumulative_count,gm.div_fact(cumulative_count,totalOcc))
                    f2.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(exonType,wefrangestr,junction,occrangestr,count_occrange,gm.div_fact(count_occrange,totalOcc), cumulative_count,gm.div_fact(cumulative_count,totalOcc)))

## NCB Cases
ATIT:Same_amino_acid
Means differeing non docing cooridnates but similar coding coordinates


# The general programming and the logs being uploaded 


python bin/creating_figures_python/Fig1_general_stats.py output/derived_data/results/objectsave_9606_0.6.pick results/

Rscript bin/plotting_figures_R/F1_plot_panelE.R results/Condition_3000Pilength_4Isf_4ExCount/panelE_Length_distribution_exonsRAW.csv results/Condition_3000Pilength_4Isf_4ExCount/
Rscript bin/plotting_figures_R/F1_plot_panelE.R results/Condition_3000Pilength_2Isf_2ExCount/panelE_Length_distribution_exonsRAW.csv results/Condition_3000Pilength_2Isf_2ExCount/

Rscript bin/plotting_figures_R/F1_plot_panelD.R results/Condition_3000Pilength_4Isf_4ExCount/Alternate_exons_WEF.csv results/Condition_3000Pilength_4Isf_4ExCount/
Rscript bin/plotting_figures_R/F1_plot_panelD.R results/Condition_3000Pilength_2Isf_2ExCount/Alternate_exons_WEF.csv results/Condition_3000Pilength_2Isf_2ExCount/

Rscript bin/plotting_figures_R/F1_plot_panelC.R results/Condition_3000Pilength_2Isf_2ExCount/Unique_Coods.csv results/Condition_3000Pilength_2Isf_2ExCount/
Rscript bin/plotting_figures_R/F1_plot_panelC.R results/Condition_3000Pilength_4Isf_4ExCount/Unique_Coods.csv results/Condition_3000Pilength_4Isf_4ExCount/


# figure2,   ss and plot
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



#fig2 
python bin/creating_figures_python/Fig2_1NCBI_file_writer_core.py output/derived_data/results/objectsave_9606_0.6.pick results/Condition_3000Pilength_2Isf_2ExCount/
python bin/creating_figures_python/Fig2_1NCBI_file_writer_core.py output/derived_data/results/objectsave_9606_0.6.pick results/Condition_3000Pilength_4Isf_4ExCount/


python bin/creating_figures_python/figureNCBandA.py results/Condition_3000Pilength_2Isf_2ExCount/
python bin/creating_figures_python/figureNCBandA.py results/Condition_3000Pilength_4Isf_4ExCount/

Rscript bin/plotting_figures_R/F2_1_N.R results/Condition_3000Pilength_2Isf_2ExCount/
Rscript bin/plotting_figures_R/F2_1_C.R results/Condition_3000Pilength_2Isf_2ExCount/
Rscript bin/plotting_figures_R/F2_1_B.R results/Condition_3000Pilength_2Isf_2ExCount/
Rscript bin/plotting_figures_R/F2_1_density_plot.R results/Condition_3000Pilength_2Isf_2ExCount/

Rscript bin/plotting_figures_R/F2_1_N.R results/Condition_3000Pilength_4Isf_4ExCount/
Rscript bin/plotting_figures_R/F2_1_C.R results/Condition_3000Pilength_4Isf_4ExCount/
Rscript bin/plotting_figures_R/F2_1_B.R results/Condition_3000Pilength_4Isf_4ExCount/
Rscript bin/plotting_figures_R/F2_1_density_plot.R results/Condition_3000Pilength_4Isf_4ExCount/




# fig4



python bin/creating_figures_python/Fig3_domains_analyse.py results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/Fig3D_exonsWise_atit.pickle results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/Fig3D_domain_annotation_atit.pickle output/genetoPDB.pickle results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/


Total_genes=6988, included cases=6482, total_aa=4667491, total_aa_dom=2095440, total_aa_dom_middle=629658, fracCoveerd_wdom=0.449, fraccoveerd_m_dom=0.135







# fig 2, all cases
(py27) paras@menten:~/exonsdrive/paras/project/protein_splicing/projectDir$ python bin/creating_figures_python/Fig2_ss_junctions_TempTjunctions.py output/derived_data/results/objectsave_9606_0.6.pick results/Condition_3000Pilength_4Isf_4ExCount/condition_genes.pick results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/ -1
(py27) paras@menten:~/exonsdrive/paras/project/protein_splicing/projectDir$ python bin/creating_figures_python/Fig2_ss_junctions_TempTjunctions.py output/derived_data/results/objectsave_9606_0.6.pick results/Condition_3000Pilength_4Isf_4ExCount/condition_genes.pick results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/ 1
(py27) paras@menten:~/exonsdrive/paras/project/protein_splicing/projectDir$ python bin/creating_figures_python/Fig2_ss_junctions_TempTjunctions.py output/derived_data/results/objectsave_9606_0.6.pick results/Condition_3000Pilength_4Isf_4ExCount/condition_genes.pick results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/ 0

Rscript bin/plotting_figures_R/F2_plot_panelB0.R results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/ _ atit
Rscript bin/plotting_figures_R/F2_plot_panelB0.R results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/ _ both
Rscript bin/plotting_figures_R/F2_plot_panelB0.R results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/ _ core


Rscript bin/plotting_figures_R/F2_plot_panelB1.R results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/ WEF atit
Rscript bin/plotting_figures_R/F2_plot_panelB1.R results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/ WEF both
Rscript bin/plotting_figures_R/F2_plot_panelB1.R results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/ WEF core

Rscript bin/plotting_figures_R/F2_plot_panelB2.R results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/ WEF atit
Rscript bin/plotting_figures_R/F2_plot_panelB2.R results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/ WEF both
Rscript bin/plotting_figures_R/F2_plot_panelB2.R results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/ WEF core

# atiCore
python bin/creating_figures_python/fraction_protein_atitcore.py output/derived_data/results/objectsave_9606_0.6.pick results/Condition_3000Pilength_4Isf_4ExCount/condition_genes.pick results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/

Rscript bin/plotting_figures_R/F0_dataGeneral.R results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/ATICoreExonsProteinFraction.csv results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/





exondomain contribution;

its importnat to know that i ma iterating firs the domains then the exons, 

so its wise to question the impcat of the domain fraction per exon, as its primarily looped, 
    so for every dfomain in question we have the reposnse here, 

so query the exon fraction per domain, we need to be sure if we want o calculate that per domain wise or the whole domain wise lets do it also


The fraction of domains when undergo split will only be considered when domain full domain is atleast 95% covered in tersn of residues, [will affect only the interafce domains of the ATI and ATT region]


iterate the genes [condition 4 isoforms and 4 exons]:
    compute pergeneDomainFraction [histogram], totalgeneDomainFraction [number]
    if domains and exonslis:
        iterate the domains:
            ietarte the exons:
                check interesection and divide it bydomainspan to getdoain span,
                if95> this domain is contained
                repeat for middle region also

                chcek containedDomains
                write the exons and domain fractions into the system, 

                if domain is 95% coevered record the instance, 

        use the has above and see if you can gain the 
        to get exons fraction out, [it will be goodf to use the domain if the fraction of whole is conatined for 95%, it will help aining proxy and middle will also be the 100pcnt, now the contratry exon fraction will be tricky and lets do it]



no  domain, isoform for gene 56
6
output/derived_data/results
['NP_064493.1', 'NP_064492.1', 'NP_064454.1', 'NP_001603.1']

gene 98 is also inetesrting
('XP_016860902.1', ['U.-2.A.1.n.1', 'T.1.A.2.0.0', 'T.1.A.3.0.0', 'T.1.A.5.0.0', 'T.1.A.6.0.0', 'T.1.G.9.0.0', 'T.2.G.10.0.0', 'T.1.A.13.0.0', 'D.-2.A.15.c.1'])
('XP_016860901.1', ['U.-2.A.1.n.1', 'T.1.A.2.0.0', 'T.1.A.3.0.0', 'T.1.A.5.0.0', 'T.1.A.6.0.0', 'T.1.G.9.0.0', 'T.1.G.10.0.0', 'T.1.A.14.0.0', 'D.-2.A.15.c.1'])
('XP_016860900.1', ['U.-2.A.1.n.1', 'T.1.A.2.0.0', 'T.1.A.3.0.0', 'T.1.A.5.0.0', 'T.1.A.6.0.0', 'T.1.G.9.0.0', 'T.1.G.10.0.0', 'D.1.A.15.c.1'])
('XP_016860895.1', ['U.-2.A.1.n.1', 'T.1.A.2.0.0', 'T.1.A.3.0.0', 'T.1.A.5.0.0', 'T.1.A.6.0.0', 'T.1.G.9.0.0', 'T.1.G.10.0.0', 'T.1.A.11.0.0'])
('XP_016860896.1', ['U.-2.A.1.n.1', 'T.1.A.2.0.0', 'T.1.A.5.0.0', 'T.1.A.6.0.0', 'T.1.G.9.0.0', 'T.1.G.10.0.0', 'T.1.A.11.0.0'])
('XP_016860898.1', ['U.-2.A.1.n.1', 'T.2.A.2.0.0', 'T.1.A.6.0.0', 'T.1.G.9.0.0', 'T.1.G.10.0.0', 'T.1.A.11.0.0'])
('NP_001307515.1', ['U.-2.A.1.0.0', 'T.1.A.2.0.0', 'T.1.A.5.0.0', 'T.1.A.6.0.0', 'T.1.G.9.0.0', 'T.1.G.10.0.0', 'D.1.A.15.0.0'])
('NP_001307516.1', ['U.-2.A.1.0.0', 'T.2.A.2.0.0', 'T.1.A.6.0.0', 'T.1.G.9.0.0', 'T.1.G.10.0.0', 'D.1.A.15.0.0'])
('XP_016860897.1', ['U.-2.A.4.0.0', 'T.1.A.5.c.1', 'T.1.A.6.0.0', 'T.1.G.9.0.0', 'T.1.G.10.0.0', 'T.1.A.11.0.0'])
('NP_001307518.1', ['T.1.A.7.n.1', 'T.2.G.9.0.0', 'T.1.G.10.0.0', 'T.1.A.14.0.0', 'D.-2.A.15.0.0'])
('NP_612457.1', ['T.1.A.7.n.1', 'T.2.G.9.0.0', 'T.1.G.10.0.0', 'D.1.A.15.0.0'])
('NP_001307519.1', ['T.1.A.7.n.1', 'T.2.G.9.0.0', 'T.1.G.10.0.0', 'T.1.A.12.0.0'])
('XP_016860899.1', ['T.1.A.7.0.0', 'T.2.G.9.0.0', 'T.1.G.10.0.0', 'T.1.A.11.0.0'])
('NP_001307517.1', ['U.-2.A.8.0.0', 'T.3.G.9.0.0', 'T.1.G.10.0.0', 'D.1.A.15.0.0'])


5328
another interesting case:
('XP_011538168.1', ['U.-2.A.1.0.0', 'D.1.F.3.0.0', 'T.1.A.4.0.0', 'D.1.G.5.0.0', 'T.1.G.6.0.0', 'T.2.G.7.0.0', 'T.1.G.8.0.0', 'T.1.G.9.0.0', 'T.1.G.10.0.0', 'T.1.G.11.0.0', 'T.1.F.12.c.1'])
('NP_001138503.1', ['U.-2.A.2.0.0', 'T.1.F.3.c.1', 'D.3.G.5.0.0', 'T.1.G.6.0.0', 'T.1.G.7.0.0', 'T.1.G.8.0.0', 'T.1.G.9.0.0', 'T.1.G.10.0.0', 'T.1.G.11.0.0', 'T.1.F.12.0.0'])
('NP_002649.1', ['U.-2.A.2.0.0', 'D.1.F.3.0.0', 'T.1.A.4.0.0', 'D.1.G.5.0.0', 'T.1.G.6.0.0', 'T.1.G.7.0.0', 'T.1.G.8.0.0', 'T.1.G.9.0.0', 'T.1.G.10.0.0', 'T.1.G.11.0.0', 'T.1.F.12.0.0'])
('NP_001306120.1', ['U.-2.A.2.0.0', 'D.-2.F.3.0.0', 'D.-2.G.5.0.0', 'T.2.G.6.0.0', 'T.1.G.7.0.0', 'T.1.G.8.0.0', 'T.1.G.9.0.0', 'T.1.G.10.0.0', 'T.1.G.11.0.0', 'T.1.F.12.0.0'])
('NP_002649.1', 'U.-2.A.2.0.0', 0)
A
('NP_002649.1', 'D.1.F.3.0.0', 19)
G
('NP_002649.1', 'T.1.A.4.0.0', 9)
A
('NP_002649.1', 'D.1.G.5.0.0', 36)
G
('NP_002649.1', 'T.1.G.6.0.0', 59)
G
('NP_002649.1', 'T.1.G.7.0.0', 30)
G
('NP_002649.1', 'T.1.G.8.0.0', 74)
G
('NP_002649.1', 'T.1.G.9.0.0', 49)
G
('NP_002649.1', 'T.1.G.10.0.0', 47)
G
('NP_002649.1', 'T.1.G.11.0.0', 50)
G
('NP_002649.1', 'T.1.F.12.0.0', 58)
G
[[19, 'G', 'D.1.F.3.0.0', 1, 0], [9, 'A', 'T.1.A.4.0.0', 1, 0], [36, 'G', 'D.1.G.5.0.0', 1, 0], [59, 'G', 'T.1.G.6.0.0', 1, 0], [30, 'G', 'T.1.G.7.0.0', 1, 0], [74, 'G', 'T.1.G.8.0.0', 1, 0], [49, 'G', 'T.1.G.9.0.0', 1, 0], [47, 'G', 'T.1.G.10.0.0', 1, 0], [50, 'G', 'T.1.G.11.0.0', 1, 0], [58, 'G', 'T.1.F.12.0.0', 1, 0]]
[((70, 151), 'PF00051.18', 1.04, 'Kringle', 'Domain'), ((179, 419), 'PF00089.26', 1.09, 'Trypsin', 'Domain')]


python bin/creating_figures_python/Fig3_domains_section.py output/derived_data/results/objectsave_9606_0.6.pick results/Condition_3000Pilength_4Isf_4ExCount/condition_genes.pick  results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/Fig3D_ pfam 0
output/derived_data/results
Total_genes=6988, included cases=6482, total_aa=2374186, total_aa_dom=981563, total_aa_dom_middle=272600, fracCoveerd_wdom=0.413, fraccoveerd_m_dom=0.115
atit

python bin/creating_figures_python/Fig3_domains_section.py output/derived_data/results/objectsave_9606_0.6.pick results/Condition_3000Pilength_4Isf_4ExCount/condition_genes.pick  results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/Fig3D_ pfam 1
#core
output/derived_data/results
Total_genes=6988, included cases=6482, total_aa=2293305, total_aa_dom=1113877, total_aa_dom_middle=357058, fracCoveerd_wdom=0.486, fraccoveerd_m_dom=0.156

python bin/creating_figures_python/Fig3_domains_section.py output/derived_data/results/objectsave_9606_0.6.pick results/Condition_3000Pilength_4Isf_4ExCount/condition_genes.pick  results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/Fig3D_ pfam -1
Total_genes=6988, included cases=6482, total_aa=4667491, total_aa_dom=2095440, total_aa_dom_middle=629658, fracCoveerd_wdom=0.449, fraccoveerd_m_dom=0.135

<!-- output/derived_data/results
Total_genes=6988, included cases=6482, total_aa=2603823, total_aa_dom=1074657, total_aa_dom_middle=303717, fracCoveerd_wdom=0.413, fraccoveerd_m_dom=0.117
atit -->
<!-- output/derived_data/results
Total_genes=6988, included cases=6482, total_aa=2063668, total_aa_dom=1020783, total_aa_dom_middle=325941, fracCoveerd_wdom=0.495, fraccoveerd_m_dom=0.158
core -->
<!-- output/derived_data/results
Total_genes=6988, included cases=6482, total_aa=4667491, total_aa_dom=2095440, total_aa_dom_middle=629658, fracCoveerd_wdom=0.449, fraccoveerd_m_dom=0.135
both -->

python bin/creating_figures_python/Fig3_domains_analyse.py results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/Fig3D_domain_exon_annotation_atit.pickle results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/Fig3D_

python bin/creating_figures_python/Fig3_domains_analyse.py results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/Fig3D_domain_exon_annotation_both.pickle results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/Fig3D_

python bin/creating_figures_python/Fig3_domains_analyse.py results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/Fig3D_domain_exon_annotation_core.pickle results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/Fig3D_


Rscript bin/plotting_figures_R/F3_consecutiveFraction.R results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/Fig3D_F3_exonsChangeFractionOnMerge_atit.csv results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/Fig3D_F3_exonsChangeFractionOnMerge_atit.pdf

Rscript bin/plotting_figures_R/F3_consecutiveFraction.R results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/Fig3D_F3_exonsChangeFractionOnMerge_both.csv results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/Fig3D_F3_exonsChangeFractionOnMerge_both.pdf

Rscript bin/plotting_figures_R/F3_consecutiveFraction.R results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/Fig3D_F3_exonsChangeFractionOnMerge_core.csv results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/Fig3D_F3_exonsChangeFractionOnMerge_core.pdf


Rscript bin/plotting_figures_R/F3_domExFraction.R results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/Fig3D_atit_perFractionDomains_domainsPerspective_ConsecCoding_M.csv results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/Fig3D_atit_perFractionDomains_domainsPerspective_ConsecCoding_M.pdf

Rscript bin/plotting_figures_R/F3_domExFraction.R results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/Fig3D_atit_perFractionDomains_domainsPerspective_ConsecCoding_W.csv results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/Fig3D_atit_perFractionDomains_domainsPerspective_ConsecCoding_W.pdf

Rscript bin/plotting_figures_R/F3_domExFraction.R results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/Fig3D_atit_perFractionDomains_domainsPerspective_normal_M.csv results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/Fig3D_atit_perFractionDomains_domainsPerspective_normal_M.pdf

Rscript bin/plotting_figures_R/F3_domExFraction.R results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/Fig3D_atit_perFractionDomains_domainsPerspective_normal_W.csv results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/Fig3D_atit_perFractionDomains_domainsPerspective_normal_W.pdf

Rscript bin/plotting_figures_R/F3_domExFraction.R results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/Fig3D_atit_perFractionDomains_exonsPerspective_ConsecCoding_M.csv results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/Fig3D_atit_perFractionDomains_exonsPerspective_ConsecCoding_M.pdf

Rscript bin/plotting_figures_R/F3_domExFraction.R results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/Fig3D_atit_perFractionDomains_exonsPerspective_ConsecCoding_W.csv results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/Fig3D_atit_perFractionDomains_exonsPerspective_ConsecCoding_W.pdf

Rscript bin/plotting_figures_R/F3_domExFraction.R results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/Fig3D_atit_perFractionDomains_exonsPerspective_normal_M.csv results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/Fig3D_atit_perFractionDomains_exonsPerspective_normal_M.pdf

Rscript bin/plotting_figures_R/F3_domExFraction.R results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/Fig3D_atit_perFractionDomains_exonsPerspective_normal_W.csv results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/Fig3D_atit_perFractionDomains_exonsPerspective_normal_W.pdf



Rscript bin/plotting_figures_R/F3_domExFraction.R results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/Fig3D_core_perFractionDomains_domainsPerspective_ConsecCoding_M.csv results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/Fig3D_core_perFractionDomains_domainsPerspective_ConsecCoding_M.pdf

Rscript bin/plotting_figures_R/F3_domExFraction.R results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/Fig3D_core_perFractionDomains_domainsPerspective_ConsecCoding_W.csv results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/Fig3D_core_perFractionDomains_domainsPerspective_ConsecCoding_W.pdf

Rscript bin/plotting_figures_R/F3_domExFraction.R results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/Fig3D_core_perFractionDomains_domainsPerspective_normal_M.csv results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/Fig3D_core_perFractionDomains_domainsPerspective_normal_M.pdf

Rscript bin/plotting_figures_R/F3_domExFraction.R results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/Fig3D_core_perFractionDomains_domainsPerspective_normal_W.csv results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/Fig3D_core_perFractionDomains_domainsPerspective_normal_W.pdf

Rscript bin/plotting_figures_R/F3_domExFraction.R results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/Fig3D_core_perFractionDomains_exonsPerspective_ConsecCoding_M.csv results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/Fig3D_core_perFractionDomains_exonsPerspective_ConsecCoding_M.pdf

Rscript bin/plotting_figures_R/F3_domExFraction.R results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/Fig3D_core_perFractionDomains_exonsPerspective_ConsecCoding_W.csv results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/Fig3D_core_perFractionDomains_exonsPerspective_ConsecCoding_W.pdf

Rscript bin/plotting_figures_R/F3_domExFraction.R results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/Fig3D_core_perFractionDomains_exonsPerspective_normal_M.csv results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/Fig3D_core_perFractionDomains_exonsPerspective_normal_M.pdf

Rscript bin/plotting_figures_R/F3_domExFraction.R results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/Fig3D_core_perFractionDomains_exonsPerspective_normal_W.csv results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/Fig3D_core_perFractionDomains_exonsPerspective_normal_W.pdf





Rscript bin/plotting_figures_R/F3_domExFraction.R results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/Fig3D_both_perFractionDomains_domainsPerspective_ConsecCoding_M.csv results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/Fig3D_both_perFractionDomains_domainsPerspective_ConsecCoding_M.pdf

Rscript bin/plotting_figures_R/F3_domExFraction.R results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/Fig3D_both_perFractionDomains_domainsPerspective_ConsecCoding_W.csv results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/Fig3D_both_perFractionDomains_domainsPerspective_ConsecCoding_W.pdf

Rscript bin/plotting_figures_R/F3_domExFraction.R results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/Fig3D_both_perFractionDomains_domainsPerspective_normal_M.csv results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/Fig3D_both_perFractionDomains_domainsPerspective_normal_M.pdf

Rscript bin/plotting_figures_R/F3_domExFraction.R results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/Fig3D_both_perFractionDomains_domainsPerspective_normal_W.csv results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/Fig3D_both_perFractionDomains_domainsPerspective_normal_W.pdf

Rscript bin/plotting_figures_R/F3_domExFraction.R results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/Fig3D_both_perFractionDomains_exonsPerspective_ConsecCoding_M.csv results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/Fig3D_both_perFractionDomains_exonsPerspective_ConsecCoding_M.pdf

Rscript bin/plotting_figures_R/F3_domExFraction.R results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/Fig3D_both_perFractionDomains_exonsPerspective_ConsecCoding_W.csv results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/Fig3D_both_perFractionDomains_exonsPerspective_ConsecCoding_W.pdf

Rscript bin/plotting_figures_R/F3_domExFraction.R results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/Fig3D_both_perFractionDomains_exonsPerspective_normal_M.csv results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/Fig3D_both_perFractionDomains_exonsPerspective_normal_M.pdf

Rscript bin/plotting_figures_R/F3_domExFraction.R results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/Fig3D_both_perFractionDomains_exonsPerspective_normal_W.csv results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/Fig3D_both_perFractionDomains_exonsPerspective_normal_W.pdf




# jan102022
python bin/creating_figures_python/Fig3_domains_section_withIntSp.py output/derived_data/results/objectsave_9606_0.6.pick results/Condition_3000Pilength_4Isf_4ExCount/condition_genes.pick  results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/Fig3D_ pfam



~/miniconda3/envs/py27/bin/python bin/creating_figures_python/Fig3_domains_analyse_withIntSp_csvwriter.py results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/Fig3D_domain_exon_annotation_super.pickle output/derived_data/results/objectsave_9606_0.6.pick results/Condition_3000Pilength_4Isf_4ExCount_onlyBetweenTfirstGexonsAndInside/Fig3D_

<!-- 625, 'containedDomain')
(3355, 'splitDomain')
(295, 'purecontain')
(3025, 'pureSplit') -->