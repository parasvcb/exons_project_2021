#!/bin/bash

Updates_pending="file fetching with oragnism ID. routine keeping_proteins_genes_raw"
echo Script name: $0
echo $# arguments
if [ "$#" -ne 3 ]; then
    echo "Illegal number of parameters, please enter the 1. genome id (9606 for human) 2. outdir  and 3 if you need fresh start 0 for no and 1 for yes"
fi

organism=$1
outdir=$2
freshStart=$3

# pid_raw_multifasta = os.path.join(source_data,'%s_ref.faa' % organism)
# ens_raw_multifasta = os.path.join(source_data,'%s_ens.faa.gz' % organism)
# swiss_raw_multifasta = os.path.join(source_data,'%s_swiss.faa.gz' % organism)

# derived_data = os.path.join(outdir1,"derived_data/")
# #-> added this time
# stride_reformat_dir = derived_data + "stride_pssm_form/"
# stride_exons_dir = derived_data + "exons_wise/stride_exons/"
# stride_reformat_dirAF = derived_data + "stride_pssm_formAF/"
# stride_exons_dirAF = derived_data + "exons_wise/stride_exonsAF/"

# pfam_der_writer = derived_data + "domains/pfam/"
# cath_der_writer = derived_data + "domains/cath/"
# ss_exons_dir = derived_data + "exons_wise/ss_exons/"
# ss_exons_dir06 = derived_data + "exons_wise/ss_exons_0.6/"
# aaseq_exons_dir = derived_data + "exons_wise/aaseq_exons/"
# pfam_der_writer = derived_data + "domains/pfam/"
# cath_der_writer = derived_data + "domains/cath/"
# faaDir = os.path.join(derived_data,'structure_data/fasta_file_pdb')
# lisDir = os.path.join(derived_data,'structure_data/res_no_pdb')
# faaDirAF = os.path.join(derived_data,'structure_data/fasta_file_pdbAF')
# lisDirAF = os.path.join(derived_data,'structure_data/res_no_pdbAF')

mktree(){
    # create directory structure in output dir (follows convention of naming on the taxid)
  tput bold
    dir="$outdir/$organism/"
    if [ "$(find $dir -type f)" ]; then
        echo "contains files"
        read -p "directory:$dir contains files, Press [Y to delete and fresh proceed] to continue:  " opt
    else
        opt='y'
        echo "empty and proceeding"
    fi
    #if any file in the stipulate directiory name, ask user for input to delete, else proceed default
    
    if [[ $opt == 'Y' ]] || [[ $opt == 'y' ]]; then
        rm -r $dir
        echo "Deleting pre-existing data and proceeding ahead .....\n"
        listSource=(source_data/pfam source_data/SSP source_data/Refseq_protein source_data/ensemble_data source_data/gene_tables source_data/swissprot_data) 
        listDerived=(derived_data/pickles derived_data/results derived_data/structure_data/PDB_structure_derived_data derived_data/structure_data/PDB_structure_derived_dataAF derived_data/stride_pssm_form/ derived_data/exons_wise/stride_exons/ derived_data/stride_pssm_formAF/ derived_data/exons_wise/stride_exonsAF/)

        listcommon=(common_files scratch)
        for j in "${listSource[@]}"; do
            echo "$dir/$j"
            mkdir -p $dir/$j
        done
        for j in "${listDerived[@]}"; do
            echo "$dir/$j"
            mkdir -p $dir/$j
        done
        for j in "${listcommon[@]}"; do
            echo "$outdir/$j"
            mkdir -p $outdir/$j
        done
    else 
    echo 'No directories were created'
    fi
    tput sgr0
    return 
}

keeping_proteins_genes_raw (){
	documentation="this routine takes arguemnet as commondata folder (fnameOut) and list of file urls, Gene_info (master), Swissprot proteins, Ensembl proteins for organism 
	9606 so far, will update it soon ->, if files are present locally, it will check the timestamp and tries to refetch them if they differs"
    fnameOut="$1"
    echo "fnameOUt $fnameOut"
    #need organism wise update and if genome updates then the name of file also, especially if the assembly is mentioned"


    file_gene_info='https://ftp.ncbi.nih.gov/gene/DATA/gene_info.gz'
    
    declare -A animals=( [9606]="UP000005640" [6239]="UP000001940" [7227]="UP000000803" [7955]="UP000000437" [10090]="UP000000589" )
    declare -A ensemblProtLinks=( [9606]="http://ftp.ensembl.org/pub/current_fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz" [6239]="http://ftp.ensembl.org/pub/current_fasta/caenorhabditis_elegans/pep/Caenorhabditis_elegans.WBcel235.pep.all.fa.gz" [7227]="http://ftp.ensembl.org/pub/current_fasta/drosophila_melanogaster/pep/Drosophila_melanogaster.BDGP6.32.pep.all.fa.gz" [7955]="http://ftp.ensembl.org/pub/current_fasta/danio_rerio/pep/Danio_rerio.GRCz11.pep.all.fa.gz" [10090]="http://ftp.ensembl.org/pub/current_fasta/mus_musculus/pep/Mus_musculus.GRCm39.pep.all.fa.gz" )

    # bash>=4, help fetch organism specific database files
    # file_uni1='https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606.fasta.gz'
    # file_uni2='https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606_additional.dat.gz'
    # file_ens='http://ftp.ensembl.org/pub/current_fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz'

    #above listed files may not wor with the mouse, fly and elegans species.
    #https://www.uniprot.org/proteomes?facets=proteome_type:1&query=(organism_id:10090)
    #UP000000589 mouse, 10090
    #UP000000803, fly, 7227
    #UP000000437, fish, 7955
    #UP000001940, worm, 6239

    file_uni1="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/${animals[$id]}/${animals[$id]}_${id}.fasta.gz"
    file_uni2="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/${animals[$id]}/${animals[$id]}_${id}_additional.fasta.gz"
    file_ens="${ensemblProtLinks[$id]}"
    
    out_file_gene_info="$fnameOut/gene_info.gz"
    out_file_uni1="$fnameOut/UP_${organism}_fasta.gz"
    out_file_uni2="$fnameOut/UP_${organism}_adfasta.gz"
    out_file_ens="$fnameOut/EN_${organism}_fasta.gz"
    
    if test -e "$out_file_gene_info"
    then zflag=(-z "$out_file_gene_info")
    else zflag=()
    fi
    echo "curl -o $out_file_gene_info ${zflag[@]} $file_gene_info"
    #curl -o "$out_file_gene_info" "${zflag[@]}" "$file_gene_info"
    #-z, --time-cond <time> (HTTP FTP) Request a file that has been modified later than the given time and date, or one that has been modified before that time. The <date expression> can be all sorts of date strings or if it doesn't match any internal ones, it is taken as a filename and tries to get the modification date (mtime) from <file> instead. See the curl_getdate(3) man pages for date expression details. Start the date expression with a dash (-) to make it request for a document that is older than the given date/time, default is a document that is newer than the specified date/time. If this option is used several times, the last one will be used.

    if test -e "$out_file_uni1"
    then zflag=(-z "$out_file_uni1")
    else zflag=()
    fi
    echo "curl -o $out_file_uni1 ${zflag[@]} $file_uni1"
    #curl -o "$out_file_uni1" "${zflag[@]}" "$file_uni1"

    if test -e "$out_file_uni2"
    then zflag=(-z "$out_file_uni2")
    else zflag=()
    fi
    echo "curl -o $out_file_uni2 ${zflag[@]} $file_uni2"
    #curl -o "$out_file_uni2" "${zflag[@]}" "$file_uni2"

    if test -e "$out_file_ens"
    then zflag=(-z "$out_file_ens")
    else zflag=()
    fi
    echo "curl -o $out_file_ens ${zflag[@]} $file_ens"
    #curl -o "$out_file_ens" "${zflag[@]}" "$file_ens"
}

getGenesfromGeneInfo(){
    fname=$1
    outname=$2
    txid=$3
    #zcat $fname | awk '$1 == "9606" {print $2}' >$outname
    python='/home/paras/miniconda3/envs/py27/bin/python'
    if [ -s $outname/genesToRetrieve_$txid.txt ]; then
        $python bin/constructing_data/downloadGenes.py $outname/genesToRetrieve_$txid.txt $outname/gene_tables/ $outname/Refseq_protein
    else
        #zcat $fname | awk -F'\t' 'BEGIN {OFS = FS};  ($1 == "9606" && $10 == "protein-coding") {print $2}' > $outname/genesToRetrieve.txt
        zcat $fname | awk -F'\t' 'BEGIN {OFS = FS};  ($1 == $txid && $10 == "protein-coding") {print $2}' > $outname/genesToRetrieve_$txid.txt
        $python bin/constructing_data/downloadGenes.py $outname/genesToRetrieve_$txid.txt $outname/gene_tables/ $outname/Refseq_protein
    fi
}

getUnAssignedProteins(){
    #awk '(NR>1) {print $2}' ../9606_1/source_data/SSP/NP_000005.2.dat.ss | head | tr -d '\n'

    #takes two args, Refseq_portineDir and dir having secondary structure
    outname=$1
    prevDir=$2
    ssPendingDir=$outname/tempProteins
    refseqDir="${outname}/source_data/Refseq_protein"
    ssStoreDir="${outname}/source_data/SSP"

    # listSource=(source_data/pfam source_data/SSP source_data/Refseq_protein source_data/ensemble_data source_data/gene_tables source_data/swissprot_data) 
    # listDerived=(derived_data/pickles derived_data/results derived_data/structure_data/PDB_structure_derived_data derived_data/structure_data/PDB_structure_derived_dataAF derived_data/stride_pssm_form/ derived_data/exons_wise/stride_exons/ derived_data/stride_pssm_formAF/ derived_data/exons_wise/stride_exonsAF/)

    python='/home/paras/miniconda3/envs/py27/bin/python'
    mkdir -p $outname/unassignedProteins
    mkdir -p $ssPendingDir
    #prog,refDir,prevss,curss,pendingDir,fname=sys.argv
    echo "bin/constructing_data/unassignedProteins.py $refseqDir $prevDir $ssStoreDir $ssPendingDir ${outname}/source_data/map_ssOldNewandRun"
    $python bin/constructing_data/unassignedProteins.py $refseqDir $prevDir $ssStoreDir $ssPendingDir "${outname}/source_data/map_ssOldNewandRun"
}

mktree
keeping_proteins_genes_raw "$outdir/common_files"
echo "Files are refreshed and downloaded, if any error, kindly try changing proxy or update the listed URLS"
echo "Next gene info will be used to screen the genes assocuiated with our organism ${organism}"

getGenesfromGeneInfo "$outdir/common_files/gene_info.gz" "$outdir/$organism/source_data/" $organism

read -p "directory:$outdir/$organism/source_data/Refseq_protein contains freshly fetched proteins, enter the loaction of alreafy conputer source having name of files in format of 'NP_000005.2.dat.ss' to compare the files" ssdirPrev
#../9606_1/source_data/SSP
getUnAssignedProteins "$outdir/$organism" "$ssdirPrev"

