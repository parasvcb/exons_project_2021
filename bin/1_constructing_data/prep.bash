#!/bin/bash

Updates_pending="file fetching with oragnism ID. routine keeping_proteins_genes_raw"
echo Script name: $0
echo $# arguments
if [ "$#" -ne 3 ]; then
    echo "Illegal number of parameters, please enter the 1. genome id (9606 for human) 2. outdir  and 3 if you need fresh start 0 for no and 1 for yes"
    exit
fi

organism=$1
outdir=$2
freshStart=$3

keeping_proteins_genes_raw (){
    # input is directory of copmmon files
	documentation="this routine takes arguemnet as commondata folder (fnameOut) and list of file urls, Gene_info (master), Swissprot proteins, Ensembl proteins for organism 
	9606 so far, will update it soon ->, if files are present locally, it will check the timestamp and tries to refetch them if they differs"
    fnameOut="$1"
    echo "fnameOUt $fnameOut"
    #need organism wise update and if genome updates then the name of file also, especially if the assembly is mentioned"

    file_gene_info='https://ftp.ncbi.nih.gov/gene/DATA/gene_info.gz'
    
    out_file_gene_info="$fnameOut/gene_info.gz"

    if test -e "$out_file_gene_info"
    then zflag=(-z "$out_file_gene_info")
    else zflag=()
    fi
    echo "curl -o $out_file_gene_info ${zflag[@]} $file_gene_info"
    # curl -o "$out_file_gene_info" "${zflag[@]}" "$file_gene_info"
    #-z, --time-cond <time> (HTTP FTP) Request a file that has been modified later than the given time and date, or one that has been modified before that time. The <date expression> can be all sorts of date strings or if it doesn't match any internal ones, it is taken as a filename and tries to get the modification date (mtime) from <file> instead. See the curl_getdate(3) man pages for date expression details. Start the date expression with a dash (-) to make it request for a document that is older than the given date/time, default is a document that is newer than the specified date/time. If this option is used several times, the last one will be used.
}

getGenesfromGeneInfo(){
    fname=$1 #"$outdir/common_files/gene_info.gz"
    outname=$2 #"$outdir/$organism/source_data/"
    txid=$3 #$organism
    #zcat $fname | awk '$1 == "9606" {print $2}' >$outname
    mkdir -p $outname/gene_tables/
    mkdir -p $outname/Refseq_protein
    python='/home/paras/miniconda3/envs/py27/bin/python'
    #zcat $fname | awk -F'\t' 'BEGIN {OFS = FS};  ($1 == "9606" && $10 == "protein-coding") {print $2}' > $outname/genesToRetrieve.txt
    
    echo "'$txid'"
    echo "$outname/genesToRetrieve_$txid.txt"
    #zcat $fname | awk -F'\t' 'BEGIN {OFS = FS};  ($1 == $txid && $10 == "protein-coding") {print $2}' > $outname/genesToRetrieve_$txid.txt
    zcat $fname | awk -F'\t' 'BEGIN {OFS = FS};  ($1 == '"$txid"' && $10 == "protein-coding") {print $2}' > $outname/genesToRetrieve_$txid.txt
    $python bin/1_constructing_data/downloadGenes.py $outname/genesToRetrieve_$txid.txt $outname/gene_tables/ $outname/Refseq_protein
    exit
}

getUnAssignedProteins(){
    # awk '(NR>1) {print $2}' ../9606_1/source_data/SSP/NP_000005.2.dat.ss | head | tr -d '\n'
    # takes two args, Refseq_portineDir and dir having secondary structure
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

keeping_proteins_genes_raw "$outdir/common_files"
# removed the concept of downloading swisspoprt and ensembl builds, just map what s retrieved from biomart, no fancy self mapping
echo "Files are refreshed and downloaded, if any error, kindly try changing proxy or update the listed URLS"
echo "Next gene info will be used to screen the genes assocuiated with our organism ${organism}"

getGenesfromGeneInfo "$outdir/common_files/gene_info.gz" "$outdir/$organism/source_data/" $organism

read -p "directory:$outdir/$organism/source_data/Refseq_protein contains freshly fetched proteins, enter the loaction of alreafy conputer source having name of files in format of 'NP_000005.2.dat.ss' to compare the files" ssdirPrev
#../9606_1/source_data/SSP
getUnAssignedProteins "$outdir/$organism" "$ssdirPrev"

#~/exonsdrive/project/protein_splicing/projectDir$ bash bin/1_constructing_data/prep.bash 9606 outdir/latest_5_org_2022_27_nov/ 1