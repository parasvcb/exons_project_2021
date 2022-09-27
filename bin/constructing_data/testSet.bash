#!/bin/bash

declare -A animals=( [9606]="UP000005640" [6239]="UP000001940" [7227]="UP000000803" [7955]="UP000000437" [10090]="UP000000589" )


file_uni1='https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606.fasta.gz'
file_uni2='https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606_additional.dat.gz'

id=9606

file_uni11="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/${animals[$id]}/${animals[$id]}_${id}.fasta.gz"
file_uni2='https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606_additional.dat.gz'

    #UP000000589 mouse, 10090
    #UP000000803, fly, 7227
    #UP000000437, fish, 7955
    #UP000001940, worm, 6239
    #UP000005640, human, 9606

echo $file_uni1
echo $file_uni11