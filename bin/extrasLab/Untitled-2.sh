#!/bin/bash
## By: Shashi Bhushan Pandit
## Date: Jan 2020

init(){
  width=$(tput cols)
  while IFS= read -r line
  do 
    len=$((width/2-${#line}/2-10))
    printf "%$((len))s%s\n" "" "$line"
  done << EOF

*************************************************************
****        Welcome to directory tree creation           ****
****        Choose purpose of directory                  ****
****          1. project                                 ****
****          2. manuscript                              ****
****          3. workspace (1st time)                    ****
*************************************************************

EOF
}

mkworksp() {
  if [[ ${#dir} -eq 0 ]]; then
   tput bold ;echo "\nNo directory name set. Give correct inputs\n"; tput sgr0
   return
  fi

  tput bold
  echo "Preparing workspace directory tree in $dir .......\n"; 
  read -p "Press [Y/N] to continue:  " opt

  if [[ $opt == 'Y' ]] || [[ $opt == 'y' ]]; then
  	echo "Creating directories .....\n"
    for j in "${list[@]}"; do
 	 echo "$dir/$j"; 
     mkdir -p $dir/$j
    done
  else 
    echo 'No directories were created'
  fi
  tput sgr0
  return 
}

mktree(){

  if [[ ${#projName} -eq 0 ]] || [[ ${#dir} -eq 0 ]]; then
    tput bold
	echo "\n The project/manuscript name or directory name not provided. Give correct inputs\n"
    tput sgr0
    return
  fi
  
  tput bold
  echo "Preparing directory tree for $projName in $dir .......\n"; 
  read -p "Press [Y/N] to continue:  " opt

  dir="$dir/$projName/"
  if [[ $opt == 'Y' ]] || [[ $opt == 'y' ]]; then
  	echo "Creating directories .....\n"
   for j in "${list[@]}"; do
 	echo "$dir/$j"
    mkdir -p $dir/$j
   done
  else 
   echo 'No directories were created'
  fi

  tput sgr0
  return 

}

while [[ $choice != 0 ]] ; do
   
   clear
   init
   cat << _EOF_
Please select appropriate option to make directory:
1. Project
2. Manuscript
3. Workspace (Only once when you login and start using a system)
0. Quit
_EOF_

  read -p "Select option [0-3] > " choice
  tput cup 16 0
  tput ed

  case $choice in
	1)  list=(bin dataset rawdata plots scratch results literature)
        read -p "Input project name :  " projName
  		read -p "Input working path (default is present directory): " dir
    	if [[ ${#dir} == 0 ]]; then
    		dir=`pwd`
		fi
        mktree
		;;
	2)  list=(finaldata finalplots rawfigures maintext)
     	read -p "Input manuscript name :  " projName
  		read -p "Input working path (default is present directory): " dir
    	if [[ ${#dir} == 0 ]]; then
    		dir=`pwd`
		fi
        mktree
		;;
    3) list=(bin projects manuscript data scratch personal)
  		read -p "Input workspace path (default is present directory): " dir
    	if [[ ${#dir} == 0 ]]; then
    		dir=`pwd`
		fi
        mkworksp
		;;
	0) break
		;;
	*) echo "Invalid entry"
		;;
   esac
   printf "\n\n Press any key to continue"
   read -n 1

done
clear
