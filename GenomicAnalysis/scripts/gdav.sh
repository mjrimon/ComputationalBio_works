#!/bin/bash
# https://devhints.io/bash Bash scripting cheatsheet -> replace first match (1 ->2 ) remove prefix, suffix ...
# Tons of stackoverflow like this one https://stackoverflow.com/questions/47302898/redirect-stdout-and-stderr-to-file-permanently-but-keep-printing-them
# script to get the results from final_project

# for conda activate to work
# https://stackoverflow.com/questions/59665231/conda-4-7-how-to-activate-env-in-bash-script
# we have to run this script in interactive mode 'bash -i' or source ~/.bash_profile

source ~/scripts/init_gdav $1

scrpath=~/scripts

timeout=10
continuemsg="Press return or wait for timeout to continue..."

echo -e "\n\n"
echo "running $scrpath/gdav_1.sh"
echo -e "\n"

source $scrpath/gdav_1.sh "NO_INIT"

echo -e "\n"
read -p "$continuemsg" -s -t $timeout


echo -e "\n\n"
echo "running $scrpath/gdav_2.sh"
echo -e "\n"

source $scrpath/gdav_2.sh "NO_INIT"

echo -e "\n"
read -p "$continuemsg" -s -t $timeout


echo -e "\n\n"
echo "running $scrpath/gdav_3.sh"
echo -e "\n"

source $scrpath/gdav_3.sh "NO_INIT"


echo -e "\n"
read -p "$continuemsg" -s -t $timeout


echo -e "\n\n"
echo "running $scrpath/gdav_4_1.sh"
echo -e "\n"

source $scrpath/gdav_4_1.sh "NO_INIT"

echo -e "\n"
read -p "$continuemsg" -s -t $timeout


echo -e "\n\n"
echo "running $scrpath/gdav_4_2.sh"
echo -e "\n"

source $scrpath/gdav_4_2.sh "NO_INIT"

echo -e "\n"
read -p "$continuemsg" -s -t $timeout


echo -e "\n\n"
echo "running $scrpath/gdav_5.sh"
echo -e "\n"

source $scrpath/gdav_5.sh "NO_INIT"

echo -e "\n"
read -p "$continuemsg" -s -t $timeout


echo -e "\n\n"
echo "running $scrpath/DEAanalysis5.r"
echo -e "\n"

$scrpath/DEAanalysis5.r

echo -e "\n"
read -p "$continuemsg" -s -t $timeout


echo -e "\n\n"
echo "running $scrpath/gdav_6.sh"
echo -e "\n"

source $scrpath/gdav_6.sh "NO_INIT"

echo -e "\n"
read -p "$continuemsg" -s -t $timeout


echo -e "\n\n"
echo "running $scrpath/gdav_7.sh"
echo -e "\n"

source $scrpath/gdav_7.sh "NO_INIT"

echo -e "\n\n"
echo "Finished all scripts"
echo -e "\n\n"



