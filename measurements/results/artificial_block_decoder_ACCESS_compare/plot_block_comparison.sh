#!/bin/bash
#
# Script for plotting block decoding speed comparison
#

if [ "$1" == "" ]; then
    echo "No arguments provided"
    echo "Try: ./plot_block_comparison.sh <population>"
    echo "E.g: ./plot_block_comparison.sh 0.1"
    exit 1
fi

dataFiles=(*.data)
# clean up old eps and data files
if [ -e "${dataFiles[0]}" ];
then
  echo "Removing old data."
  ls | grep -P ".*[0-9].data" | xargs -d"\n" rm
fi

p=$1

echo "Generating command files for p=$p!"
escaped_pop="0\\.${p:2}"
stripped_pop="0${p:2}"

op_file="./block_access_${stripped_pop}.data"

grep " ${escaped_pop}" ./block_decoder_compare.txt > $op_file

gnuplot -e "p='${stripped_pop}'" plot_block_comparison.gnu
