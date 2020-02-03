#!/bin/bash

split_input=$1
split_output1=$2
split_output2=$3

if [ ! -s $split_input ]; then
	echo "split_alis.bash: First input file not exists or is empty" >&2
fi

cut -f 1 $split_input > $split_output1
cut -f 2 $split_input > $split_output2
