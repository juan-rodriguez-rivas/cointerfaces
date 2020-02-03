#!/bin/bash

# TODO: Input options handling

raw_file1=$1
raw_file2=$2
just_profile1=$3
just_profile2=$4

if [ ! -s $raw_file1 ]; then
	echo "just_profile.bash: First input file not exists or is empty" >&2
fi

if [ ! -s $raw_file2 ]; then
	echo "just_profile.bash: Second input file not exists or is empty" >&2
fi

awk ' {if($1~/^>/){print} else{gsub(/[a-z\.]/, ""); print }}' $raw_file1 > $just_profile1

awk ' {if($1~/^>/){print} else{gsub(/[a-z\.]/, ""); print }}' $raw_file2 > $just_profile2
