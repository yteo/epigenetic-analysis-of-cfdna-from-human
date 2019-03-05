#!/bin/bash
file=$1
outfolder=$2
### count how many repeats a read has mapped to
sort ${file}|uniq -c| sed 's/^[ \t]*//;s/[ \t]*$//' |awk -F ' ' '{print $2"\t"$1}' > ${outfolder}/allreads_repeat.txt
