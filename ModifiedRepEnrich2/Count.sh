#!/bin/bash
file=$1
## count the occurence of reads mapping to each consensus location. collapse paired end reads
cut -f1 -d "/" ${file}|sort -u|cut -f9|sort|uniq -c| sed 's/^[ \t]*//;s/[ \t]*$//' |awk -F ' ' '{print $2"\t"$1}' > ${file%_consensuscoord.bed}_consensuscoord_count.txt

## count how many L1HS is present at each consensus position
cut -f4,9 ${file}|sort -u|cut -f2|sort|uniq -c| sed 's/^[ \t]*//;s/[ \t]*$//' |awk -F ' ' '{print $2"\t"$1}' > ${file%_consensuscoord.bed}_consensuscoord_num.txt
