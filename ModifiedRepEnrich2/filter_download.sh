# download Repbase_consensus.fa from /server/RepBase/protected/RepBase22.11.fasta (combination of pseudo.ref, humanrep.ref and humansub.ref)
awk '/^>/{s=++d".fasta"} {print > s}' Repbase_consensus.fa
mv *.fasta ./Fasta/
cd Fasta
for file in *.fasta
do
dos2unix $file
done

for file in *
do
if [ -f "$file" ]
then 
newname=`head -1 $file|cut -d ">" -f2|cut -d$'\t' -f1|sed 'y/[]/__/'`
mv $file "$newname".fa
fi
done
mkdir Present
for file in ../../RepEnrich2_setup_hg19/*.fa; do a=$(echo $file|cut -d "/" -f4); mv $a ./Present; done
