#!/bin/bash 

source ~/.bash_profile

# Sort barcode bams by spot ID
for f in *_bcs.bam; do echo "Sorting: ${f%%.*}"; samtools sort -n $f ${f%%.*}_nameSorted; done 
# -m 10240000000

# Sort BED files by spot ID
for f in *_output/*_Final.bed; do b=${f##*/}; echo "Sorting: ${b%%.*}"; sort -k4 -u $f > ${b%%.*}_nameSorted.bed; done

# Convert sorted barcode BAM to SAM
for f in *_nameSorted.bam; do echo "Converting: ${f%%.*}"; samtools view $f > ${f%%.*}.sam; done 

# Preform the join for each species! 
while read s; do 
	echo "Joining ${s} on mm9"; 
	LANG=en_EN join -a1 -1 4 -2 1 -o 1.4 1.1 1.2 1.3 1.6 1.7 1.8 1.9 1.10 1.11 1.12 2.10 \
		<(LANG=en_EN sort -k4 ${s}_mm9_Final_nameSorted.bed) \
		<(LANG=en_EN sort -k1 ${s}_bcs_nameSorted.sam) \
       | tr ' ' '\t' \
       | awk -F $'\t' 'BEGIN {OFS = FS} {$2="MOUSE_"$2;$9="GENE:"$9;$10="ACC:"$10;$6="MOUSE_"$6;$13="BC:"substr($12,0,12);$14="UMI:"substr($12,13,8);$15="FEAT:CODING";print}' \
       > ${s}_mm9_Joined.bed; 
done < samples.txt

while read s; do
        echo "Joining ${s} on hg19";
        LANG=en_EN join -a1 -1 4 -2 1 -o 1.4 1.1 1.2 1.3 1.6 1.7 1.8 1.9 1.10 1.11 1.12 2.10 \
                <(LANG=en_EN sort -k4 ${s}_hg19_Final_nameSorted.bed) \
                <(LANG=en_EN sort -k1 ${s}_bcs_nameSorted.sam) \
       | tr ' ' '\t' \
       | awk -F $'\t' 'BEGIN {OFS = FS} {$2="HUMAN_"$2;$9="GENE:"$9;$10="ACC:"$10;$6="HUMAN_"$6;$13="BC:"substr($12,0,12);$14="UMI:"substr($12,13,8);$15="FEAT:CODING";print}' \
       > ${s}_hg19_Joined.bed;
done < samples.txt

# combine the joins
while read samples; do echo $samples; cat ${samples}_mm9_Joined.bed ${samples}_hg19_Joined.bed > ${samples}_FinalJoined.bed; done <samples.txt


# garbage collector
rm *_bcs_nameSorted.sam *_Final_nameSorted.bed *_hg19_Joined.bed *_mm9_Joined.bed

