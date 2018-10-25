#!/bin/bash
#script for tp1 of metagenomics

# - Cleaning/Trimming
#--------------------
#--fastq_filter FILENAME      filter and truncate sequences in FASTQ file
#--fastq_ascii  INT           FASTQ input quality score ASCII base char (33)

# - Merging
#----------
#Paired-end reads merging
#  --fastq_mergepairs FILENAME merge paired-end reads into one sequence

#--fastq_maxdiffs INT         maximum number of different bases in overlap (10)
#--fastq_minovlen             minimum length of overlap between reads (10)


#Cleaning/Trimming (Quality control 20) and Merging (use of Alientrimer)
touch all_merged.txt
R1_files="../fastq/*R1.fastq.gz"
for file in $R1_files
do
	echo $file
      	zcat $file > temp_R1.fastq
	R2_file=$(echo $file | sed "s:R1:R2:g")
	echo $R2_file
	zcat $R2_file > temp_R2.fastq
	#Trimming reads
	java -jar AlienTrimmer_0.4.0/src/AlienTrimmer.jar -if temp_R1.fastq -ir temp_R2.fastq -c AlienTrimmer_0.4.0/alienTrimmerPF8contaminants.fasta  -of fwdtrim_file -or revtrim_file
	#Merging
	./vsearch --fastq_mergepairs fwdtrim_file --reverse revtrim_file --fastq_minovlen 40 --fastq_maxdiffs 15 --fastqout $file."merged" --label_suffix  
	cat $file.merged >> all_merged.txt 
done


#formating the headlines


#Dereplication (fullength on both strands)



