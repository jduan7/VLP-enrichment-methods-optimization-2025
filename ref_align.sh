#!/bin/bash

# This bash script maps the reads of all samples to all reference viral genomes

SUBJECT="250617_nextseq"
ROOT_DIR="/media/lorax/users/jiayi/HVP/viruscocktail"

REF_DIR="/home/jiayi/ncbi_2.14/bin/viruscocktail"
REF_GENOME="viruscocktail_genomes.fna"

W_DIR="$ROOT_DIR/$SUBJECT"
RAW_DIR="$ROOT_DIR/$SUBJECT/reads_VC_4"

THREADS=12

#genome = sample.fna
#raw reads = sample_R1.fastq.gz; sample_R2.fastq.gz

source activate /home/jiayi/miniconda3/envs/new_env

# Call the Python script (which looks for all the fastq.gz files in the read directory to get the sample names)
output=$(python3 /media/lorax/users/jiayi/HVP/viruscocktail/240927_nextseq/readname_v2.py "$RAW_DIR") #the python script looks for target directories under the working directory given here

# Make a directory for reference alignment
#mkdir $W_DIR/ref_alignment

# Go to the sunbeam output directory. By looping through all the sample names (obtained through using the readname.py), map all reference genomes against the reads of each sample
while IFS= read -r dir; do
	cd $W_DIR/ref_alignment
	mkdir $dir
	cd $dir

	source activate /home/jiayi/miniconda3/envs/new_env
	#read quality filtering
	fastp --in1 $RAW_DIR/"$dir"_R1_001.fastq.gz --in2 $RAW_DIR/"$dir"_R2_001.fastq.gz --out1 trimmed_"$dir"_1.fastq.gz --out2 trimmed_"$dir"_2.fastq.gz -j fastp_stats.json -h fastp_stats.html -q 30 --thread $THREADS
	#change env
	source deactivate /home/jiayi/miniconda3/envs/new_env

	source activate /home/jiayi/miniconda3/envs/aln_env
	#copy the reference genomes to the working directory and index the genome
	cp $REF_DIR/$REF_GENOME ./
	bwa index $REF_GENOME
	#align reads to genome
	bwa mem -M $REF_GENOME trimmed_"$dir"_1.fastq.gz trimmed_"$dir"_2.fastq.gz > ref_"$dir".sam
	#convert SAM to BAM
	samtools view -S -b -o ref_"$dir".bam ref_"$dir".sam
	#sort BAM file
	samtools sort -o ref_sorted_"$dir".bam ref_"$dir".bam
	#index BAM file
	samtools index ref_sorted_"$dir".bam
	#count alignment depth
	#samtools mpileup -A -a -Q 0 -o pileup_"$dir".txt -f "$dir"_virus_sequences.fna sorted_"$dir".bam
	#alignment depth
	#head -n 1 pileup_"$dir".txt
	source deactivate /home/jiayi/miniconda3/envs/aln_env

	source activate /home/jiayi/miniconda3/envs/pysam-env
	#calculate mapped reads of good quality
	python $ROOT_DIR/coverage.py ref_sorted_"$dir".bam $REF_GENOME > $ROOT_DIR/$SUBJECT/"$dir"_ref.read
	source deactivate /home/jiayi/miniconda3/envs/pysam-env
	#add a column showing the raw read file that the alignment was run on
	awk -v new_col_name="raw_read" -v new_col="$dir" 'BEGIN {OFS="\t"} NR==1 {print $0, new_col_name; next} {print $0, new_col}' $ROOT_DIR/$SUBJECT/"$dir"_ref.read > temp && mv temp $ROOT_DIR/$SUBJECT/"$dir"_ref.read
done <<< "$output"

#merge the coverage tables together
awk 'NR==1{print; next} FNR > 1' $ROOT_DIR/$SUBJECT/*_ref.read > $ROOT_DIR/$SUBJECT/ref_alignment/"$SUBJECT"_ref_cov.tsv
#rm -r $ROOT_DIR/$SUBJECT/*.read
