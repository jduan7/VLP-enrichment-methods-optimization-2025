#!/bin/bash

#This script runs through all the samples that have viral contigs have map their respective raw reads to the contigs

ROOT_DIR="/media/lorax/users/jiayi/HVP/viruscocktail"
SUBJECT="250617_nextseq"
W_DIR="$ROOT_DIR/$SUBJECT/megahit_cenote/sunbeam_output"
RAW_DIR="$ROOT_DIR/$SUBJECT/reads"
THREADS=12

#make a directory list that contains all the samples that have viral contigs
find "$W_DIR/virus/cenote_taker" -type f -name "*_virus_sequences.fna" -printf "%h\n" | awk -F'/' '{print $NF}' | sort -u > $W_DIR/vcontig_dir.txt

#genome = sample.fna
#raw reads = sample_R1.fastq.gz; sample_R2.fastq.gz

source activate /home/jiayi/miniconda3/envs/new_env

# Go to the sunbeam output directory. By looping through all the sample names (obtained through using the readname.py), map all reference genomes against the reads of each sample
while IFS= read -r dir; do
        cd $W_DIR/virus/cenote_taker/$dir
        mkdir vcontig_alignment
        cd vcontig_alignment

        source activate /home/jiayi/miniconda3/envs/new_env
        #read quality filtering
        fastp --in1 $RAW_DIR/"$dir"_R1_001.fastq.gz --in2 $RAW_DIR/"$dir"_R2_001.fastq.gz --out1 trimmed_"$dir"_1.fastq.gz --out2 trimmed_"$dir"_2.fastq.gz -j fastp_stats.json -h fastp_stats.html -q 30 --thread $THREADS
        #change env
        source deactivate /home/jiayi/miniconda3/envs/new_env

        source activate /home/jiayi/miniconda3/envs/aln_env
        #copy the reference genomes to the working directory and index the genome
        cp $W_DIR/virus/cenote_taker/$dir/"$dir"_virus_sequences.fna ./
        bwa index "$dir"_virus_sequences.fna
        #align reads to genome
        bwa mem -M "$dir"_virus_sequences.fna trimmed_"$dir"_1.fastq.gz trimmed_"$dir"_2.fastq.gz > vcontig_"$dir".sam
        #convert SAM to BAM
        samtools view -S -b -o vcontig_"$dir".bam vcontig_"$dir".sam
        #sort BAM file
        samtools sort -o vcontig_sorted_"$dir".bam vcontig_"$dir".bam
        #index BAM file
        samtools index vcontig_sorted_"$dir".bam

        source deactivate /home/jiayi/miniconda3/envs/aln_env

        source activate /home/jiayi/miniconda3/envs/pysam-env
        #calculate mapped reads of good quality
        python /media/lorax/users/jiayi/HVP/coverage.py vcontig_sorted_"$dir".bam "$dir"_virus_sequences.fna > $W_DIR/"$dir"_vcontig.read
        source deactivate /home/jiayi/miniconda3/envs/pysam-env
        #add a column showing the raw read file that the alignment was run on
        awk -v new_col_name="raw_read" -v new_col="$dir" 'BEGIN {OFS="\t"} NR==1 {print $0, new_col_name; next} {print $0, new_col}' $W_DIR/"$dir"_vcontig.read > temp && mv temp $W_DIR/"$dir"_vcontig.read
done < $W_DIR/vcontig_dir.txt

#merge the coverage tables together
awk 'NR==1{print; next} FNR > 1' $W_DIR/*_vcontig.read > $W_DIR/"$SUBJECT"_vcontig_cov.tsv
#rm -r $W_DIR/*_vcontig.read

