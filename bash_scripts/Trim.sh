#this script is adapted from Findlay Finseth's processing script
#   originally found at: carnation.dbs.umt.edu:/home/findley_finseth/findley_carnation/Reseq/fastqs_to_bam_PE_mem.sh
# This script takes two fastq files and makes them into bams
# This script takes ZIPPED fastq files
# This script takes PAIRED END fastqs that need to be matched	

### CHECK IF ARGUMENTS INCLUDED ON COMMAND LINE

if [ -z $1 ] && [ -z $3 ] && [ -z $4 ]
then
    echo ""
    echo "Usage: fastqs_to_bam.sh /file/path/prefix /path/to/genome.fa [readGroup] [threads]"
    echo ""
    exit
fi

### CHECK IF A SAMPLE FILE WAS PROVIDED

filename=$1
if [ -z $filename ]
then
    echo "You must specify a path to the file prefix: /path/to/sample[.1.fastq.gz] (omit stuff in brackets). Exiting..."
    exit
fi

### CHECK IF A READ GROUP IS SPECIFIED. IF NOT, DEFAULT TO FILE PREFIX
rdgrp=$3
if [ -z $rdgrp ]
then
    rdgrp=$(basename $filename)
fi

### CHECK IF THE NUMBER OF THREADS IS SPECIFIED. IF NOT, DEFAULT OT 1
t=$4
if [ -z $t ]
then
    t=1
fi

echo $filename
echo $rdgrp
echo $t

 	gunzip ${filename}.1.fastq.gz
	gunzip ${filename}.2.fastq.gz
	#bunzip2 ${filename}.1.fastq.bz2
	#bunzip2 ${filename}.2.fastq.bz2

	# remove adapters and low quality seq
	echo "Trimming low quality reads and adapters from ${filename}"
	java -jar /home/thom_nelson/opt/Trimmomatic-0.35/trimmomatic-0.35.jar PE \
	     -threads $t -phred33 ${filename}.1.fastq ${filename}.2.fastq \
	     ${filename}.1.paired.fastq ${filename}.1.unpaired.fastq \
	     ${filename}.2.paired.fastq ${filename}.2.unpaired.fastq \
	     ILLUMINACLIP:/home/thom_nelson/resources/Illumina/Many.TruSeq.PE.fa:2:20:10:4 \
	     LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36  
	     
