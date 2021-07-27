#!/bin/bash
for name in *.fastq.gz
do
	gunzip ${name
	java -jar /home/thom_nelson/opt/Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 ${name}.1.fastq
