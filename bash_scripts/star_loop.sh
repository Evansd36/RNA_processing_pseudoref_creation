#!/bin/bash

maindir=./
names=AHQF1i_LD_11

for name in ${names}

do
~/STAR-2.7.10a/bin/Linux_x86_64_static/STAR \
--runThreadN 6 \
--runMode alignReads \
--outFilterMultimapNmax 1 \
--outFileNamePrefix ./${name} \
--genomeDir ${maindir} \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript ID  \
--sjdbGTFtagExonParentGene Parent \
--sjdbGTFfile ${maindir}/Mguttatus_256_v2.0.gene.gff3 \
--quantMode GeneCounts \
--readFilesIn ${maindir}/${name}.1.paired.fastq ${maindir}/${name}.2.paired.fastq
done
