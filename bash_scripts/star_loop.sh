#!/bin/bash

for name in AHQF1A_YB_18 AHQF1A_YB_19 AHQF1A_YB_3 AHQNK_YB_12 AHQNK_YB_9 AHQNK_topM_2 AHQNK_topM_8 AHQN_topM_6 AHQT2_YB_4 AHQT_YB_3 AHQT_topM_4 AHQT_topM_5 AHQT_topM_6
do
~/STAR-2.7.9a/bin/Linux_x86_64_static/STAR \
--runThreadN 6 \
--runMode alignReads \
--outFilterMultimapNmax 1 \
--outFileNamePrefix ./${name} \
--genomeDir ./ \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript ID  \
--sjdbGTFtagExonParentGene Parent \
--sjdbGTFfile ./Mguttatus_256_v2.0.gene_exons.gff3 \
--quantMode GeneCounts \
--readFilesIn ./paired/${name}.1.paired.fastq ./paired/${name}.2.paired.fastq
done
