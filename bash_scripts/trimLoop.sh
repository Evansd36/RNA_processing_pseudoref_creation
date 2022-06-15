maindir=/Users/lilafishman/Desktop/cberg_YNP_rna 

num=AHQF1x_SD_33

for n in ${num}	 
do echo n; bash ${maindir}/RNA_processing/bash_scripts/Trim.sh ${maindir}/${n} ${maindir}/${n}
done