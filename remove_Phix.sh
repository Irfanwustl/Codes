bsub -P Costas_Koumenis_P01_Nektaria_2018_10 -n 8 -q normal -o phix51.output -e phix51.error bowtie2 --un FGC1472_s_1.noPhiX.fastq -x /project/ibilab/library/bowtie2/phix/phix -U ./FGC1472_s_1.fastq.gz -S FGC1472_s_1.phix.sam
bsub -P Costas_Koumenis_P01_Nektaria_2018_10 -n 8 -q normal -o phix52.output -e phix52.error bowtie2 --un FGC1472_s_2.noPhiX.fastq -x /project/ibilab/library/bowtie2/phix/phix -U ./FGC1472_s_2.fastq.gz -S FGC1472_s_2.phix.sam
bsub -P Costas_Koumenis_P01_Nektaria_2018_10 -n 8 -q normal -o phix53.output -e phix53.error bowtie2 --un FGC1472_s_3.noPhiX.fastq -x /project/ibilab/library/bowtie2/phix/phix -U ./FGC1472_s_3.fastq.gz -S FGC1472_s_3.phix.sam
bsub -P Costas_Koumenis_P01_Nektaria_2018_10 -n 8 -q normal -o phix54.output -e phix54.error bowtie2 --un FGC1472_s_4.noPhiX.fastq -x /project/ibilab/library/bowtie2/phix/phix -U ./FGC1472_s_4.fastq.gz -S FGC1472_s_4.phix.sam

bsub -P Costas_Koumenis_P01_Nektaria_2018_10 -n 8 -q normal -o phix3.output -e phix3.error bowtie2 --un FGC1377_s_3.noPhiX.fastq -x /project/ibilab/library/bowtie2/phix/phix -U ./FGC1377_s_3.fastq.gz -S FGC1377_s_3.phix.sam
bsub -P Costas_Koumenis_P01_Nektaria_2018_10 -n 8 -q normal -o phix3.output -e phix3.error bowtie2 --un FGC1338_s_3.noPhiX.fastq -x /project/ibilab/library/bowtie2/phix/phix -U ./FGC1338_s_3.fastq.gz -S FGC1338_s_3.phix.sam

cut -f 4,5,13,14,15 B2/AAA-StudyInfo.txt |sort -u |grep -v NULL|awk '{if($5 ne "")print $3"_s_"$4}'
cut -f 4,5,13,14,15 B2/AAA-StudyInfo.txt |sort -u |grep -v NULL|awk '{if($5 ne "")print $3"_s_"$4}'|sort -u > sample_list_B2_2.txt
cut -f 4,5,13,14,15 B2/AAA-StudyInfo.txt |sort -u |grep -v NULL|awk '{if($5 ne "")print $2"\t"$5"\t"$1}'|sort -u > sample_list_B2_1.txt
