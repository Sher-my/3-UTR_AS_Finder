###!/bin/bash
###This script is used to explore read that match multiple 3UTR
###Specie:danRer10
###SRR9854035/SRR9854037=12hpf=paired

wkd=/disk/yt/zebrafish/ ###set work dir
cd ${wkd}

mkdir 02clean_data
mkdir 04unmap
cat SRR_Acc_List.txt | while read i
do
        trim_galore --paired --quality 20 --length 20 --gzip --fastqc -o 02clean_data/ 01raw_data/${i}_1.fastq.gz 01raw_data/${i}_2.fastq.gz
        hisat2 --un-conc 04unmap/${i}_unmap.fq -x 00ref/zebrafish_CDS_hi -1 02clean_data/${i}_1_val_1.fq.gz -2 02clean_data/${i}_1_val_1.fq.gz -S ${i}_useless.sam -p 16
done

mkdir 05blastn
cat SRR_Acc_List.txt | while read i
do
        pear -f 04unmap/${i}_unmap.1.fq -r 04unmap/${i}_unmap.2.fq -o 04unmap/${i}
        blastn -db 00ref/zebrafish_3UTR_bl -query 04unmap/${i}.assembled.fastq -out 05blastn/${i}.txt -outfmt ' 6 qseqid sseqid qstart qend sstart send qseq sseq' -max_hsps 1
        cp 04unmap/${i}.assembled.fastq 05blastn/${i}.txt
done

cd 05blastn
python ./zebrafish_1.py
cd ../

cat SRR_Acc_List.txt | while read i
do
        blastn -db 00ref/zebrafish_3UTR_bl -query 05blastn/${i}_SraUp.fasta -out 05blastn/${i}_Up.txt -outfmt ' 6 qseqid sseqid qstart qend sstart send qseq sseq' -max_hsps 1
        blastn -db 00ref/zebrafish_3UTR_bl -query 05blastn/${i}_SraDown.fasta -out 05blastn/${i}_Down.txt -outfmt ' 6 qseqid sseqid qstart qend sstart send qseq sseq' -max_hsps 1
done

awk '/^>/&&NR>1{print "";}{printf "%s",/^>/?$0"\t":$0}' 00ref/danRer10_3UTR2.fa > 05blastn/zebrafish_3UTR.txt

cd 05blastn
python ./zebrafish_2.py
cd ../

cat samtools.txt | while read i
do
        bash ${i}
done
###The final three files contain the same read that matches at least three locations. 
