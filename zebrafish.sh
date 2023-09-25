#!/bin/bash
###This script is used to explore read that match multiple 3UTR
###specie:danRer10
###data:SRR9854015=24hpf=single;SRR9854035/SRR9854037=12hpf=paired

wkd=/disk/yt/zebrafish/ ###set work path
cd ${wkd}

mkdir 00ref
cd 00ref
###download RNA from https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/035/GCF_000002035.5_GRCz10/
gunzip GCF_000002035.5_GRCz10_rna_from_genomic.fna.gz
cat GCF_000002035.5_GRCz10_rna_from_genomic.fna | grep "^>" | grep "gbkey=rRNA" | awk '{print $1}'|sed 's/>//g' > rRNA_id.list
seqkit grep -f rRNA_id.list GCF_000002035.5_GRCz10_rna_from_genomic.fna > rRNA.fa
###download CDS and UTR3 from http://genome.ucsc.edu/cgi-bin/hgTables
sed 's/ /+/g' danRer10_3UTR.fa >> danRer10_3UTR2.fa
bowtie-build rRNA.fa rRNA
hisat2-build -p 16 danRer10_CDS.fa danRer10_CDS_hi
makeblastdb -in danRer10_3UTR2.fa -dbtype nucl -out danRer10_3UTR_bl
cd ../

mkdir 01raw_data
cd 01raw_data
prefetch --option-file ../SRR_Acc_List.txt
cat ../SRR_Acc_List.txt | while read i
do
        fastq-dump --split-3 --gzip ${i}/${i}.sra
done
cd ../

###mv fq to 01raw_data/

echo "SRR9854015" > Single_List.txt
echo "SRR9854035" > Paired_List.txt
echo "SRR9854037" >> Paired_List.txt

mkdir 02clean_data
mkdir 03de_rRNA
mkdir 04unmap
touch 03de_rRNA/rRNA_report
cat Paired_List.txt | while read i
do
        trim_galore --paired --quality 20 --length 20 --gzip --fastqc -o 02clean_data/ 01raw_data/${i}_1.fq 01raw_data/${i}_2.fq
        echo 'rRNA_sample:'${i} >> 03de_rRNA/rRNA_report
        bowtie -q -p 16 --un 03de_rRNA/${i}_derRNA.fq -x 00ref/rRNA -1 02clean_data/${i}_1_trimmed.fq.gz -2 02clean_data/${i}_2_trimmed.fq.gz -S 03de_rRNA/aligned_reads_mates.sam >> 03de_rRNA/rRNA_report
        hisat2 --un-conc 04unmap/${i}_unmap.fq -x 00ref/danRer10_cds_hi -U 03de_rRNA/${i}_derRNA.fq -S ${i}_useless.sam -p 16 >> 04unmap/umaplog
done

cat Single_List.txt | while read i
do
        trim_galore --quality 20 --length 20 --gzip --fastqc -o 02clean_data/ 01raw_data/${i}.fq
        echo 'rRNA_sample:'${i} >> 03de_rRNA/rRNA_report
        bowtie -q -p 12 --un 03de_rRNA/${i}_derRNA.fq -x 00ref/rRNA -U 02clean_data/${i}_trimmed.fq.gz -S 03de_rRNA/aligned_reads_mates.sam >> 03de_rRNA/rRNA_report
        hisat2 --un-conc 04unmap/${i}_unmap.fq -x 00ref/danRer10_cds_hi -U 03de_rRNA/${i}_derRNA.fq -S ${i}_useless.sam -p 16 >> 04unmap/umaplog
done

mkdir 05blastn
cat SRR_Acc_List.txt | while read i
do
        blastn -db 00ref/danRer10_3UTR_bl -query 04unmap/${i}_unmap.fq -out 05blastn/${i}.txt -outfmt ' 6 qseqid sseqid qstart qend sstart send qseq sseq' -max_hsps 1
        cp 04unmap/${i}_unmap.fq mkdir 05blastn/${i}.txt
done

cd 05blastn
python ./zebrafish_1.py
cd ../

cat SRR_Acc_List.txt | while read i
do
        blastn -db 00ref/danRer10_3UTR_bl -query 05blastn/${i}_SraUp.fasta -out 05blastn/${i}_Up.txt -outfmt ' 6 qseqid sseqid qstart qend sstart send qseq sseq' -max_hsps 1
        blastn -db 00ref/danRer10_3UTR_bl -query 05blastn/${i}_SraDown.fasta -out 05blastn/${i}_Down.txt -outfmt ' 6 qseqid sseqid qstart qend sstart send qseq sseq' -max_hsps 1
done

awk '/^>/&&NR>1{print "";}{printf "%s",/^>/?$0"\t":$0}' 00ref/danRer10_3UTR2.fa > 05blastn/danRer10_3UTR.txt

cd 05blastn
python ./zebrafish_2.py
cd ../

cat samtools.txt | while read i
do
        bash ${i}
done
###The final three files contain the same read that matches at least three locations. 
