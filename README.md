# 3-UTR_AS_Finder
Finding alternative splicing(AS) events of 3'UTR for Danio rerio during 12-24hpf.
# step1 Preprocessing
## Enter your work dir
wkd = /disk/yt/zebrafish/<br>
cd ${wkd}<br>
mkdir 00ref<br>
cd 00ref
## Download specie RNA information from [ncbi](web:https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/035/GCF_000002035.5_GRCz10/)
gunzip GCF_000002035.5_GRCz10_rna_from_genomic.fna.gz<br>
cat GCF_000002035.5_GRCz10_rna_from_genomic.fna | grep "^>" | grep "gbkey=rRNA" | awk '{print $1}'|sed 's/>//g' > rRNA_id.list<br>
seqkit grep -f rRNA_id.list GCF_000002035.5_GRCz10_rna_from_genomic.fna > zebrafish_rRNA.fa
## Download specie CDS and 3'UTR information from [ucsc](web:http://genome.ucsc.edu/cgi-bin/hgTables)
sed 's/ /+/g' danRer10_3UTR.fa >> danRer10_3UTR2.fa<br>
bowtie-build zebrafish_rRNA.fa zebrafish_rRNA<br>
hisat2-build -p 16 danRer10_CDS.fa zebrafish_CDS_hi<br>
makeblastdb -in danRer10_3UTR2.fa -dbtype nucl -out zebrafish_3UTR_bl<br>
cd ../
## Download SRA file from ncbi
mkdir 01raw_data<br>
cd 01raw_data<br>
prefetch --option-file ../SRR_Acc_List.txt<br>
cat ../SRR_Acc_List.txt | while read i<br>
do<br>
        fastq-dump --split-3 --gzip ${i}/${i}.sra<br>
done<br>
cd ../
# step2 bash zebrafish.sh
