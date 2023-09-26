# 3-UTR_AS_Finder
Finding alternative splicing(AS) events of 3'UTR for Danio rerio during 12-24hpf.
# step1 Preprocessing
## Enter your work dir
wkd = /disk/yt/zebrafish/
cd ${wkd}
mkdir 00ref
cd 00ref
## Download specie RNA information from [ncbi](web:https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/035/GCF_000002035.5_GRCz10/)
gunzip GCF_000002035.5_GRCz10_rna_from_genomic.fna.gz
cat GCF_000002035.5_GRCz10_rna_from_genomic.fna | grep "^>" | grep "gbkey=rRNA" | awk '{print $1}'|sed 's/>//g' > rRNA_id.list
seqkit grep -f rRNA_id.list GCF_000002035.5_GRCz10_rna_from_genomic.fna > zebrafish_rRNA.fa
## Download specie CDS and 3'UTR information from [ucsc](web:http://genome.ucsc.edu/cgi-bin/hgTables)
sed 's/ /+/g' danRer10_3UTR.fa >> danRer10_3UTR2.fa
bowtie-build zebrafish_rRNA.fa zebrafish_rRNA
hisat2-build -p 16 danRer10_CDS.fa zebrafish_CDS_hi
makeblastdb -in danRer10_3UTR2.fa -dbtype nucl -out zebrafish_3UTR_bl
cd ../
# step2 bash zebrafish.sh
