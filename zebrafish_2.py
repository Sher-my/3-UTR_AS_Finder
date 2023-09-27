###!/usr/bin/python
###This script is used to find where read upstream and downstream unmatched sequences to match
###Specie:danRer10

import pandas as pd
import time
import argparse

parser = argparse.ArgumentParser()
parser.description='Useage: python zebrafish_2.py -sp zebrafish -sa SRR9854035 -wd /disk/yt/zebrafish/05blastn/'
parser.add_argument("-sp", "--specie", dest='sp', help="Enter specie like zebrafish", type=str, default="zebrafish")
parser.add_argument("-sa", "--sample", dest='sa', help="Enter sample like SRR9854015", type=str, default="SRR9854015")
parser.add_argument("-wd", "--workdir", dest='wd', help="Enter workdir like /disk/yt/zebrafish/05blastn/", type=str, default="/disk/yt/zebrafish/05blastn/")
args = parser.parse_args()
specie = args.sp
sample = args.sa
path = args.wd
start_time = time.time()

path_name = path + sample
blastn_up = pd.read_csv(path_name + '_Up.txt', delimiter='\t', header=None, names=['qsid', 'sseqid2', 'qstart2', 'qend2', 'sstart2', 'send2', 'qseq2', 'sseq2'])
blastn_down = pd.read_csv(path_name + '_Down.txt', delimiter='\t', header=None, names=['qsid', 'sseqid3', 'qstart3', 'qend3', 'sstart3', 'send3', 'qseq3', 'sseq3'])
blastn_file = pd.read_csv(path_name + '._4.csv', encoding='utf-8')
df_mix1 = blastn_file.merge(blastn_up, on='qsid', how='left')
index = df_mix1['qstart2'].notnull()
df_mix1 = df_mix1[index]
df_mix2 = blastn_file.merge(blastn_down, on='qsid', how='left')
index2 = df_mix2['qstart3'].notnull()
df_mix2 = df_mix2[index2]
UTR3_txt2 = open(path + 'zebrafish_3UTR.txt', 'r').readlines()
UTR3_out2 = open(path + 'UTR3_1.txt', 'w')
for line in UTR3_txt2:
    if '>' in line:
        UTR3_out2.write(line.replace('>', '\n'))
UTR3_out2.close()
UTR3_txt3 = pd.read_csv(path + 'UTR3_1.txt', delimiter='\t', header=None, names=['UTR', 'full_seq'])
UTR3_csv = UTR3_txt3.to_csv(path + 'UTR3_1.csv', encoding='utf-8', index=False)

filter_data1 = []
filter_data2 = []
filter_sum = []
blastn_qsid = df_mix1['qsid']
blastn_qsid2 = df_mix2['qsid']
blastn_sseqid = df_mix1['2']
blastn_sseqid2 = df_mix1['sseqid2']
blastn_sseqid3 = df_mix2['2']
blastn_sseqid4 = df_mix2['sseqid3']
blastn_qstart = df_mix1['3']
blastn_qstart2 = df_mix1['qstart2']
blastn_qstart3 = df_mix2['3']
blastn_qstart4 = df_mix2['qstart3']
blastn_sstart = df_mix1['5']
blastn_sstart2 = df_mix1['sstart2']
blastn_sstart3 = df_mix2['5']
blastn_sstart4 = df_mix2['sstart3']
blastn_qend = df_mix1['4']
blastn_qend2 = df_mix1['qend2']
blastn_qend3 = df_mix2['4']
blastn_qend4 = df_mix2['qend3']
blastn_send = df_mix1['6']
blastn_send2 = df_mix1['send2']
blastn_send3 = df_mix2['6']
blastn_send4 = df_mix2['send3']

blastn_qsid_value = list(map(str, blastn_qsid))
blastn_qsid2_value = list(map(str, blastn_qsid2))
blastn_sseqid_value = list(map(str, blastn_sseqid))
blastn_sseqid2_value = list(map(str, blastn_sseqid2))
blastn_sseqid3_value = list(map(str, blastn_sseqid3))
blastn_sseqid4_value = list(map(str, blastn_sseqid4))
blastn_qstart_value = list(map(int, blastn_qstart))
blastn_qstart2_value = list(map(int, blastn_qstart2))
blastn_qstart3_value = list(map(int, blastn_qstart3))
blastn_qstart4_value = list(map(int, blastn_qstart4))
blastn_sstart_value = list(map(int, blastn_sstart))
blastn_sstart2_value = list(map(int, blastn_sstart2))
blastn_sstart3_value = list(map(int, blastn_sstart3))
blastn_sstart4_value = list(map(int, blastn_sstart4))
blastn_qend_value = list(map(int, blastn_qend))
blastn_qend2_value = list(map(int, blastn_qend2))
blastn_qend3_value = list(map(int, blastn_qend3))
blastn_qend4_value = list(map(int, blastn_qend4))
blastn_send_value = list(map(int, blastn_send))
blastn_send2_value = list(map(int, blastn_send2))
blastn_send3_value = list(map(int, blastn_send3))
blastn_send4_value = list(map(int, blastn_send4))

for i in range(0, len(blastn_sseqid)):
    try:
        if blastn_send_value[i]-blastn_sstart_value[i] > 0 and blastn_send2_value[i]-blastn_sstart2_value[i] > 0 and blastn_sstart_value[i]-blastn_send2_value[i] > 0 and blastn_sseqid_value[i] == blastn_sseqid2_value[i] and blastn_qstart_value[i] - blastn_qend2_value[i] == 1:
            filter_data1.append(str(blastn_qsid_value[i]) + '\t' + str(blastn_qstart_value[i]) + '~' + str(blastn_qend_value[i]) + '\t' + str(blastn_sstart_value[i]) + '\t' + str(blastn_send_value[i]) + '\t' + 'forward' + '\t' + '1')
            filter_data1.append(str(blastn_qsid_value[i]) + '\t' + str(blastn_qstart2_value[i]) + '~' + str(blastn_qend2_value[i]) + '\t' + str(blastn_sstart2_value[i]) + '\t' + str(blastn_send2_value[i]) + '\t' + 'forward' + '\t' + '1')
            filter_sum.append(str(blastn_sseqid_value[i]) + '\t' + str(blastn_send2_value[i]) + '\t' + str(blastn_sstart_value[i]))
        if blastn_send4_value[i] - blastn_sstart4_value[i] > 0 and blastn_send3_value[i] - blastn_sstart3_value[i] > 0 and blastn_sstart4_value[i] - blastn_send3_value[i] > 0 and blastn_sseqid3_value[i] == blastn_sseqid4_value[i] and blastn_qstart4_value[i] + blastn_qend3_value[i] - blastn_qend3_value[i] == 1:
            filter_data2.append(str(blastn_qsid2_value[i]) + '\t' + str(blastn_qstart3_value[i]) + '~' + str(blastn_qend3_value[i]) + '\t' + str(blastn_sstart3_value[i]) + '\t' + str(blastn_send3_value[i]) + '\t' + 'forward' + '\t' + '1')
            filter_data2.append(str(blastn_qsid2_value[i]) + '\t' + str(blastn_qstart4_value[i]+blastn_qend3_value[i]) + '~' + str(blastn_qend4_value[i]+blastn_qend3_value[i]) + '\t' + str(blastn_sstart4_value[i]) + '\t' + str(blastn_send4_value[i]) + '\t' + 'forward' + '\t' + '1')
            filter_sum.append(str(blastn_sseqid3_value[i]) + '\t' + str(blastn_send3_value[i]) + '\t' + str(blastn_sstart4_value[i]))
    except:
        print('err')
header1 = ['UTR', 'reference', 'sstart', 'send', 'strand', 'drection']
header2 = ['UTR', 'sstart', 'send']
file_csv1 = pd.DataFrame(data=filter_data1)
file_csv2 = pd.DataFrame(data=filter_data2)
file_csv3 = pd.DataFrame(data=filter_sum)

file_csv1.to_csv(path_name + '.U_2.csv', mode='w', index=False, encoding='utf-8')
file_csv2.to_csv(path_name + '.D_2.csv', mode='w', index=False, encoding='utf-8')
file_csv3.to_csv(path + 'samtools_0.csv', mode='w', index=False, encoding='utf-8')

new_csv1 = pd.read_csv(path_name + '.U_2.csv', delimiter='\t', header=None, names=header1, index_col=0)
new_csv2 = pd.read_csv(path_name + '.D_2.csv', delimiter='\t', header=None, names=header1, index_col=0)
new_csv3 = pd.read_csv(path + 'samtools_0.csv', delimiter='\t', header=None, names=header2, index_col=0)
UTR3_csv_read = pd.read_csv(path + 'UTR3_1.csv', encoding='utf-8')
df_mix4 = new_csv3.merge(UTR3_csv_read, on='UTR', how='left')
index4 = df_mix4['UTR'].notnull()
df_mix4 = df_mix4[index4]
UTR_data3 = df_mix4.fillna(0)
intron_UTR = UTR_data3['UTR']
SAMTOOLS = []
samtools_slurm = open(path + 'samtools.txt', 'a+')
intron_start_site = UTR_data3['sstart']
intron_send_site = UTR_data3['send']
intron_pre_end_site = UTR_data3['full_seq']
intron_start_site_v = list(map(int, intron_start_site))
intron_send_site_v = list(map(int, intron_send_site))
intron_UTR_v = list(map(str, intron_UTR))
intron_pre_end_site_v =list(map(str, intron_pre_end_site))
for i in range(0, len(intron_start_site)):
    try:
        if intron_start_site_v[i] > 200:
            SAMTOOLS.append('samtools faidx danRer10_3UTR2.fa' + ' ' + str(intron_UTR_v[i]) + ':' + str(int(intron_start_site_v[i] - 200)) + '-' + str(int(intron_start_site[i])) + ' ' + '>>' + ' ' + 'Up.txt')
        else:
            SAMTOOLS.append('samtools faidx danRer10_3UTR2.fa' + ' ' + str(intron_UTR_v[i]) + ':' + '1' + '-' + str(int(intron_start_site[i])) + ' ' + '>>' + ' ' + 'Up.txt')
        if len(intron_pre_end_site_v[i]) - intron_send_site_v[i] > 200:
            SAMTOOLS.append('samtools faidx danRer10_3UTR2.fa' + ' ' + str(intron_UTR_v[i]) + ':' + str(int(intron_send_site_v[i])) + '-' + str(int(intron_send_site[i] + 200)) + ' ' + '>>' + ' ' + 'Down.txt')
        else:
            SAMTOOLS.append('samtools faidx danRer10_3UTR2.fa' + ' ' + str(intron_UTR_v[i]) + ':' + str(int(intron_send_site_v[i])) + '-' + str(int(len(intron_pre_end_site_v[i]))) + ' ' + '>>' + ' ' + 'Down.txt')
        SAMTOOLS.append('samtools faidx danRer10_3UTR2.fa' + ' ' + str(intron_UTR_v[i]) + ':' + str(int(intron_start_site_v[i])) + '-' + str(int(intron_send_site_v[i])) + ' ' + '>>' + ' ' + 'intron.txt')
    except:
        print('err')
for line in SAMTOOLS:
    samtools_slurm.write(line + '\n')
samtools_slurm.close()
