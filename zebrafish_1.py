###!/usr/bin/python
###This script is used to find read sequences that do not match upstream or downstream
###Specie:danRer10

import pandas as pd
import numpy as np
import time
import argparse

parser = argparse.ArgumentParser()
parser.description='Useage: python zebrafish_1.py -sp zebrafish -sa SRR9854035 -wd /disk/yt/zebrafish/05blastn/'
parser.add_argument("-sp", "--specie", dest='sp', help="Enter specie like zebrafish", type=str, default="zebrafish")
parser.add_argument("-sa", "--sample", dest='sa', help="Enter sample like SRR9854015", type=str, default="SRR9854015")
parser.add_argument("-wd", "--workdir", dest='wd', help="Enter workdir like /disk/yt/zebrafish/05blastn/", type=str, default="/disk/yt/zebrafish/05blastn/")
args = parser.parse_args()
specie = args.sp
sample = args.sa
path = args.wd
start_time = time.time()

path_name = str(path + sample)
sra_seq = {}
sra_file_in = open(path_name + '.txt', 'r').readlines()
for i in range(0, len(sra_file_in)):
    id = sra_file_in[i].strip()
    seq = sra_file_in[i + 1].strip()
    sra_seq[id] = seq

blastn_txt = pd.read_csv(path_name + '.txt', delimiter='\t', header=None, names=['qseqid', 'sseqid', 'qstart', 'qend', 'sstart', 'send', 'qseq', 'sseq'])
qseqid = blastn_txt['qseqid']
sstart = blastn_txt['sstart']
send = blastn_txt['send']
sstart_v = list(map(float, sstart))
send_v = list(map(float, send))
es_value = np.array(send)-np.array(sstart)
blastn_txt['es'] = es_value
sra_l = []
for key in sra_seq.keys():
    if key == id:
        sra = sra_seq[sra]
        sra_l.append(sra)
blastn_txt['sra'] = sra_l
blastn_txt.to_csv(path_name + '_0.csv', index=False, encoding='utf-8') #"es" means send-sstart, "_0" means the data filtered es
print('Adding es and sra use %.2f seconds' % (time.time() - start_time))

blastn_file = pd.read_csv(path_name + '_0.csv')
blastn_data = np.array(blastn_file)
filter_seq = []
for item in blastn_data:
    es_value = item[8]
    if es_value > 0:
        filter_seq.append(item)
header = ['qseqid', 'sseqid', 'qstart', 'qend', 'sstart', 'send', 'qseq', 'sseq', 'es', 'sra_seq']
file_csv = pd.DataFrame(columns=header, data=filter_seq)
file_csv.to_csv(path_name + '.es_filter.csv', mode='w', index=False, encoding='utf-8')
print('Filtering es use %.2f seconds' % (time.time() - start_time))

full_u = []
full_d = []
UD_s = []
blastn_csv = pd.read_csv(path_name + '_1.csv', encoding='utf-8')
qstart = blastn_csv['qstart']
qend = blastn_csv['qend']
sra = blastn_csv['sra_seq']
qstart_v = np.array(list(map(int, qstart)))
qend_v = np.array(list(map(int, qend)))
sra_v = np.array(list(map(str,sra)))
for i in range(0, len(qstart_v)):
    us = str(sra_v[i])[:int(qstart_v[i] - 1)] #us=up_seq
    ds = str(sra_v[i])[int(qend_v[i]):]
    full_u.append(us)
    full_d.append(ds)
for us, ds in zip(full_u, full_d):
    us_ds = us, ds
    UD_s.append(us_ds)

new_header = ['Up', 'Down']
blastn_result = pd.DataFrame(columns=new_header, data=UD_s)
blastn_result.to_csv(path_name + '_2.csv', mode='w', index=False, encoding='utf-8')

UD_csv = pd.read_csv(path_name + '_2.csv', encoding='utf-8')
UD_csv_reset = UD_csv.reset_index()
blastn_add_full_seq_data = pd.read_csv(path_name + '._1.csv', encoding='utf-8')
UD_csv2 = pd.concat([blastn_add_full_seq_data, UD_csv_reset], axis=1, ignore_index=True)
UD_csv2.to_csv(path_name + '_3.csv', index=False, encoding='utf-8', sep=',')


UD_csv2 = pd.read_csv(path_name + '_3.csv', encoding='utf-8', on_bad_lines='skip')
with open(path_name + '_SraUp.fa', 'a+', encoding='utf-8') as Up:
    for line in UD_csv2.values:
        Up.write(('>' + str(line[1]) + str(line[2]) + '\t' + 'length=150' + '\n' + str(line[12]) + '\n'))
with open(path_name + '_SraDown.fa', 'a+', encoding='utf-8') as Down:
    for line in UD_csv2.values:
        Down.write(('>' + str(line[1]) + str(line[2]) + '\t' + 'length=150' + '\n' + str(line[13]) + '\n'))
with open(path_name + '.qsid.txt', 'a+', encoding='utf-8') as qsid:
    for line in UD_csv2.values:
        qsid.write((str(line[1]) + str(line[2]) + '\n'))
Up.close()
Down.close()
qsid.close()
print('Obtaining up/down seq use %.2f seconds' % (time.time() - start_time))

qsid_txt = pd.read_csv(path_name + '.qsid.txt',  delimiter=' ', header=None, names=['qsid'])
qsid_csv = qsid_txt.to_csv(path_name + '.qsid.csv', encoding='utf-8', index=False)
blastn1_UD = pd.read_csv(path_name + '_3.csv', encoding='utf-8')
plus_qsid = pd.read_csv(path_name + '.qsid.csv', encoding='utf-8')
blastn1_UD['qsid'] = plus_qsid['qsid']
blastn1_UD.to_csv(path_name + '_4.csv', encoding='utf-8', index=False)
