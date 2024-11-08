import polars as pl
import argparse
from Bio import SeqIO
from pathlib import Path
from collections import Counter

# Parse
parser = argparse.ArgumentParser()
parser.add_argument('--input')
parser.add_argument('--output')
args = parser.parse_args()

seq_id_list = []
sci_name_list = []
seq_len_list = []
DE_num_list = []
DE_ratio_list = []

for record in SeqIO.parse(args.input, 'fasta'):
	sci_name = record.id.split('_')[0] + '_' + record.id.split('_')[1]
	seq_len = len(record.seq)
	cnt = Counter(str(record.seq))
	DE_num = cnt['D'] + cnt['E']
	
	seq_id_list.append(record.id)
	sci_name_list.append(sci_name)
	seq_len_list.append(seq_len)
	DE_num_list.append(DE_num)
	DE_ratio_list.append(DE_num / seq_len)

df = pl.DataFrame({'seq_id': seq_id_list, 'sci_name': sci_name_list, 'seq_len': seq_len_list, 'DE_num': DE_num_list, 'DE_ratio': DE_ratio_list})
df.write_csv(args.output, separator='\t')
