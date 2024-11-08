import polars as pl
import argparse
from Bio import SeqIO
from collections import defaultdict
from pathlib import Path

# Parse
parser = argparse.ArgumentParser()
parser.add_argument('--species')
parser.add_argument('--fasta')
parser.add_argument('--output_dir')
args = parser.parse_args()

species_list_df = pl.read_csv(args.species, separator='\t')
input_fasta = args.fasta
target_organism_id_set = set(species_list_df.get_column('organism_id'))
organism_id2record_dict = defaultdict(list)

# Select proteomes of each Actinopterygii species
for record in SeqIO.parse(input_fasta, 'fasta'):
	organism_id = record.description.split('\t')[1]
	if organism_id in target_organism_id_set:
		organism_id2record_dict[organism_id].append(record)

# Output (for BUSCO inputs)
output_dir = Path(args.output_dir)
if not output_dir.exists():
	output_dir.mkdir()
for organism_id, seq_record_list in organism_id2record_dict.items():
	SeqIO.write(seq_record_list, f'{output_dir}/{organism_id}.fa', 'fasta')