import polars as pl
import argparse
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
from time import sleep
from pathlib import Path

# Parse
parser = argparse.ArgumentParser()
parser.add_argument('--og_id')
parser.add_argument('--sp_list')
parser.add_argument('--gene_list')
parser.add_argument('--gene_protid')
parser.add_argument('--output_protein')
parser.add_argument('--output_nucleotide')
parser.add_argument('--email')
args = parser.parse_args()

# Input
og_id = args.og_id
sp_list_df = pl.read_csv(args.sp_list, separator='\t')
organism_id_list = list(sp_list_df.get_column('organism_id'))
gene_list_df = pl.read_csv(args.gene_list, separator='\t').filter(pl.col('organism_id').is_in(organism_id_list)).filter(pl.col('OG_id').is_in([og_id]))
gene_protid_df = pl.read_csv(args.gene_protid, separator='\t', has_header=False, new_columns=['gene_id', 'organism_id', 'protein_id', 'synonyms', 'uniprot_id', 'ensembl_id', 'ncbi_gid', 'description']).filter(pl.col('organism_id').is_in(organism_id_list))
gene_list_df = gene_list_df.join(sp_list_df, on='organism_id', how='left').join(gene_protid_df, on='gene_id', how='left').select(['gene_id', 'sci_name', 'protein_id', 'description'])

# Download protein and nucleotide sequences
def get_protein_and_nucleotide_sequences_from_protein_id(protein_id, email):
	# Set your email (NCBI requires this)
	Entrez.email = email
	
	# Fetch the protein record using the protein ID
	handle = Entrez.efetch(db='protein', id=protein_id, rettype='gb', retmode='text')
	protein_record = SeqIO.read(handle, 'genbank')
	handle.close()
	
	# Extract the nucleotide accession number and the coordinates for the CDS
	for feature in protein_record.features:
		if feature.type == 'CDS':
			cds_feature = feature
			break
	else:
		raise ValueError('No CDS feature found in protein record')
	
	nucleotide_acc = cds_feature.qualifiers['coded_by'][0].split(':')[0]
	nucleotide_coords = cds_feature.qualifiers['coded_by'][0].split(':')[1]
	start = nucleotide_coords.split('..')[0]
	end = nucleotide_coords.split('..')[1]
	if start[0] == '<':
		start = start[1:]
	start = int(start)
	if end[0] == '>':
		end = end[1:]
	end = int(end)
	
	# Fetch the nucleotide record using the nucleotide accession number
	handle = Entrez.efetch(db="nucleotide", id=nucleotide_acc, rettype="gb", retmode="text")
	nucleotide_record = SeqIO.read(handle, "genbank")
	handle.close()
	
	# Extract the protein and nucleotide sequences
	protein_sequence = protein_record.seq
	nucleotide_sequence = nucleotide_record.seq[start-1:end]
	return protein_sequence, nucleotide_sequence

prot_rec_list = []
nuc_rec_list = []
for row in gene_list_df.iter_rows(named=True):
	prot_seq, nuc_seq = get_protein_and_nucleotide_sequences_from_protein_id(row['protein_id'], args.email)
	sleep(1)
	seq_id = row['sci_name'] + '_' + row['protein_id']
	prot_rec_list.append(SeqRecord(prot_seq, id=seq_id, description=row['description']))
	nuc_rec_list.append(SeqRecord(nuc_seq, id=seq_id, description=row['description']))

# Save
output_dir = Path(args.output_protein).parent
if not output_dir.exists():
	output_dir.mkdir()
SeqIO.write(prot_rec_list, args.output_protein, 'fasta')
SeqIO.write(nuc_rec_list, args.output_nucleotide, 'fasta')