import argparse
from ete3 import Tree
from pathlib import Path
from Bio import SeqIO
import subprocess

# Parse
parser = argparse.ArgumentParser()
parser.add_argument('--species_tree')
parser.add_argument('--gene_tree')
parser.add_argument('--aln')
parser.add_argument('--output')
args = parser.parse_args()

output_dir = Path(args.output).parent.parent.parent.parent
with open(output_dir/'mapping_file.link', 'w') as f:
    for rec in SeqIO.parse(args.aln, format='fasta'):
        spl = rec.id.split('_')
        species = spl[0] + '_' + spl[1]
        f.write(f'{rec.id} {species}\n')

with open(output_dir/'families_generax.txt', 'w') as f:
    f.write('[FAMILIES]\n')
    f.write(f'- {output_dir.name}\n')
    f.write(f'starting_gene_tree = {args.gene_tree}\n')
    f.write(f'alignment = {args.aln}\n')
    f.write(f'mapping = {output_dir}/mapping_file.link\n')
    f.write('subst_model = LG+G\n')

subprocess.run([f'generax -f {output_dir}/families_generax.txt -s {args.species_tree} -r UndatedDL -p {output_dir}/generax --strategy SPR'], shell=True)