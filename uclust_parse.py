'''
Script to parse uclust files.
The output is a text-separated file. The first field
gives you the RECORD TYPE (S, H, C, or N).

S and C represents centroids sequences (eg. if cluster has a single protein, C and S will be the same).

If you want to get a non-redundant dataset, just take all records with C. Else, take all H (hits) and
(S) centroid secundary sequences.

Passarelli
'''

#!/usr/local/bin/python3
import os
from Bio import SeqIO
import subprocess as sp
import argparse

parser = argparse.ArgumentParser(
    description="This script uses the uclust output to get single copy orthologous genes",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
)


# Adding arguments
parser.add_argument('-d', '--directory',
                    help="directory with all fasta files to be clustered in FASTA extension", required=True)
parser.add_argument('-id', '--identity',
                    help="Identity threshold to cluster proteins", type=float, default=0.4, required=True)

args = parser.parse_args()  # Storing arguments in args variable.


# Preparing files for uclust.

all_files = os.listdir(args.directory)
fasta_files = []

for file in all_files:
    if file.endswith('.fasta'):
        fasta_files.append(file)

# Moving for desired directory
os.chdir(args.directory)

# Genome counter to detect paralogs
genome_number = 0

# Dictionary to keep information about sequences

seqs_info = {}  # {'genome1': {'filename':[genes]}}
seqs_for_uclust = {}  # {'genome1_gene1: seq'}

for filename in fasta_files:
    print(f'processing file {filename}')
    genome_number += 1
    gene_number = 1

    with open(filename, 'r') as f:
        all_original_headers = []
        for record in SeqIO.parse(f, 'fasta'):
            all_original_headers.append(record.id)
            record.id = f'Gene{gene_number}'  # Changing
            seqs_for_uclust[f'genome{genome_number}_gene{gene_number}'] = str(
                record.seq)
            gene_number += 1

        seqs_info[f'genome{genome_number}'] = {filename: all_original_headers}

print('writing seqs_info.txt.uc')

# Writing information about genomes
#{'genome1': {'filename':[genes]}}
with open('seqs_info.txt.uc', 'a') as info:
    for k, v in seqs_info.items():
        for filename in v.keys():
            info.write(f'{k}\t{filename}\n')

print('writing seqs_for_uclust.fasta')

# write fasta file for uclust input
with open('seqs_for_uclust.fasta.uc', 'a') as w:
    for k, v in seqs_for_uclust.items():
        w.write(f'>{k}\n{v}\n')

# running uclust
sp.Popen('uclust --sort seqs_for_uclust.fasta.uc --output seqs_sorted.fasta.uc',
         shell=True).wait()
print('sequences sorted')


print('running uclust')
sp.Popen(
    f'uclust --input seqs_sorted.fasta.uc --uc results.uc --id {args.identity}', shell=True).wait()

# Parsing uclust output
with open('results.uc', 'r') as fh:
    dictionary = {}
    while True:
        hit = fh.readline().strip().split()
        if len(hit) == 0:
            break

        if hit[0] != 'C' and hit[0] != '#':

            if hit[1] not in dictionary:
                dictionary[hit[1]] = []
                dictionary[hit[1]].append(hit[8])
            else:
                dictionary[hit[1]].append(hit[8])

# Removing paralogs
print('removing paralogs')
dict_no_repeats = {}

for k, v in dictionary.items():
    prefix = []
    for i in v:
        prefix.append(i.split('_')[0])

    if len(prefix) == len(set(prefix)):  # discarding key with prefix repetitions
        dict_no_repeats[k] = v


# Getting core genes

core = dict(filter(lambda key_value: len(key_value[1]) == len(fasta_files), dict_no_repeats.items()))


# Creating and moving to output directory
try:
    os.mkdir('Output_uclust')
except FileExistsError:
    pass

os.chdir('Output_uclust')

# writing all files
print(core)
for k, v in core.items():
    with open(f'{k}.fasta', 'a') as w:
        for i in v:
            w.write(f'>{i.split("_")[0]}\n{seqs_for_uclust[i]}\n')

print('Analysis is DONE')
