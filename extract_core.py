#!/usr/local/bin/python3

''' Script to extract sequences from core genome determined by roary and the protein clustered file
from CD-HIT'''

from Bio import SeqIO
import argparse
import os


parser = argparse.ArgumentParser(
    description="This script uses the clustered_proteins output from roary to extract the core\
    genome.",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
)

# Adding arguments
parser.add_argument('-c', '--clusters',
                    help="input clustered proteins file.", required=True)
parser.add_argument('-n', '--strains',
                    help="Number of strains analysed in the study", type=int, required=True)
parser.add_argument('-d', '--directory',
                    help="Name of the directory containing the ffn files", required=True, type=str)
parser.add_argument('-o', '--output',
                    help="Name of the directory containing the output files", required=True, type=str)


args = parser.parse_args()  # Storing arguments in args variable.


class Cluster:

    def __init__(self, filename):
        self.filename = filename
        pass

    def readClusters(self):
        names = []
        ORFS = []

        with open(self.filename, 'r') as file:
            while True:
                hits = file.readline().rstrip().split()
                # Final of the file.
                if len(hits) == 0:
                    break
                # [:-1] to remove : in the end of orf name.
                names.append(hits[0][:-1])
                ORFS.append(hits[1:])
        return dict(zip(names, ORFS))

    def get_core_unique(self, dictionary):
        # filtering the keys with exact number of strains used in the analysis.
        # key_value[1] refers to the values in the .items() module.
        # The common way to remove duplicates is to convert the lists in sets. The block is doing all at the same time.

        # Removing duplicates
        # unique sequences can be obtained based on prokka id (MHLJKLJO_02609).
        # if there is a repetition of prefix of 8 letters, there are two proteins from the same genome.

        dict_no_repeats = {}

        for k, v in dictionary.items():
            prefix = []
            for i in v:
                prefix.append(i[:8])

            if len(prefix) == len(set(prefix)):  # discarding key with prefix repetitions
                dict_no_repeats[k] = v

        # Getting the core

        core = filter(lambda key_value: len(
            key_value[1]) == args.strains, dict_no_repeats.items())

        return dict(core)


if __name__ == '__main__':
    # Generating a dictionary with the cluster name and all orfs associated with it.
    # The file used in the argument --input is stored at args.inpuy

    # Creating the object
    clusters = Cluster(args.clusters)
    dictio = clusters.readClusters()
    # It is a dictionary with core genes.
    core = clusters.get_core_unique(dictio)
    print('Clusters do core-genoma lidos com sucesso.')

    # Output folder

    os.mkdir(args.output)

    # Reading fasta files to get core sequences:
    for k, v in core.items():
        dict_temp = {}
        for filename in os.listdir(args.directory):
            with open(os.path.join(args.directory, filename), 'r') as f:
                for record in SeqIO.parse(f, 'fasta'):
                    if record.id in v:
                        # split to get the name of the strain and the sequence associated.
                        dict_temp[os.path.splitext(filename)[0].split(
                            '/')[-1]] = record.seq
        for strain, seqs in dict_temp.items():
            with open(f'{args.output}/{k}.fasta', 'a') as f:  # appending mode
                f.write(f'>{strain}\n{seqs}\n')
