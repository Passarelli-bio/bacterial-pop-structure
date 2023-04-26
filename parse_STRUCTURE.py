#!/usr/local/bin/python3

import sys
import os
import argparse
import pandas as pd


parser = argparse.ArgumentParser(
    description="This script parses STRUCTURE output files",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
)

parser.add_argument('-d', '--directory_run',
                    help="directory with results from for a given K tested", type=str)

parser.add_argument('-f', '--file_with_best_K',
                    help="file containg the run associated with the best predicted K", type=str)

parser.add_argument('-t', '--type_of_operation',
                    help="options: k_info --> get table of K and ln P(D) from all iterations;\
                                    ancestrality_info --> get a table of ancestrality and groups predicted", required=True, type=str)


args = parser.parse_args()


class Structure:

    def __init__(self, directory):
        self.directory = directory
        self.filename = [
            filename for filename in os.listdir(args.directory_run)]

    # Dictionary to store K and all associated Estimated Ln Prob of Data.
    def get_k_dict(self):

        k_info = {}

        for filename in self.filename:
            with open(os.path.join(args.directory_run, filename), 'r') as f:
                lines = f.readlines()

                for l in lines:
                    if 'populations assumed' in l:
                        K = l.rstrip().lstrip().split()[0]
                    elif l.startswith('Estimated Ln Prob of Data'):
                        # 30th field to get specifically the number.
                        value = l.rstrip()[30:]
                if int(K) not in k_info:
                    k_info[int(K)] = []
                    k_info[int(K)].append(float(value))
                else:
                    k_info[int(K)].append(float(value))

        return k_info

    # It requires a dictionary from get_k_dict() module.
    def write_k_info_csv(self, outname='k_info_runs.txt'):
        df = pd.DataFrame(structure.get_k_dict()).sort_index(axis=1)
        df.to_csv(outname, sep='\t', index=False)

    def get_ancestrality_dic(self):

        with open(args.file_with_best_K, 'r') as f:
            lines = f.readlines()
            ancestry_clusters = []

            # getting the index where this information begins.
            start_index = lines.index('Inferred ancestry of individuals:\n')
            end_index = lines.index(
                'Estimated Allele Frequencies in each cluster\n')
            for l in lines[start_index + 2:end_index - 2]:  # skkiping unecessary lines
                ancestry_clusters.append(l.strip().split())

            # Creating the dictionary
            dictionary = {'Id_number': [], 'Strain': [],
                          'CC1': [], 'CC2': [], 'CC3': [], 'CC4': [], 'CC5': [], 'CC6': [], 'CC7': []}

            for i in ancestry_clusters:
                dictionary['Id_number'].append(i[0])
                dictionary['Strain'].append(i[1])
                dictionary['CC1'].append(float(i[4]))
                dictionary['CC2'].append(float(i[5]))
                dictionary['CC3'].append(float(i[6]))
                dictionary['CC4'].append(float(i[7]))
                dictionary['CC5'].append(float(i[8]))
                dictionary['CC6'].append(float(i[9]))
                dictionary['CC7'].append(float(i[10]))
        return dictionary

    def write_ancestrality_info_csv(self, outname='ancestrality_info_runs.txt'):
        df = pd.DataFrame(structure.get_ancestrality_dic())
        df.to_csv(outname, sep='\t', index=False)


if __name__ == '__main__':
    # Writing a table with K-associated posterior probabilities
    if args.type_of_operation == 'k_info':
        structure = Structure(args.directory_run)
        structure.write_k_info_csv('k_runs_info.txt')

    elif args.type_of_operation == 'ancestrality_info':
        structure = Structure(args.file_with_best_K)
        # structure.write_ancestrality_info_csv()
        # write_ancestrality_info_csv()
        structure.write_ancestrality_info_csv()
