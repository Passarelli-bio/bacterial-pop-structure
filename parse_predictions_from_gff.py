#!/usr/local/bin/python3

# ./parse_predictions_from_gff.py --help

import argparse

# List with plasmids ID generated in R.

parser = argparse.ArgumentParser(
    description="This script parses gff files based on a list of contigs",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
)

parser.add_argument('-g', '--gff_file',
                    help="gff file to be parsed", type=str, required=True)

parser.add_argument('-p', '--plasmid_file',
                    help="file containg the name of contigs", type=str, required=True)


args = parser.parse_args()


# Parsing plasmids ids:
plasmids_id = {}  # {'Strain' : ['plasmid_Contig']}

with open(args.plasmid_file, 'r') as fh:
    while True:

        line = fh.readline().strip().split()

        if len(line) == 0:
            break

        if line[1] not in plasmids_id:
            plasmids_id[line[1]] = []
            plasmids_id[line[1]].append(line[2])
        else:
            plasmids_id[line[1]].append(line[2])


# Parsing gff file
gff_info = {}  # {'Contig': {'ORF':['start', 'end']}

with open(args.gff_file, 'r') as fh:
    print(f'analysing {args.gff_file}')
    while True:
        hit = fh.readline().strip().split()

        if len(hit) == 0:
            break

        if hit[0].startswith('gnl'):

            if hit[2] == 'CDS':
                for i in hit:
                    if i.startswith('ID='):
                        temp_id = i.split(';')[0][3:]  # id with ORF name
            else:
                continue

            if hit[0] not in gff_info:
                gff_info[hit[0]] = {}

                if temp_id not in gff_info[hit[0]]:
                    gff_info[hit[0]][temp_id] = [hit[3], hit[4]]
            else:
                if temp_id not in gff_info[hit[0]]:
                    gff_info[hit[0]][temp_id] = [hit[3], hit[4]]


with open(f'{args.gff_file.split(".")[0]}.plasmid_info.txt', 'a') as w:
    for strain, contigs in plasmids_id.items():
        if args.gff_file.split(".")[0] == strain:
            for c in contigs:
                for orfs in gff_info[c].keys():
                    w.write(f'{strain}\t{c}\t{orfs}\n')
