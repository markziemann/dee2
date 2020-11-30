import csv
import argparse

parser = argparse.ArgumentParser(description='Generate SRP to SRX mapping file')
parser.add_argument('--accessions', help='Path to SRA_Accessions.tab file', default='./SRA_Accessions.tab')
parser.add_argument('--dest', help='Destination for output file', default='./SRP_to_SRX.tsv')

args = parser.parse_args()

mapping = {}
with open(args.accessions, 'r', newline='', encoding='utf8') as f:
    reader = csv.DictReader(f, delimiter='\t', quoting=csv.QUOTE_NONE)
    for row in reader:
        if row['Accession'][2] == "X" and row['Status'] == 'live':
            if not row['Study'] in mapping:
                mapping.update({row['Study']: row['Accession']})

with open(args.dest, 'w', encoding='utf8') as f:
    for key, value in mapping.items():
        f.write(f'{key}\t{value}\n')
