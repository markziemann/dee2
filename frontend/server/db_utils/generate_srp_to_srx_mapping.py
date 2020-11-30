import csv

mapping = {}
with open(r'C:\Users\James\PycharmProjects\dee2\SRA_Accessions.tab', 'r', newline='', encoding='utf8') as f:
    reader = csv.DictReader(f, delimiter='\t', quoting=csv.QUOTE_NONE)
    for row in reader:
        if row['Accession'][:3] == "SRX":
            if not row['Accession'] in mapping:
                mapping.update({row['Study']: row['Accession']})

with open('SRP_to_SRX.tsv', 'w', encoding='utf8') as f:
    for key, value in mapping.items():
        f.write(f'{key}\t{value}\n')
