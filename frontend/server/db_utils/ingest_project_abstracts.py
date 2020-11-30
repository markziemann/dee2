import argparse
import csv

from elasticsearch import Elasticsearch
from elasticsearch.helpers import bulk

parser = argparse.ArgumentParser(description='Ingest project abstracts into ElasticSearch')
parser.add_argument('--abstracts', help='TSV file containing abstracts', default='abstracts.tsv')

args = parser.parse_args()

client = Elasticsearch()


def create_index(client, header):
    """Creates an index in Elasticsearch if one isn't already there."""
    return client.indices.create(
        index='projects',
        ignore=400,
        body=
        {
            "mappings": {
                "properties": dict(
                    map(lambda field: (field, {"type": "search_as_you_type"}), header)
                )
            }
        }
    )


def actions(headers, tsv_file):
    for row in tsv_file:
        yield dict(zip(headers, row), **{'_index': 'projects'})


with open(args.abstracts, 'r', newline='', encoding='utf8') as f:
    tsv_reader = csv.reader(f, dialect='excel-tab')
    headers = next(tsv_reader)
    create_index(client, headers)
    result = bulk(client, actions(headers, tsv_reader))

print(result)
