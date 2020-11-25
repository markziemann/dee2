# Ingest metadata into Elasticsearch

import csv
import itertools
import os
from ..config import SEARCH_AS_YOU_TYPE_FIELDS
from elasticsearch import Elasticsearch
from elasticsearch.helpers import bulk

METADATA_DIR = os.environ['DEE2_METADATA_DIR']

metadata_file_names = list(filter(lambda file: file.endswith('.cut'), os.listdir(f'{METADATA_DIR}')))


def iter_metadata_files():
    for file_path in metadata_file_names:
        print(f"Indexing {file_path}")
        with open(f'{METADATA_DIR}/{file_path}', encoding='utf8') as file:
            yield file_path, iter(csv.reader(file, dialect='excel-tab'))


client = Elasticsearch()

properties = dict(
    map(lambda field: (field, {"type": "search_as_you_type"}), SEARCH_AS_YOU_TYPE_FIELDS)
)


def species_name(metadata_path):
    assert '_' in metadata_path
    # returns all characters until an underscore
    return ''.join(itertools.takewhile(lambda character: character != '_', metadata_path))


def create_index(client, index_name):
    """Creates an index in Elasticsearch if one isn't already there."""
    return client.indices.create(
        index=index_name,
        ignore=400,
        body=
        {
            "mappings": {
                "properties": properties
            }
        }
    )


def actions():
    for file_path, tsv_file in iter_metadata_files():
        headers = next(tsv_file)
        for row in tsv_file:
            yield dict(zip(headers, row), **{'_index': species_name(file_path)})


index_results = map(lambda metadata_path: create_index(client, species_name(metadata_path)), metadata_file_names)
result = bulk(client, actions())

print(result)
