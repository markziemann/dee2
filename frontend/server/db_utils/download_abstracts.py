import argparse
import asyncio
import glob
import re
from functools import partial
from itertools import takewhile, zip_longest
from lxml import etree
import csv
import warnings
from aiohttp.client import ClientSession

parser = argparse.ArgumentParser(description='Download Abstracts from NCBI SRA database for Dee2 SRPs')
parser.add_argument('--mapping', help='TSV file that maps SRPs to SRXs', default='SRP_to_SRX.tsv')
parser.add_argument('--data', help='Folder of Dee2 processed SRPs')
parser.add_argument('--dest', help='Folder to save data', default='./')

args = parser.parse_args()

HOST = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

E_SEARCH = "esearch.fcgi"

E_FETCH = "efetch.fcgi"

DEFAULT_PARAMS = {"db": "sra"}

TEST_SRA_PROJECTS = ['/dee2_data/bulk/scerevisiae/SRP038992_GSE55400.zip',
                     '/dee2_data/bulk/scerevisiae/SRP004431_GSE25107.zip',
                     '/dee2_data/bulk/scerevisiae/SRP073391_GSE80357.zip',
                     '/dee2_data/bulk/scerevisiae/SRP059949_GSE70378.zip',
                     '/dee2_data/bulk/scerevisiae/SRP020556_GSE45803.zip',
                     '/dee2_data/bulk/scerevisiae/SRP068940_GSE77257.zip',
                     '/dee2_data/bulk/scerevisiae/SRP092767_GSE89601.zip',
                     '/dee2_data/bulk/scerevisiae/SRP058910_GSE69414.zip',
                     '/dee2_data/bulk/scerevisiae/SRP193133_NA.zip', '/dee2_data/bulk/scerevisiae/SRP124652_NA.zip',
                     '/dee2_data/bulk/scerevisiae/SRP202200_GSE133136.zip',
                     '/dee2_data/bulk/scerevisiae/SRP178164_NA.zip', '/dee2_data/bulk/scerevisiae/ERP002319_NA.zip',
                     '/dee2_data/bulk/scerevisiae/SRP212214_GSE133457.zip',
                     '/dee2_data/bulk/scerevisiae/SRP166836_GSE121762.zip',
                     '/dee2_data/bulk/scerevisiae/SRP047413_GSE61663.zip',
                     '/dee2_data/bulk/scerevisiae/SRP193882_GSE130332.zip',
                     '/dee2_data/bulk/scerevisiae/SRP083773_NA.zip',
                     '/dee2_data/bulk/scerevisiae/SRP059799_GSE70191.zip']

SRA_PROJECT_RE = r'(?<!^\/)(?<=\/)(\w+)_(\w+)(?=\.zip)'

if args.data:
    PATHS = glob.glob(args.data + '/**/*.zip', recursive=True)
else:
    PATHS = TEST_SRA_PROJECTS

SRP_to_SRX = {}
with open(args.mapping, 'r', encoding='utf8') as f:
    SRP_to_SRX.update({k: v for k, v in map(lambda row: row.split('\t'), f)})


def iter_sra_geo(project_paths):
    return map(lambda x: x[0], map(partial(re.findall, SRA_PROJECT_RE), project_paths))


def take(iterator, n):
    """Return n elements from iterator"""
    return list(map(lambda x: x[1], takewhile(lambda x: x[0] < n, enumerate(iterator))))

def text(obj):
    return obj.text

async def downloader(session, host, tsv_writer, sra_geo_iter):
    while True:
        sra_geo_params = take(sra_geo_iter, 200)
        if not sra_geo_params:
            # Nothing left to download
            return

        print(f"Downloading SRA chunk beginning with: {sra_geo_params[0]}")

        # Example using requests
        # r = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi", params={"db": "sra", "id":"SRX096035"})

        srxs = ','.join(map(lambda x: SRP_to_SRX[x[0]], sra_geo_params))
        rsp = await session.post(host+E_FETCH, params=dict(DEFAULT_PARAMS, **{'id': srxs}))
        if not rsp.status == 200:
            warnings.warn(f"Failed to download abstract chunk beginning with {sra_geo_params[0]}")
            continue  # < - Or raise error'
        root = etree.fromstring(await rsp.text())
        species = map(text, root.xpath('EXPERIMENT_PACKAGE/SAMPLE/SAMPLE_NAME/SCIENTIFIC_NAME'))
        abstract = map(text, root.xpath('EXPERIMENT_PACKAGE/STUDY/DESCRIPTOR/STUDY_ABSTRACT'))
        for ((sra, geo), species, abstract) in zip(sra_geo_params, species, abstract):
            tsv_writer.writerow([sra, geo, species, abstract])


async def run(host, dest, project_paths):
    async with ClientSession() as session:
        sra_geo_iter = iter_sra_geo(project_paths)
        with open(f'{dest}/abstracts.tsv', 'w', encoding='utf8') as out_file:
            header = ['SRP', "GSE", "SPECIES", "ABSTRACT"]
            tsv_writer = csv.writer(out_file, dialect='excel-tab')
            tsv_writer.writerow(header)
            await asyncio.gather(*[downloader(session, host, tsv_writer, sra_geo_iter) for _ in range(1)])
        print("Completed!")


loop = asyncio.get_event_loop()

loop.run_until_complete(run(HOST, args.dest, PATHS))
