import asyncio
import io
import zipfile

import aiohttp
from aiohttp import web
from elasticsearch import AsyncElasticsearch

from config import SEARCH_AS_YOU_TYPE_FIELDS
from helpers import *

client = AsyncElasticsearch()
routes = web.RouteTableDef()


@routes.get('/')
async def index(request):
    return web.FileResponse('./dist/improved_search.html')


@routes.get('/api/download/')
async def download(request):
    resp = web.StreamResponse()
    resp.content_type = 'application/zip'
    resp.enable_chunked_encoding()
    await resp.prepare(request)
    zip_files = []

    params = request.query.items()
    if len(params) > 100:
        # TODO: provide a nice error msg
        return

    async with aiohttp.ClientSession() as session:
        for key, value in params:
            # Split up the request params into a MultiDict
            params = [('org', key)] + [('x', srr) for srr in value.split(',')]
            async with session.get('http://dee2.io/cgi-bin/request.sh', params=params) as response:
                zip_files.append((key, await response.read()))

    zip_buffer = io.BytesIO()
    with zipfile.ZipFile(zip_buffer, "a", zipfile.ZIP_DEFLATED, False) as zip_file:
        for (species, data) in zip_files:
            with zipfile.ZipFile(io.BytesIO(data)) as zf:
                for file_name in zf.namelist():
                    zip_file.writestr(f"{species}/{file_name}", zf.open(file_name).read())
                    await asyncio.sleep(.0)

    await resp.write(zip_buffer.getvalue())
    await resp.write_eof()


# /simple_query_search/?searchString=a&perPage=20&offset=0"

@routes.get('/api/simple_query_search/')
@query_params('searchMode', 'searchString', 'offset', 'perPage')
async def search(request, /, searchMode, searchString, offset=0, perPage=20):
    """Seacrh Elasticsearch"""

    if searchMode == 'Strict':
        search_response = await client.search(
            {"from": offset,
             "size": perPage,
             "query": {"simple_query_string": {"query": f"{searchString}",
                                               "analyze_wildcard": "true",
                                               "fields": SEARCH_AS_YOU_TYPE_FIELDS,
                                               }
                       },
             }
        )
        search_hits = get_hit_count(search_response)
        search_results = get_data(get_hits(search_response))
        return web.json_response({'hits': search_hits, 'rows': search_results})
    elif searchMode == 'Fuzzy':
        search_response = await client.search(
            {"from": offset,
             "size": perPage,
             "query": {
                 "multi_match": {
                     "query": f"{searchString}",
                     "fields": SEARCH_AS_YOU_TYPE_FIELDS,
                     "fuzziness": "AUTO"
                 }
             }}
        )
        search_hits = get_hit_count(search_response)
        search_results = get_data(get_hits(search_response))
        return web.json_response({'hits': search_hits, 'rows': search_results})

    else:
        return web.json_response({'Reason': 'searchMode query parameter not one of "Strict" or "Fuzzy"'}
                                 , status=400
                                 )


@routes.get('/api/search_as_you_type/{search_string}')
@url_params('search_string')
async def search_as_you_type(request, /, search_string):
    search_string = search_string
    return web.json_response({"suggestions":
        extract_relevant_terms(
            await client.search(
                {
                    "query": {
                        "multi_match": {
                            "query": search_string,
                            "type": "bool_prefix",
                            "fields": SEARCH_AS_YOU_TYPE_FIELDS
                        }

                    },
                    # Filters out non-relevant data
                    "_source": SEARCH_AS_YOU_TYPE_FIELDS,
                }
            ),
            search_string)
    })


app = web.Application()
app.add_routes(routes)
app.add_routes([web.static('/', './dist')])

# r = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi", params={"db": "sra", "id":"SRX096035"})
