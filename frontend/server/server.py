import asyncio
import io
import zipfile

import aiohttp
import elasticsearch as es
from elasticsearch import AsyncElasticsearch, Elasticsearch

from config import SEARCH_AS_YOU_TYPE_FIELDS
from helpers import *

client = AsyncElasticsearch()
cat_client = es.client.CatClient(Elasticsearch())

routes = web.RouteTableDef()

RUN_INDICES = [index['index'] for index in cat_client.indices(format='json')]
PROJECT_INDICES = ['projects']

for index in PROJECT_INDICES:
    RUN_INDICES.remove(index)

RUN_INDICES = ','.join(RUN_INDICES)
PROJECT_INDICES = ','.join(PROJECT_INDICES)


def level_to_indices_and_fields(level):
    if level == 'Projects':
        return PROJECT_INDICES, ['SRP', "GSE", "SPECIES", "ABSTRACT"]
    else:
        return RUN_INDICES, SEARCH_AS_YOU_TYPE_FIELDS


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

@routes.get('/api/search/')
@query_params('level', 'mode', 'query', 'offset', 'per_page')
async def search(request, /, level, mode, query, offset=0, per_page=20):
    """Seacrh Elasticsearch"""
    indices, fields = level_to_indices_and_fields(level)
    if mode == 'Strict':
        search_response = await client.search(
            {"from": offset,
             "size": per_page,
             "query": {"simple_query_string": {"query": f"{query}",
                                               "analyze_wildcard": "true",
                                               "fields": fields,
                                               }
                       },
             }, index=indices
        )
        search_hits = get_hit_count(search_response)
        search_results = get_data(get_hits(search_response))
        return web.json_response({'hits': search_hits, 'rows': search_results})
    elif mode == 'Fuzzy':
        search_response = await client.search(
            {"from": offset,
             "size": per_page,
             "query": {
                 "multi_match": {
                     "query": f"{query}",
                     "fields": fields,
                     "fuzziness": "AUTO"
                 }
             }}, index=indices
        )
        search_hits = get_hit_count(search_response)
        search_results = get_data(get_hits(search_response))
        return web.json_response({'hits': search_hits, 'rows': search_results})

    else:
        return web.json_response({'Reason': 'searchMode query parameter not one of "Strict" or "Fuzzy"'}
                                 , status=400
                                 )


@routes.get('/api/suggestions/')
@query_params('level', 'query')
async def suggestions(request, /, level, query):
    indices, fields = level_to_indices_and_fields(level)
    return web.json_response({"suggestions":
        extract_relevant_terms(
            await client.search(
                {
                    "query": {
                        "multi_match": {
                            "query": query,
                            "type": "bool_prefix",
                            "fields": fields
                        }

                    },
                    # Filters out non-relevant data
                    "_source": fields,
                }, index=indices
            ),
            query)
    })


app = web.Application()
app.add_routes(routes)
app.add_routes([web.static('/', './dist')])
