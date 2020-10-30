from aiohttp import web
from elasticsearch import AsyncElasticsearch

client = AsyncElasticsearch()
routes = web.RouteTableDef()

SEARCH_AS_YOU_TYPE_FIELDS = ['SRAStudy', 'SRR_accession', 'SRX_accession', 'SRS_accession', 'SRP_accession',
                             'Sample_name', 'Library_name', 'ScientificName']


def get_hits(search_results: dict) -> list:
    return dict.get(search_results, 'hits', {}).get('hits', [])


def get_data(hits: list) -> list:
    return list(map(lambda hit: dict.get(hit, '_source', {}), hits))


def get_hit_count(search_results: dict) -> int:
    return dict.get(search_results, 'hits', {}).get('total', {}).get('value', 0)


@routes.get('/')
async def index(request):
    return web.FileResponse('./dist/index.html')


@routes.get('/simple_query_search/{search_string}')
async def search(request):
    """Seach using Elasticsearch simple_query_string API"""
    search_response = await client.search(
        {"query": {"simple_query_string": {"query": f"{request.match_info['search_string']}",
                                           "analyze_wildcard": "true",
                                           "fields": SEARCH_AS_YOU_TYPE_FIELDS,
                                           }
                   },
         # "_source": SEARCH_AS_YOU_TYPE_FIELDS,
         }
    )
    search_hits = get_hit_count(search_response)
    search_results = get_data(get_hits(search_response))
    return web.json_response({'hits': search_hits, 'rows': search_results})


@routes.get('/fuzzy_search/{search_string}')
async def search(request):
    """Seach using Elasticsearch simple_query_string API"""
    search_response = await client.search(
        {"query": {
            "multi_match": {
                "query": f"{request.match_info['search_string']}",
                "fields": SEARCH_AS_YOU_TYPE_FIELDS,
                "fuzziness": "AUTO"
            }
        }}
    )
    search_hits = get_hit_count(search_response)
    search_results = get_data(get_hits(search_response))
    return web.json_response({'hits': search_hits, 'rows': search_results})


def extract_relevant_terms(search_response, search_string):
    # Must Remove duplicates
    # Must be case insensitive
    relevant_terms = set()
    for item in get_hits(search_response):
        for field in item['_source'].values():
            if search_string.lower() in field.lower():
                relevant_terms.add(field)
    return list(relevant_terms)


@routes.get('/search_as_you_type/{search_string}')
async def search(request):
    search_string = request.match_info['search_string']
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

web.run_app(app, host="localhost", port=8080)
