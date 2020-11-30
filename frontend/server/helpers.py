from aiohttp import web


def bad_request(reason):
    web.json_response({'Reason': str(reason)}, status=400)


def query_params(*parameters):
    def wrapper(route):
        def inner(request):
            try:
                return route(request, **{parameter: request.query.get(parameter, None)
                                         for parameter in parameters if parameter in request.query})

            except KeyError:
                return bad_request(f'Missing query parameters {"".join(set(parameters).difference(parameters))}')

            except Exception as e:
                return bad_request(e)

        return inner

    return wrapper


def url_params(*parameters):
    def wrapper(route):
        def inner(request):
            try:
                return route(request, **{parameter: request.match_info.get(parameter, None)
                                         for parameter in parameters if parameter in request.match_info})

            except KeyError:
                return bad_request('Malformed URL')

            except Exception as e:
                return bad_request(e)

        return inner

    return wrapper


def get_hits(search_results: dict) -> list:
    return dict.get(search_results, 'hits', {}).get('hits', [])


def extract_index_and_data(hit):
    try:
        data = hit['_source']
        data.update({"index": hit['_index']})
    except KeyError:
        return {}
    else:
        return data


def get_data(hits: list) -> list:
    return list(map(extract_index_and_data, hits))


def get_hit_count(search_results: dict) -> int:
    return dict.get(search_results, 'hits', {}).get('total', {}).get('value', 0)


def extract_relevant_terms(search_response, search_string):
    # Must Remove duplicates
    # Must be case insensitive
    relevant_terms = set()
    for item in get_hits(search_response):
        for field in item['_source'].values():
            if search_string.lower() in field.lower():
                relevant_terms.add(field)
    return list(relevant_terms)
