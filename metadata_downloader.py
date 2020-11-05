import asyncio

from aiohttp.client import ClientSession
from lxml import etree

URL = "http://dee2.io/metadata/"


async def downloader(session, url_dest):
    while True:
        try:
            url, dest = next(url_dest)
        except StopIteration:
            return
        print(f"Downloading: {url}")
        rsp = await session.get(url)
        if not rsp.status == 200:
            continue  # < - Or raise error
        with open(dest, 'wb') as f:
            print(f"Writing: {dest}")
            f.write(await rsp.read())


async def run(url, destination):
    async with ClientSession() as session:
        files = await generate_download_urls(url, session)
        print(f"Files: {files}")
        dests = (f'{destination}/{file}' for file in files)
        urls = (url + file for file in files)
        url_dest = zip(urls, dests)
        await asyncio.gather(*[downloader(session, url_dest) for _ in range(1)])
        print("Completed!")


loop = asyncio.get_event_loop()


async def generate_download_urls(url, session):
    response = await session.get(url)
    tree = etree.fromstring(await response.read(), etree.HTMLParser())
    return [element.text for element in tree.xpath("/html/body/table/tr[position()>4]/td[2]/a")]


loop.run_until_complete(run(URL, "./metadata"))
