from aiohttp import web
from server import app

web.run_app(app, host="localhost", port=8080)