"""FastAPI application for ec-molsubtype web interface."""

from __future__ import annotations

import uuid
from collections import OrderedDict
from pathlib import Path

from fastapi import FastAPI
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates

from .. import __version__
from .routes.classify import router as classify_router
from .routes.methods import router as methods_router

WEB_DIR = Path(__file__).parent
TEMPLATE_DIR = WEB_DIR / "templates"
STATIC_DIR = WEB_DIR / "static"

MAX_RESULTS = 500


class ResultStore:
    """Simple in-memory result store with FIFO eviction."""

    def __init__(self, maxsize: int = MAX_RESULTS):
        self._store: OrderedDict[str, dict] = OrderedDict()
        self._maxsize = maxsize

    def put(self, data: dict) -> str:
        rid = uuid.uuid4().hex[:12]
        if len(self._store) >= self._maxsize:
            self._store.popitem(last=False)
        self._store[rid] = data
        return rid

    def get(self, rid: str) -> dict | None:
        return self._store.get(rid)


def create_app() -> FastAPI:
    app = FastAPI(
        title="ec-molsubtype",
        description="Molecular subtyping of endometrial cancer",
        version=__version__,
    )

    app.mount("/static", StaticFiles(directory=str(STATIC_DIR)), name="static")

    templates = Jinja2Templates(directory=str(TEMPLATE_DIR))
    templates.env.globals["version"] = __version__

    app.state.templates = templates
    app.state.results = ResultStore()

    app.include_router(classify_router)
    app.include_router(methods_router)

    return app


# Module-level instance for uvicorn CLI (e.g., Cloud Run: uvicorn ec_molsubtype.web.app:app)
app = create_app()
