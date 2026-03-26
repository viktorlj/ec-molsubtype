"""CLI & installation documentation page route."""

from __future__ import annotations

from fastapi import APIRouter, Request
from fastapi.responses import HTMLResponse

from ... import __version__

router = APIRouter()


@router.get("/cli", response_class=HTMLResponse)
async def cli_docs(request: Request):
    templates = request.app.state.templates
    return templates.TemplateResponse(request, "cli.html", {
        "version": __version__,
    })
