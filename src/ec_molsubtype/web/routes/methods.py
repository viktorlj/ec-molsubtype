"""Methods page route."""

from __future__ import annotations

import json
from pathlib import Path

from fastapi import APIRouter, Request
from fastapi.responses import HTMLResponse

router = APIRouter()

DATA_DIR = Path(__file__).resolve().parent.parent.parent / "data"


@router.get("/methods", response_class=HTMLResponse)
async def methods(request: Request):
    templates = request.app.state.templates

    pole_data = json.loads((DATA_DIR / "pole_pathogenic.json").read_text())
    tp53_data = json.loads((DATA_DIR / "tp53_pathogenic.json").read_text())

    return templates.TemplateResponse(request, "methods.html", {
        "pole_tier1": pole_data["tier1"],
        "pole_tier2": pole_data["tier2"],
        "pole_edm": pole_data["exonuclease_domain"],
        "tp53_hotspots": tp53_data["hotspots"],
        "tp53_dbd": tp53_data["dna_binding_domain"],
    })
