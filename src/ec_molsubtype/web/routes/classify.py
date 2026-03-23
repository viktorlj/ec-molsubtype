"""Classification routes."""

from __future__ import annotations

import json
from pathlib import Path

from fastapi import APIRouter, File, Form, Request, UploadFile
from fastapi.responses import HTMLResponse, StreamingResponse

from ...classify import classify_sample
from ...io import load_maf_variants_from_content, parse_maf_content
from ...models import SampleInput, SampleMetadata
from ...report import result_to_dict
from ..pdf import generate_pdf

router = APIRouter()

DEMO_DIR = Path(__file__).resolve().parent.parent.parent.parent.parent / "demo"


@router.get("/", response_class=HTMLResponse)
async def index(request: Request):
    templates = request.app.state.templates
    return templates.TemplateResponse(request, "index.html")


@router.post("/classify", response_class=HTMLResponse)
async def classify(
    request: Request,
    maf_file: UploadFile = File(...),
    metadata_file: UploadFile | None = File(None),
    sample_id: str = Form(""),
    tmb: str = Form(""),
    msi_pct: str = Form(""),
    fraction_genome_altered: str = Form(""),
    msi_threshold: str = Form("20.0"),
    mmr_ihc_mlh1: str = Form(""),
    mmr_ihc_msh2: str = Form(""),
    mmr_ihc_msh6: str = Form(""),
    mmr_ihc_pms2: str = Form(""),
    p53_ihc: str = Form(""),
):
    templates = request.app.state.templates
    results_store = request.app.state.results

    try:
        maf_content = (await maf_file.read()).decode("utf-8")
        variants = load_maf_variants_from_content(maf_content)

        # Build metadata from JSON sidecar or form fields
        if metadata_file and metadata_file.filename:
            meta_content = (await metadata_file.read()).decode("utf-8")
            meta_data = json.loads(meta_content)
            # Form IHC fields override JSON if set (form is more recent)
            for field, val in [("mmr_ihc_mlh1", mmr_ihc_mlh1), ("mmr_ihc_msh2", mmr_ihc_msh2),
                               ("mmr_ihc_msh6", mmr_ihc_msh6), ("mmr_ihc_pms2", mmr_ihc_pms2),
                               ("p53_ihc", p53_ihc)]:
                if val.strip():
                    meta_data[field] = val.strip()
            metadata = SampleMetadata(**meta_data)
        else:
            sid = sample_id.strip() or _infer_sample_id(maf_content, maf_file.filename)
            metadata = SampleMetadata(
                sample_id=sid,
                tmb=_parse_float(tmb),
                msi_pct=_parse_float(msi_pct),
                fraction_genome_altered=_parse_float(fraction_genome_altered),
                mmr_ihc_mlh1=mmr_ihc_mlh1.strip() or None,
                mmr_ihc_msh2=mmr_ihc_msh2.strip() or None,
                mmr_ihc_msh6=mmr_ihc_msh6.strip() or None,
                mmr_ihc_pms2=mmr_ihc_pms2.strip() or None,
                p53_ihc=p53_ihc.strip() or None,
            )

        sample = SampleInput(metadata=metadata, variants=variants)
        threshold = float(msi_threshold) if msi_threshold else 20.0
        result = classify_sample(sample, msi_threshold=threshold)

        result_data = result_to_dict(result)
        rid = results_store.put(result_data)

        return templates.TemplateResponse(request, "results.html", {
            "r": result_data,
            "result_id": rid,
        })

    except Exception as e:
        return templates.TemplateResponse(request, "error.html", {
            "error": str(e),
        })


@router.get("/demo/{sample_name}", response_class=HTMLResponse)
async def demo(request: Request, sample_name: str):
    templates = request.app.state.templates
    results_store = request.app.state.results

    maf_path = DEMO_DIR / f"{sample_name}.maf"
    json_path = DEMO_DIR / f"{sample_name}.json"

    if not maf_path.exists():
        return templates.TemplateResponse(request, "error.html", {
            "error": f"Demo sample '{sample_name}' not found.",
        })

    try:
        variants = load_maf_variants_from_content(maf_path.read_text())

        if json_path.exists():
            metadata = SampleMetadata(**json.loads(json_path.read_text()))
        else:
            metadata = SampleMetadata(sample_id=sample_name)

        sample = SampleInput(metadata=metadata, variants=variants)
        result = classify_sample(sample, msi_threshold=20.0)
        result_data = result_to_dict(result)
        rid = results_store.put(result_data)

        return templates.TemplateResponse(request, "results.html", {
            "r": result_data,
            "result_id": rid,
            "is_demo": True,
        })
    except Exception as e:
        return templates.TemplateResponse(request, "error.html", {
            "error": str(e),
        })


@router.get("/results/{result_id}/pdf")
async def download_pdf(request: Request, result_id: str):
    result_data = request.app.state.results.get(result_id)
    if result_data is None:
        return HTMLResponse("Result not found or expired.", status_code=404)

    pdf_bytes = generate_pdf(result_data)
    sid = result_data.get("sample_id", "sample")
    return StreamingResponse(
        iter([pdf_bytes]),
        media_type="application/pdf",
        headers={"Content-Disposition": f'attachment; filename="{sid}_subtype_report.pdf"'},
    )


@router.get("/results/{result_id}/json")
async def result_json(request: Request, result_id: str):
    result_data = request.app.state.results.get(result_id)
    if result_data is None:
        return HTMLResponse("Result not found or expired.", status_code=404)
    return result_data


def _parse_float(val: str) -> float | None:
    if not val or not val.strip():
        return None
    try:
        return float(val.strip())
    except ValueError:
        return None


def _infer_sample_id(maf_content: str, filename: str | None) -> str:
    try:
        df = parse_maf_content(maf_content)
        if "Tumor_Sample_Barcode" in df.columns:
            return str(df["Tumor_Sample_Barcode"][0])
    except Exception:
        pass
    if filename:
        return Path(filename).stem
    return "unknown"
