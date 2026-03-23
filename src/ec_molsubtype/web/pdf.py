"""PDF report generation using fpdf2."""

from __future__ import annotations

import io
from datetime import datetime

from fpdf import FPDF

def _sanitize(text: str) -> str:
    """Replace Unicode characters that built-in PDF fonts can't render."""
    return (text
            .replace("\u2014", "--")   # em dash
            .replace("\u2013", "-")    # en dash
            .replace("\u2018", "'")    # left single quote
            .replace("\u2019", "'")    # right single quote
            .replace("\u201c", '"')    # left double quote
            .replace("\u201d", '"')    # right double quote
            .replace("\u2265", ">=")   # >=
            .replace("\u2264", "<=")   # <=
            .replace("\u2192", "->")   # arrow
            .replace("\u2022", "*")    # bullet
            .replace("\u00b1", "+/-")  # plus-minus
            )


SUBTYPE_COLORS = {
    "POLEmut": (46, 125, 50),
    "MMRd": (21, 101, 192),
    "p53abn": (198, 40, 40),
    "NSMP": (97, 97, 97),
}


def generate_pdf(result: dict) -> bytes:
    """Generate a PDF classification report from a result dict."""
    pdf = FPDF()
    pdf.set_auto_page_break(auto=True, margin=20)
    pdf.add_page()
    pw = pdf.w - pdf.l_margin - pdf.r_margin  # printable width

    # --- Header ---
    pdf.set_font("Helvetica", "B", 16)
    pdf.cell(pw, 10, "Endometrial Cancer Molecular Subtype Report", new_x="LMARGIN", new_y="NEXT")
    pdf.set_font("Helvetica", "", 9)
    pdf.set_text_color(100, 100, 100)
    pdf.cell(pw, 5, f"Generated {datetime.now().strftime('%Y-%m-%d %H:%M')}  |  ec-molsubtype v{result.get('version', '0.1.0')}", new_x="LMARGIN", new_y="NEXT")
    pdf.ln(4)

    # --- Primary result ---
    subtype = result.get("primary_subtype", "")
    confidence = result.get("confidence", "")
    sample_id = result.get("sample_id", "")
    figo = result.get("figo_molecular_annotation", "")

    pdf.set_font("Helvetica", "", 11)
    pdf.set_text_color(0, 0, 0)
    pdf.cell(pw, 7, f"Sample: {sample_id}", new_x="LMARGIN", new_y="NEXT")

    r, g, b = SUBTYPE_COLORS.get(subtype, (0, 0, 0))
    pdf.set_font("Helvetica", "B", 20)
    pdf.set_text_color(r, g, b)
    pdf.cell(0, 12, subtype, new_x="END")
    pdf.set_font("Helvetica", "", 12)
    pdf.set_text_color(80, 80, 80)
    pdf.cell(0, 12, f"   Confidence: {confidence}", new_x="LMARGIN", new_y="NEXT")
    if figo:
        pdf.set_font("Helvetica", "", 10)
        pdf.cell(pw, 6, f"FIGO 2023 annotation: {figo}", new_x="LMARGIN", new_y="NEXT")
    pdf.ln(3)

    # --- Classification path ---
    _section_header(pdf, pw, "Classification Path")
    path = result.get("classification_path", [])
    for step in path:
        marker = "+" if step["result"] == "positive" else "-"
        line = f"  [{marker}] Step {step['step']}: {step['test']} = {step['result']}"
        if step.get("variant"):
            line += f"  ({_sanitize(step['variant'])})"
        pdf.set_font("Courier", "", 9)
        pdf.set_text_color(46, 125, 50) if step["result"] == "positive" else pdf.set_text_color(150, 150, 150)
        pdf.cell(pw, 5, line, new_x="LMARGIN", new_y="NEXT")
    pdf.set_text_color(0, 0, 0)
    pdf.ln(2)

    # --- Multiple classifier ---
    mc = result.get("multiple_classifier", {})
    if mc.get("is_multiple"):
        _section_header(pdf, pw, "Multiple Classifier")
        pdf.set_font("Helvetica", "", 9)
        feats = ", ".join(mc.get("secondary_features", []))
        pdf.cell(pw, 5, f"Secondary features detected: {feats}", new_x="LMARGIN", new_y="NEXT")
        if mc.get("tp53_variant"):
            pdf.cell(pw, 5, f"  TP53: {mc['tp53_variant']}", new_x="LMARGIN", new_y="NEXT")
        if mc.get("mmr_evidence"):
            pdf.cell(pw, 5, f"  MMR: {mc['mmr_evidence']}", new_x="LMARGIN", new_y="NEXT")
        if mc.get("pole_variant"):
            pdf.cell(pw, 5, f"  POLE: {mc['pole_variant']}", new_x="LMARGIN", new_y="NEXT")
        pdf.ln(2)

    # --- Secondary evidence ---
    sec = result.get("secondary_evidence", {})
    if sec:
        _section_header(pdf, pw, "Secondary Evidence")
        _evidence_row(pdf, pw, "TMB", sec.get("tmb"))
        _evidence_row(pdf, pw, "CNA burden", sec.get("cna_burden"))
        _evidence_row(pdf, pw, "Signatures", sec.get("signatures"))
        _evidence_row(pdf, pw, "Substitution spectrum", sec.get("substitution_profile"))
        pdf.ln(2)

    # --- Clinical notes ---
    notes = result.get("clinical_notes", [])
    if notes:
        _section_header(pdf, pw, "Clinical Notes")
        pdf.set_font("Helvetica", "", 9)
        for note in notes:
            pdf.multi_cell(pw, 5, _sanitize(f"  - {note}"))
        pdf.ln(2)

    # --- Flags ---
    flags = result.get("flags", [])
    if flags:
        _section_header(pdf, pw, "Flags")
        pdf.set_font("Helvetica", "", 9)
        pdf.set_text_color(180, 80, 0)
        for flag in flags:
            pdf.multi_cell(pw, 5, _sanitize(f"  ! {flag}"))
        pdf.set_text_color(0, 0, 0)
        pdf.ln(2)

    # --- Disclaimer ---
    pdf.ln(5)
    pdf.set_font("Helvetica", "I", 8)
    pdf.set_text_color(120, 120, 120)
    pdf.multi_cell(pw, 4,
        "DISCLAIMER: This report is generated by ec-molsubtype for research and clinical "
        "decision support only. Final molecular classification must be confirmed by a "
        "molecular pathologist. MMR status from sequencing should be verified by IHC. "
        "POLE VUS require additional evidence review. TP53 sequencing does not replace p53 IHC "
        "(~92% concordance). See the Methods page for full interpretation guidance."
    )

    buf = io.BytesIO()
    pdf.output(buf)
    return buf.getvalue()


def _section_header(pdf: FPDF, pw: float, title: str) -> None:
    pdf.set_font("Helvetica", "B", 11)
    pdf.set_text_color(0, 0, 0)
    pdf.cell(pw, 7, title, new_x="LMARGIN", new_y="NEXT")
    pdf.set_draw_color(200, 200, 200)
    pdf.line(pdf.l_margin, pdf.get_y(), pdf.l_margin + pw, pdf.get_y())
    pdf.ln(2)


def _evidence_row(pdf: FPDF, pw: float, label: str, item: dict | None) -> None:
    if item is None:
        return
    concordant = item.get("concordant")
    details = item.get("details", {})
    desc = details.get("description", "") if isinstance(details, dict) else str(details)
    value = item.get("value")

    pdf.set_font("Helvetica", "B", 9)
    pdf.set_text_color(0, 0, 0)
    pdf.cell(35, 5, f"  {label}:", new_x="END")

    if concordant is True:
        pdf.set_text_color(46, 125, 50)
        status = "concordant"
    elif concordant is False:
        pdf.set_text_color(198, 40, 40)
        status = "discordant"
    else:
        pdf.set_text_color(150, 150, 150)
        status = "N/A"

    pdf.set_font("Helvetica", "", 9)
    val_str = f"{value:.1f}" if isinstance(value, (int, float)) and value is not None else ""
    line = f"{val_str}  [{status}]" if val_str else f"[{status}]"
    pdf.cell(30, 5, line, new_x="END")

    pdf.set_text_color(100, 100, 100)
    pdf.set_font("Helvetica", "", 8)
    remaining = pw - 65
    if desc and len(desc) < 100:
        pdf.cell(remaining, 5, _sanitize(desc[:80]), new_x="LMARGIN", new_y="NEXT")
    else:
        pdf.cell(remaining, 5, "", new_x="LMARGIN", new_y="NEXT")
