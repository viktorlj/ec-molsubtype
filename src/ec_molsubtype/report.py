"""Output formatting for classification results."""

from __future__ import annotations

import json
from pathlib import Path

from .models import ClassificationResult


def result_to_json(result: ClassificationResult, indent: int = 2) -> str:
    """Serialize a ClassificationResult to JSON string."""
    return result.model_dump_json(indent=indent)


def result_to_dict(result: ClassificationResult) -> dict:
    """Convert a ClassificationResult to a dictionary (JSON-safe, enums as strings)."""
    return result.model_dump(mode="json")


def write_json_report(result: ClassificationResult, output_path: str | Path) -> None:
    """Write a single classification result to a JSON file."""
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(result_to_json(result))


def results_to_tsv(results: list[ClassificationResult]) -> str:
    """Convert a list of results to TSV format (summary table)."""
    headers = [
        "sample_id",
        "primary_subtype",
        "confidence",
        "multiple_classifier",
        "secondary_features",
        "tmb",
        "cna_burden",
        "n_flags",
    ]

    lines = ["\t".join(headers)]

    for r in results:
        tmb_val = ""
        if r.secondary_evidence.tmb and r.secondary_evidence.tmb.value is not None:
            tmb_val = f"{r.secondary_evidence.tmb.value:.1f}"

        cna_val = ""
        if r.secondary_evidence.cna_burden and r.secondary_evidence.cna_burden.value is not None:
            cna_val = f"{r.secondary_evidence.cna_burden.value:.3f}"

        row = [
            r.sample_id,
            r.primary_subtype.value,
            r.confidence.value,
            str(r.multiple_classifier.is_multiple),
            ";".join(r.multiple_classifier.secondary_features),
            tmb_val,
            cna_val,
            str(len(r.flags)),
        ]
        lines.append("\t".join(row))

    return "\n".join(lines)


def write_tsv_report(results: list[ClassificationResult], output_path: str | Path) -> None:
    """Write batch results to a TSV file."""
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(results_to_tsv(results))


def format_human_readable(result: ClassificationResult) -> str:
    """Format a classification result for human reading."""
    lines = [
        f"{'='*60}",
        f"Sample: {result.sample_id}",
        f"{'='*60}",
        f"",
        f"Primary Subtype: {result.primary_subtype.value}",
        f"Confidence: {result.confidence.value}",
        f"",
    ]

    # Classification path
    lines.append("Classification Path:")
    for step in result.classification_path:
        marker = "+" if step.result == "positive" else "-"
        line = f"  Step {step.step}: {step.test} [{marker}] {step.result}"
        if step.variant:
            line += f" ({step.variant})"
        lines.append(line)

    lines.append("")

    # Multiple classifier
    if result.multiple_classifier.is_multiple:
        mc = result.multiple_classifier
        lines.append(f"Multiple Classifier: Yes")
        lines.append(f"  Secondary features: {', '.join(mc.secondary_features)}")
        if mc.tp53_variant:
            lines.append(f"  TP53: {mc.tp53_variant}")
        if mc.mmr_evidence:
            lines.append(f"  MMR: {mc.mmr_evidence}")
        if mc.pole_variant:
            lines.append(f"  POLE: {mc.pole_variant}")
        lines.append("")

    # Clinical notes
    if result.clinical_notes:
        lines.append("Clinical Notes:")
        for note in result.clinical_notes:
            lines.append(f"  - {note}")
        lines.append("")

    # Flags
    if result.flags:
        lines.append("Flags:")
        for flag in result.flags:
            lines.append(f"  ! {flag}")
        lines.append("")

    lines.append(f"{'='*60}")
    return "\n".join(lines)
