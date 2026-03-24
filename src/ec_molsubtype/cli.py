"""Typer CLI for ec-molsubtype."""

from __future__ import annotations

from pathlib import Path

import typer

from .classify import classify_sample
from .io import load_sample, load_sample_metadata_tsv, load_maf_variants
from .models import SampleInput, SampleMetadata
from .pole import check_pole_variant
from .models import Variant
from .report import (
    format_human_readable,
    write_json_report,
    write_tsv_report,
)

app = typer.Typer(
    name="ec-molsubtype",
    help="Molecular subtyping of endometrial cancer from panel sequencing data.",
)


@app.command()
def classify(
    maf_path: Path = typer.Argument(..., help="Path to MAF file"),
    metadata: Path | None = typer.Option(None, "--metadata", "-m", help="Path to JSON metadata sidecar"),
    seg: Path | None = typer.Option(None, "--seg", help="Path to CBS segmentation file (.seg) for FGA computation"),
    output: Path | None = typer.Option(None, "--output", "-o", help="Output file path (JSON)"),
    msi_threshold: float = typer.Option(20.0, "--msi-threshold", help="MSI percentage threshold for MSI-H"),
    human_readable: bool = typer.Option(False, "--human", "-H", help="Print human-readable output"),
) -> None:
    """Classify a single sample from a MAF file."""
    sample = load_sample(maf_path, metadata_path=metadata)

    # Compute FGA from SEG file if provided and not already in metadata
    if seg and sample.metadata.fraction_genome_altered is None:
        from .cna import compute_fga_from_seg

        fga = compute_fga_from_seg(seg.read_text(), sample_id=sample.metadata.sample_id)
        sample.metadata.fraction_genome_altered = fga
        typer.echo(f"FGA computed from SEG file: {fga:.3f}")

    result = classify_sample(sample, msi_threshold=msi_threshold)

    if output:
        write_json_report(result, output)
        typer.echo(f"Result written to {output}")
    elif human_readable:
        typer.echo(format_human_readable(result))
    else:
        typer.echo(result.model_dump_json(indent=2))


@app.command("classify-batch")
def classify_batch(
    input_dir: Path = typer.Option(..., "--input-dir", help="Directory containing MAF files"),
    metadata_tsv: Path | None = typer.Option(None, "--metadata-tsv", help="TSV with sample metadata"),
    output: Path = typer.Option("results.tsv", "--output", "-o", help="Output TSV file"),
    msi_threshold: float = typer.Option(20.0, "--msi-threshold", help="MSI percentage threshold"),
) -> None:
    """Classify a batch of samples."""
    maf_files = sorted(input_dir.glob("*.maf"))
    if not maf_files:
        typer.echo(f"No .maf files found in {input_dir}", err=True)
        raise typer.Exit(1)

    metadata_map: dict[str, SampleMetadata] = {}
    if metadata_tsv:
        for meta in load_sample_metadata_tsv(metadata_tsv):
            metadata_map[meta.sample_id] = meta

    results = []
    for maf_path in maf_files:
        sample_id = maf_path.stem
        variants = load_maf_variants(maf_path)

        meta = metadata_map.get(sample_id, SampleMetadata(sample_id=sample_id))
        sample = SampleInput(metadata=meta, variants=variants)
        result = classify_sample(sample, msi_threshold=msi_threshold)
        results.append(result)
        typer.echo(f"  {sample_id}: {result.primary_subtype.value} ({result.confidence.value})")

    write_tsv_report(results, output)
    typer.echo(f"\nBatch results written to {output}")


@app.command("check-pole")
def check_pole(
    variant_str: str = typer.Argument(..., help="POLE variant (e.g., 'p.P286R')"),
    tmb: float | None = typer.Option(None, "--tmb", help="TMB value (mut/Mb)"),
) -> None:
    """Inspect a POLE variant for pathogenicity."""
    v = Variant(hugo_symbol="POLE", hgvsp_short=variant_str)
    result = check_pole_variant(v)

    typer.echo(f"Variant: POLE {variant_str}")
    typer.echo(f"Tier: {result.tier.value}")
    typer.echo(f"Pathogenic: {result.is_pathogenic}")
    typer.echo(f"In exonuclease domain: {result.tier != 'not_edm'}")

    if result.needs_secondary_evidence:
        typer.echo("\n⚠ VUS in exonuclease domain — secondary evidence required:")
        typer.echo("  - TMB >100 mut/Mb")
        typer.echo("  - SBS10a/SBS10b signature")
        typer.echo("  - C>A fraction >20%")

        if tmb is not None:
            if tmb > 100:
                typer.echo(f"\n  ✓ TMB {tmb:.1f} mut/Mb supports POLEmut")
            else:
                typer.echo(f"\n  ✗ TMB {tmb:.1f} mut/Mb does NOT support POLEmut")


@app.command()
def serve(
    host: str = typer.Option("0.0.0.0", "--host", help="Bind host"),
    port: int = typer.Option(None, "--port", help="Bind port (default: $PORT or 8000)"),
    reload: bool = typer.Option(False, "--reload", help="Auto-reload on changes"),
) -> None:
    """Start the web interface."""
    import os

    import uvicorn

    from .web.app import create_app

    if port is None:
        port = int(os.environ.get("PORT", "8000"))
    web_app = create_app()
    typer.echo(f"Starting ec-molsubtype web interface at http://{host}:{port}")
    uvicorn.run(web_app, host=host, port=port)


if __name__ == "__main__":
    app()
