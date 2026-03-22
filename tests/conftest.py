"""Shared test fixtures for ec-molsubtype tests."""

from __future__ import annotations

import pytest

from ec_molsubtype.models import (
    SampleInput,
    SampleMetadata,
    Variant,
)


# --- Variant factories ---


def make_variant(
    hugo_symbol: str = "TP53",
    hgvsp_short: str = "p.R175H",
    variant_classification: str = "Missense_Mutation",
    t_alt_count: int = 50,
    t_ref_count: int = 50,
    clinvar_classification: str = "",
    **kwargs,
) -> Variant:
    return Variant(
        hugo_symbol=hugo_symbol,
        hgvsp_short=hgvsp_short,
        variant_classification=variant_classification,
        t_alt_count=t_alt_count,
        t_ref_count=t_ref_count,
        clinvar_classification=clinvar_classification,
        **kwargs,
    )


# --- POLE variants ---


@pytest.fixture
def pole_p286r():
    return make_variant(hugo_symbol="POLE", hgvsp_short="p.P286R")


@pytest.fixture
def pole_v411l():
    return make_variant(hugo_symbol="POLE", hgvsp_short="p.V411L")


@pytest.fixture
def pole_tier2_f367s():
    return make_variant(hugo_symbol="POLE", hgvsp_short="p.F367S")


@pytest.fixture
def pole_vus():
    """A POLE variant in the exonuclease domain but not in curated lists."""
    return make_variant(hugo_symbol="POLE", hgvsp_short="p.A300T")


@pytest.fixture
def pole_outside_edm():
    """A POLE variant outside the exonuclease domain."""
    return make_variant(hugo_symbol="POLE", hgvsp_short="p.R100H")


# --- TP53 variants ---


@pytest.fixture
def tp53_r175h():
    return make_variant(hugo_symbol="TP53", hgvsp_short="p.R175H")


@pytest.fixture
def tp53_r248w():
    return make_variant(hugo_symbol="TP53", hgvsp_short="p.R248W")


@pytest.fixture
def tp53_truncating():
    return make_variant(
        hugo_symbol="TP53",
        hgvsp_short="p.R196*",
        variant_classification="Nonsense_Mutation",
    )


@pytest.fixture
def tp53_p72r():
    """Benign polymorphism."""
    return make_variant(hugo_symbol="TP53", hgvsp_short="p.P72R")


@pytest.fixture
def tp53_low_vaf():
    """Low VAF TP53 variant (subclonal)."""
    return make_variant(
        hugo_symbol="TP53",
        hgvsp_short="p.R273H",
        t_alt_count=2,
        t_ref_count=98,
    )


# --- MMR variants ---


@pytest.fixture
def msh6_frameshift():
    return make_variant(
        hugo_symbol="MSH6",
        hgvsp_short="p.F1088fs",
        variant_classification="Frame_Shift_Del",
    )


@pytest.fixture
def mlh1_nonsense():
    return make_variant(
        hugo_symbol="MLH1",
        hgvsp_short="p.R100X",
        variant_classification="Nonsense_Mutation",
    )


@pytest.fixture
def msh2_single_het():
    return make_variant(
        hugo_symbol="MSH2",
        hgvsp_short="p.A636P",
        variant_classification="Missense_Mutation",
        clinvar_classification="Pathogenic",
    )


# --- Sample fixtures ---


@pytest.fixture
def polemut_sample():
    """Complete POLEmut sample."""
    return SampleInput(
        metadata=SampleMetadata(
            sample_id="POLE_SAMPLE",
            tmb=250.0,
            msi_pct=2.0,
            fraction_genome_altered=0.05,
            signature_weights={"SBS10a": 0.45, "SBS10b": 0.22, "SBS1": 0.05},
        ),
        variants=[
            make_variant(hugo_symbol="POLE", hgvsp_short="p.P286R"),
            make_variant(hugo_symbol="TP53", hgvsp_short="p.R248W"),
        ],
    )


@pytest.fixture
def mmrd_sample():
    """Complete MMRd sample via MSI markers."""
    return SampleInput(
        metadata=SampleMetadata(
            sample_id="MMRD_SAMPLE",
            tmb=45.0,
            msi_pct=55.0,
            fraction_genome_altered=0.15,
            signature_weights={"SBS6": 0.30, "SBS15": 0.15, "SBS1": 0.10},
        ),
        variants=[
            make_variant(
                hugo_symbol="MSH6",
                hgvsp_short="p.F1088fs",
                variant_classification="Frame_Shift_Del",
            ),
        ],
    )


@pytest.fixture
def p53abn_sample():
    """Complete p53abn sample."""
    return SampleInput(
        metadata=SampleMetadata(
            sample_id="P53ABN_SAMPLE",
            tmb=5.0,
            msi_pct=0.0,
            fraction_genome_altered=0.55,
            signature_weights={"SBS1": 0.35, "SBS5": 0.30},
        ),
        variants=[
            make_variant(hugo_symbol="TP53", hgvsp_short="p.R175H"),
        ],
    )


@pytest.fixture
def nsmp_sample():
    """Complete NSMP sample."""
    return SampleInput(
        metadata=SampleMetadata(
            sample_id="NSMP_SAMPLE",
            tmb=3.0,
            msi_pct=0.0,
            fraction_genome_altered=0.08,
            signature_weights={"SBS1": 0.40, "SBS5": 0.35},
        ),
        variants=[
            make_variant(
                hugo_symbol="PTEN",
                hgvsp_short="p.R130Q",
                variant_classification="Missense_Mutation",
            ),
        ],
    )
