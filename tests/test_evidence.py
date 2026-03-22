"""Tests for secondary evidence and confidence scoring."""

import pytest
from ec_molsubtype.evidence import compute_confidence, compute_secondary_evidence
from ec_molsubtype.models import ConfidenceLevel, EvidenceItem, MolecularSubtype, SecondaryEvidence
from ec_molsubtype.tmb import assess_tmb, compute_tmb
from ec_molsubtype.cna import assess_cna
from ec_molsubtype.signatures import assess_signatures
from tests.conftest import make_variant


class TestComputeTmb:
    def test_basic(self):
        variants = [make_variant() for _ in range(100)]
        tmb = compute_tmb(variants, panel_size_mb=1.5)
        assert tmb == pytest.approx(100 / 1.5)

    def test_excludes_silent(self):
        variants = [
            make_variant(variant_classification="Missense_Mutation"),
            make_variant(variant_classification="Silent"),
        ]
        tmb = compute_tmb(variants, panel_size_mb=1.0)
        assert tmb == 1.0

    def test_zero_panel_raises(self):
        with pytest.raises(ValueError):
            compute_tmb([], panel_size_mb=0.0)


class TestAssessTmb:
    def test_polemut_concordant(self):
        result = assess_tmb(250.0, MolecularSubtype.POLEmut)
        assert result.concordant is True

    def test_polemut_discordant(self):
        result = assess_tmb(5.0, MolecularSubtype.POLEmut)
        assert result.concordant is False

    def test_mmrd_concordant(self):
        result = assess_tmb(40.0, MolecularSubtype.MMRd)
        assert result.concordant is True

    def test_none_tmb(self):
        result = assess_tmb(None, MolecularSubtype.POLEmut)
        assert result.concordant is None


class TestAssessCna:
    def test_p53abn_high_cna(self):
        result = assess_cna(0.55, MolecularSubtype.p53abn)
        assert result.concordant is True

    def test_p53abn_low_cna(self):
        result = assess_cna(0.1, MolecularSubtype.p53abn)
        assert result.concordant is False


class TestAssessSignatures:
    def test_polemut_concordant(self):
        weights = {"SBS10a": 0.4, "SBS10b": 0.2, "SBS1": 0.05}
        result = assess_signatures(weights, MolecularSubtype.POLEmut)
        assert result.concordant is True

    def test_polemut_no_sbs10(self):
        weights = {"SBS1": 0.5, "SBS5": 0.3}
        result = assess_signatures(weights, MolecularSubtype.POLEmut)
        assert result.concordant is False

    def test_none_weights(self):
        result = assess_signatures(None, MolecularSubtype.POLEmut)
        assert result.concordant is None


class TestComputeConfidence:
    def test_all_concordant_high(self):
        evidence = SecondaryEvidence(
            tmb=EvidenceItem(value=250.0, concordant=True),
            cna_burden=EvidenceItem(value=0.05, concordant=True),
            signatures=EvidenceItem(concordant=True),
        )
        assert compute_confidence(ConfidenceLevel.HIGH, evidence) == ConfidenceLevel.HIGH

    def test_some_discordant_moderate(self):
        evidence = SecondaryEvidence(
            tmb=EvidenceItem(value=5.0, concordant=False),
            cna_burden=EvidenceItem(value=0.05, concordant=True),
            signatures=EvidenceItem(concordant=True),
        )
        assert compute_confidence(ConfidenceLevel.HIGH, evidence) == ConfidenceLevel.MODERATE

    def test_low_stays_low(self):
        evidence = SecondaryEvidence(
            tmb=EvidenceItem(value=250.0, concordant=True),
        )
        assert compute_confidence(ConfidenceLevel.LOW, evidence) == ConfidenceLevel.LOW

    def test_no_evidence_keeps_primary(self):
        evidence = SecondaryEvidence()
        assert compute_confidence(ConfidenceLevel.HIGH, evidence) == ConfidenceLevel.HIGH
