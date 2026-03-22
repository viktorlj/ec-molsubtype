"""Tests for MMR deficiency assessment."""

import pytest
from ec_molsubtype.models import ConfidenceLevel
from ec_molsubtype.mmr import assess_mmr, classify_msi_status, find_mmr_variants
from tests.conftest import make_variant


class TestClassifyMsiStatus:
    def test_msi_h(self):
        assert classify_msi_status(45.0) == "MSI-H"
        assert classify_msi_status(20.0) == "MSI-H"

    def test_msi_l(self):
        assert classify_msi_status(10.0) == "MSI-L"
        assert classify_msi_status(0.5) == "MSI-L"

    def test_mss(self):
        assert classify_msi_status(0.0) == "MSS"

    def test_custom_threshold(self):
        assert classify_msi_status(15.0, threshold=10.0) == "MSI-H"
        assert classify_msi_status(15.0, threshold=20.0) == "MSI-L"


class TestAssessMmr:
    def test_msi_h_high_confidence(self):
        """MSI-H from marker panel → MMRd high confidence."""
        result = assess_mmr([], msi_pct=55.0)
        assert result.is_mmrd is True
        assert result.confidence == ConfidenceLevel.HIGH
        assert result.msi_status == "MSI-H"

    def test_msi_h_with_mmr_mutations(self, msh6_frameshift):
        """MSI-H plus supporting MMR mutations."""
        result = assess_mmr([msh6_frameshift], msi_pct=45.0)
        assert result.is_mmrd is True
        assert result.confidence == ConfidenceLevel.HIGH
        assert len(result.mmr_mutations) > 0

    def test_mss_biallelic_mmr(self, msh6_frameshift, mlh1_nonsense):
        """MSS but biallelic MMR → moderate with discordance."""
        # Two hits in different genes — not biallelic in same gene
        # Need two hits in the same gene for biallelic
        v1 = make_variant(
            hugo_symbol="MSH6",
            hgvsp_short="p.F1088fs",
            variant_classification="Frame_Shift_Del",
        )
        v2 = make_variant(
            hugo_symbol="MSH6",
            hgvsp_short="p.R200X",
            variant_classification="Nonsense_Mutation",
        )
        result = assess_mmr([v1, v2], msi_pct=5.0)
        assert result.is_mmrd is True
        assert result.confidence == ConfidenceLevel.MODERATE
        assert result.is_biallelic is True

    def test_no_msi_biallelic_mmr(self):
        """No MSI data + biallelic MMR → moderate."""
        v1 = make_variant(
            hugo_symbol="MLH1",
            hgvsp_short="p.R100X",
            variant_classification="Nonsense_Mutation",
        )
        v2 = make_variant(
            hugo_symbol="MLH1",
            hgvsp_short="p.W200fs",
            variant_classification="Frame_Shift_Del",
        )
        result = assess_mmr([v1, v2], msi_pct=None)
        assert result.is_mmrd is True
        assert result.confidence == ConfidenceLevel.MODERATE

    def test_no_msi_single_het(self, msh2_single_het):
        """No MSI data + single het → NOT MMRd, flagged."""
        result = assess_mmr([msh2_single_het], msi_pct=None)
        assert result.is_mmrd is False
        assert any("IHC recommended" in f for f in result.flags)

    def test_mss_no_mutations(self):
        """MSS + no MMR mutations → MMR proficient."""
        result = assess_mmr([], msi_pct=0.0)
        assert result.is_mmrd is False

    def test_no_data_no_mutations(self):
        """No MSI data, no MMR mutations → proficient."""
        result = assess_mmr([], msi_pct=None)
        assert result.is_mmrd is False

    def test_msi_status_discordance_flag(self):
        """Flag when computed MSI status differs from override."""
        result = assess_mmr([], msi_pct=55.0, msi_status_override="MSS")
        assert result.is_mmrd is True
        assert any("discordance" in f.lower() for f in result.flags)

    def test_custom_threshold(self):
        result = assess_mmr([], msi_pct=15.0, msi_threshold=10.0)
        assert result.is_mmrd is True


class TestFindMmrVariants:
    def test_filters_mmr_genes(self):
        variants = [
            make_variant(hugo_symbol="MSH6", hgvsp_short="p.F1088fs"),
            make_variant(hugo_symbol="TP53", hgvsp_short="p.R175H"),
            make_variant(hugo_symbol="MLH1", hgvsp_short="p.R100X"),
        ]
        mmr = find_mmr_variants(variants)
        assert len(mmr) == 2
        assert all(v.hugo_symbol in ("MSH6", "MLH1") for v in mmr)
