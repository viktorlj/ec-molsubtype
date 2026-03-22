"""Tests for TP53 variant interpretation."""

import pytest
from ec_molsubtype.tp53 import assess_tp53, check_tp53_variant, get_pathogenic_tp53
from tests.conftest import make_variant


class TestCheckTp53Variant:
    def test_hotspot_r175h(self, tp53_r175h):
        result = check_tp53_variant(tp53_r175h)
        assert result.is_pathogenic is True
        assert result.is_hotspot is True
        assert "R175H" in result.variant_str

    def test_hotspot_r248w(self, tp53_r248w):
        result = check_tp53_variant(tp53_r248w)
        assert result.is_pathogenic is True
        assert result.is_hotspot is True

    def test_all_hotspots(self):
        for var_str in ["p.R175H", "p.G245S", "p.R248W", "p.R248Q", "p.R249S", "p.R273H", "p.R273C", "p.R282W"]:
            v = make_variant(hugo_symbol="TP53", hgvsp_short=var_str)
            result = check_tp53_variant(v)
            assert result.is_pathogenic is True, f"{var_str} should be pathogenic"
            assert result.is_hotspot is True, f"{var_str} should be hotspot"

    def test_truncating(self, tp53_truncating):
        result = check_tp53_variant(tp53_truncating)
        assert result.is_pathogenic is True
        assert result.is_truncating is True
        assert result.is_hotspot is False

    def test_benign_p72r(self, tp53_p72r):
        result = check_tp53_variant(tp53_p72r)
        assert result.is_pathogenic is False
        assert result.is_benign_polymorphism is True

    def test_low_vaf_filtered(self, tp53_low_vaf):
        result = check_tp53_variant(tp53_low_vaf, min_vaf=0.05)
        assert result.is_pathogenic is False
        assert any("low VAF" in n for n in result.clinical_notes)

    def test_clinvar_pathogenic(self):
        v = make_variant(
            hugo_symbol="TP53",
            hgvsp_short="p.Y220C",
            clinvar_classification="Likely_pathogenic",
        )
        result = check_tp53_variant(v)
        assert result.is_pathogenic is True

    def test_non_tp53_gene(self):
        v = make_variant(hugo_symbol="POLE", hgvsp_short="p.P286R")
        result = check_tp53_variant(v)
        assert result.is_pathogenic is False
        assert result.variant_str == ""

    def test_in_dna_binding_domain(self):
        v = make_variant(hugo_symbol="TP53", hgvsp_short="p.R175H")
        result = check_tp53_variant(v)
        assert result.in_dna_binding_domain is True

    def test_outside_dna_binding_domain(self):
        v = make_variant(hugo_symbol="TP53", hgvsp_short="p.R342H")
        result = check_tp53_variant(v)
        assert result.in_dna_binding_domain is False


class TestAssessTp53:
    def test_multiple_tp53(self):
        variants = [
            make_variant(hugo_symbol="TP53", hgvsp_short="p.R175H"),
            make_variant(hugo_symbol="TP53", hgvsp_short="p.R248W"),
        ]
        results = assess_tp53(variants)
        assert len(results) == 2

    def test_no_tp53(self):
        variants = [make_variant(hugo_symbol="PTEN", hgvsp_short="p.R130Q")]
        results = assess_tp53(variants)
        assert len(results) == 0


class TestGetPathogenicTp53:
    def test_hotspot_preferred(self):
        variants = [
            make_variant(hugo_symbol="TP53", hgvsp_short="p.R175H"),
            make_variant(
                hugo_symbol="TP53",
                hgvsp_short="p.Q100X",
                variant_classification="Nonsense_Mutation",
            ),
        ]
        results = assess_tp53(variants)
        best = get_pathogenic_tp53(results)
        assert best is not None
        assert best.is_hotspot is True

    def test_none_when_no_pathogenic(self):
        variants = [make_variant(hugo_symbol="TP53", hgvsp_short="p.P72R")]
        results = assess_tp53(variants)
        best = get_pathogenic_tp53(results)
        assert best is None
