"""Tests for POLE EDM assessment."""

import pytest
from ec_molsubtype.models import PoleTier, Variant
from ec_molsubtype.pole import assess_pole, check_pole_variant, get_best_pole_result
from tests.conftest import make_variant


class TestCheckPoleVariant:
    def test_tier1_p286r(self, pole_p286r):
        result = check_pole_variant(pole_p286r)
        assert result.is_pathogenic is True
        assert result.tier == PoleTier.TIER1
        assert "P286R" in result.variant_str
        assert result.needs_secondary_evidence is False

    def test_tier1_v411l(self, pole_v411l):
        result = check_pole_variant(pole_v411l)
        assert result.is_pathogenic is True
        assert result.tier == PoleTier.TIER1

    def test_tier1_all_hotspots(self):
        for var_str in ["p.P286R", "p.V411L", "p.S297F", "p.A456P", "p.S459F"]:
            v = make_variant(hugo_symbol="POLE", hgvsp_short=var_str)
            result = check_pole_variant(v)
            assert result.is_pathogenic is True, f"{var_str} should be tier1"
            assert result.tier == PoleTier.TIER1, f"{var_str} should be tier1"

    def test_tier2_f367s(self, pole_tier2_f367s):
        result = check_pole_variant(pole_tier2_f367s)
        assert result.is_pathogenic is True
        assert result.tier == PoleTier.TIER2

    def test_tier2_all(self):
        for var_str in ["p.F367S", "p.L424I", "p.L424V", "p.M444K", "p.D368N", "p.P286H"]:
            v = make_variant(hugo_symbol="POLE", hgvsp_short=var_str)
            result = check_pole_variant(v)
            assert result.is_pathogenic is True, f"{var_str} should be tier2"
            assert result.tier == PoleTier.TIER2, f"{var_str} should be tier2"

    def test_vus_in_edm(self, pole_vus):
        result = check_pole_variant(pole_vus)
        assert result.is_pathogenic is False
        assert result.tier == PoleTier.VUS
        assert result.needs_secondary_evidence is True

    def test_outside_edm(self, pole_outside_edm):
        result = check_pole_variant(pole_outside_edm)
        assert result.is_pathogenic is False
        assert result.tier == PoleTier.NOT_EDM
        assert result.needs_secondary_evidence is False

    def test_non_pole_gene(self):
        v = make_variant(hugo_symbol="TP53", hgvsp_short="p.R175H")
        result = check_pole_variant(v)
        assert result.is_pathogenic is False
        assert result.tier == PoleTier.NOT_EDM

    def test_unparseable_hgvsp(self):
        v = make_variant(hugo_symbol="POLE", hgvsp_short="weird_format")
        result = check_pole_variant(v)
        assert result.is_pathogenic is False


class TestAssessPole:
    def test_multiple_pole_variants(self):
        variants = [
            make_variant(hugo_symbol="POLE", hgvsp_short="p.P286R"),
            make_variant(hugo_symbol="POLE", hgvsp_short="p.F367S"),
            make_variant(hugo_symbol="TP53", hgvsp_short="p.R175H"),
        ]
        results = assess_pole(variants)
        assert len(results) == 2

    def test_no_pole_variants(self):
        variants = [make_variant(hugo_symbol="TP53", hgvsp_short="p.R175H")]
        results = assess_pole(variants)
        assert len(results) == 0


class TestGetBestPoleResult:
    def test_tier1_wins(self):
        variants = [
            make_variant(hugo_symbol="POLE", hgvsp_short="p.P286R"),  # tier1
            make_variant(hugo_symbol="POLE", hgvsp_short="p.F367S"),  # tier2
        ]
        results = assess_pole(variants)
        best = get_best_pole_result(results)
        assert best is not None
        assert best.tier == PoleTier.TIER1

    def test_empty(self):
        assert get_best_pole_result([]) is None
