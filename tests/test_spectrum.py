"""Tests for substitution spectrum analysis."""

from ec_molsubtype.models import MolecularSubtype, Variant
from ec_molsubtype.spectrum import (
    SBS_CLASSES,
    assess_spectrum,
    compute_spectrum,
)


def _snv(ref: str, alt: str, vc: str = "Missense_Mutation") -> Variant:
    """Create a minimal SNV variant."""
    return Variant(
        hugo_symbol="GENE",
        reference_allele=ref,
        tumor_seq_allele2=alt,
        variant_classification=vc,
    )


def _indel(vc: str = "Frame_Shift_Del") -> Variant:
    """Create a minimal indel variant."""
    return Variant(
        hugo_symbol="GENE",
        reference_allele="A",
        tumor_seq_allele2="-",
        variant_classification=vc,
    )


class TestComputeSpectrum:
    def test_empty_variants(self):
        result = compute_spectrum([])
        assert result.n_snvs == 0
        assert result.n_indels == 0
        assert result.indel_fraction == 0.0
        assert all(v == 0.0 for v in result.fractions.values())

    def test_single_c_to_a(self):
        variants = [_snv("C", "A")]
        result = compute_spectrum(variants)
        assert result.n_snvs == 1
        assert result.counts["C>A"] == 1
        assert result.fractions["C>A"] == 1.0

    def test_complement_g_to_t_maps_to_c_to_a(self):
        """G>T on + strand = C>A on pyrimidine strand."""
        variants = [_snv("G", "T")]
        result = compute_spectrum(variants)
        assert result.counts["C>A"] == 1

    def test_complement_a_to_c_maps_to_t_to_g(self):
        """A>C on + strand = T>G on pyrimidine strand."""
        variants = [_snv("A", "C")]
        result = compute_spectrum(variants)
        assert result.counts["T>G"] == 1

    def test_all_six_classes(self):
        variants = [
            _snv("C", "A"),  # C>A
            _snv("C", "G"),  # C>G
            _snv("C", "T"),  # C>T
            _snv("T", "A"),  # T>A
            _snv("T", "C"),  # T>C
            _snv("T", "G"),  # T>G
        ]
        result = compute_spectrum(variants)
        assert result.n_snvs == 6
        for cls in SBS_CLASSES:
            assert result.counts[cls] == 1
            assert abs(result.fractions[cls] - 1 / 6) < 0.001

    def test_indels_counted(self):
        variants = [
            _snv("C", "A"),
            _indel("Frame_Shift_Del"),
            _indel("Frame_Shift_Ins"),
            _indel("In_Frame_Del"),
        ]
        result = compute_spectrum(variants)
        assert result.n_snvs == 1
        assert result.n_indels == 3
        assert abs(result.indel_fraction - 0.75) < 0.001

    def test_non_coding_ignored(self):
        variants = [
            _snv("C", "A"),
            Variant(
                hugo_symbol="GENE",
                reference_allele="C",
                tumor_seq_allele2="T",
                variant_classification="Intron",
            ),
            Variant(
                hugo_symbol="GENE",
                reference_allele="C",
                tumor_seq_allele2="T",
                variant_classification="3'UTR",
            ),
        ]
        result = compute_spectrum(variants)
        assert result.n_snvs == 1
        assert result.n_total_variants == 1

    def test_silent_counted_in_spectrum(self):
        """Silent variants are SNVs and contribute to the spectrum."""
        variants = [_snv("C", "T", vc="Silent")]
        result = compute_spectrum(variants)
        assert result.n_snvs == 1
        assert result.counts["C>T"] == 1


class TestAssessSpectrum:
    def _pole_like_variants(self, n: int = 100) -> list[Variant]:
        """Generate a POLE-like spectrum: high C>A, some T>G, very low C>G."""
        variants = []
        # 35% C>A
        for _ in range(int(n * 0.35)):
            variants.append(_snv("C", "A"))
        # 30% C>T
        for _ in range(int(n * 0.30)):
            variants.append(_snv("C", "T"))
        # 10% T>G
        for _ in range(int(n * 0.10)):
            variants.append(_snv("T", "G"))
        # 15% T>C
        for _ in range(int(n * 0.15)):
            variants.append(_snv("T", "C"))
        # 10% T>A
        for _ in range(int(n * 0.10)):
            variants.append(_snv("T", "A"))
        # 0% C>G → POLE hallmark
        # 2 indels out of ~102 total = ~2% indel fraction
        variants.append(_indel())
        variants.append(_indel())
        return variants

    def _mmrd_like_variants(self, n: int = 50) -> list[Variant]:
        """Generate an MMRd-like spectrum: high indel fraction."""
        variants = []
        # Mix of SNVs
        for _ in range(int(n * 0.4)):
            variants.append(_snv("C", "T"))
        for _ in range(int(n * 0.3)):
            variants.append(_snv("T", "C"))
        for _ in range(int(n * 0.3)):
            variants.append(_snv("C", "A"))
        # High indel fraction: 15 indels out of 65 total = 23%
        for _ in range(15):
            variants.append(_indel())
        return variants

    def _p53abn_like_variants(self) -> list[Variant]:
        """Generate a p53abn-like spectrum: low TMB, no distinctive features."""
        return [
            _snv("C", "T"),
            _snv("C", "T"),
            _snv("C", "A"),
            _snv("T", "C"),
            _snv("C", "G"),
        ]

    def test_pole_spectrum_concordant(self):
        variants = self._pole_like_variants()
        result = assess_spectrum(variants, MolecularSubtype.POLEmut)
        assert result.concordant is True
        assert result.checks["C>A"]["passed"] is True
        assert result.checks["T>G"]["passed"] is True
        assert result.checks["C>G"]["passed"] is True
        assert result.checks["indel_fraction"]["passed"] is True

    def test_pole_spectrum_discordant_low_ca(self):
        """Non-POLE spectrum should be discordant with POLEmut."""
        variants = [_snv("C", "T")] * 50 + [_snv("C", "A")] * 5
        result = assess_spectrum(variants, MolecularSubtype.POLEmut)
        assert result.concordant is False
        assert result.checks["C>A"]["passed"] is False

    def test_mmrd_spectrum_concordant(self):
        variants = self._mmrd_like_variants()
        result = assess_spectrum(variants, MolecularSubtype.MMRd)
        assert result.concordant is True
        assert result.checks["indel_fraction"]["passed"] is True

    def test_mmrd_spectrum_discordant_low_indels(self):
        """Low indel fraction should be discordant with MMRd."""
        variants = [_snv("C", "T")] * 50 + [_indel()]  # 1/51 = 2%
        result = assess_spectrum(variants, MolecularSubtype.MMRd)
        assert result.concordant is False

    def test_insufficient_snvs(self):
        variants = [_snv("C", "A")] * 5
        result = assess_spectrum(variants, MolecularSubtype.POLEmut)
        assert result.concordant is None
        assert "Insufficient" in result.details

    def test_no_profile_for_nsmp(self):
        """NSMP has no substitution profile thresholds."""
        variants = self._p53abn_like_variants() * 3
        result = assess_spectrum(variants, MolecularSubtype.NSMP)
        assert result.concordant is None

    def test_custom_min_snvs(self):
        variants = [_snv("C", "A")] * 5
        result = assess_spectrum(variants, MolecularSubtype.POLEmut, min_snvs=3)
        assert result.concordant is not None  # now has enough
