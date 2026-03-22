"""End-to-end classification tests."""

import pytest
from ec_molsubtype.classify import classify_sample
from ec_molsubtype.models import (
    ConfidenceLevel,
    MolecularSubtype,
    SampleInput,
    SampleMetadata,
)
from tests.conftest import make_variant


class TestClassifyPOLEmut:
    def test_tier1_polemut(self, polemut_sample):
        result = classify_sample(polemut_sample)
        assert result.primary_subtype == MolecularSubtype.POLEmut
        assert result.confidence == ConfidenceLevel.HIGH

    def test_polemut_overrides_tp53(self, polemut_sample):
        """POLEmut + p53abn → still POLEmut (hierarchy)."""
        result = classify_sample(polemut_sample)
        assert result.primary_subtype == MolecularSubtype.POLEmut
        assert result.multiple_classifier.is_multiple is True
        assert "p53abn" in result.multiple_classifier.secondary_features

    def test_polemut_overrides_mmrd(self):
        """POLEmut + MMRd → still POLEmut."""
        sample = SampleInput(
            metadata=SampleMetadata(
                sample_id="POLE_MMR",
                tmb=200.0,
                msi_pct=55.0,
                signature_weights={"SBS10a": 0.4, "SBS10b": 0.2},
            ),
            variants=[
                make_variant(hugo_symbol="POLE", hgvsp_short="p.V411L"),
                make_variant(
                    hugo_symbol="MSH6",
                    hgvsp_short="p.F1088fs",
                    variant_classification="Frame_Shift_Del",
                ),
            ],
        )
        result = classify_sample(sample)
        assert result.primary_subtype == MolecularSubtype.POLEmut

    def test_pole_vus_not_auto_classified(self):
        """POLE VUS should NOT be classified as POLEmut."""
        sample = SampleInput(
            metadata=SampleMetadata(sample_id="POLE_VUS", tmb=150.0),
            variants=[make_variant(hugo_symbol="POLE", hgvsp_short="p.A300T")],
        )
        result = classify_sample(sample)
        assert result.primary_subtype != MolecularSubtype.POLEmut
        assert any("VUS" in f for f in result.flags)


class TestClassifyMMRd:
    def test_msi_h_mmrd(self, mmrd_sample):
        result = classify_sample(mmrd_sample)
        assert result.primary_subtype == MolecularSubtype.MMRd
        assert result.confidence in (ConfidenceLevel.HIGH, ConfidenceLevel.MODERATE)

    def test_mmrd_overrides_tp53(self):
        """MMRd + p53abn → still MMRd."""
        sample = SampleInput(
            metadata=SampleMetadata(
                sample_id="MMR_P53",
                tmb=40.0,
                msi_pct=45.0,
                signature_weights={"SBS6": 0.3},
            ),
            variants=[
                make_variant(hugo_symbol="TP53", hgvsp_short="p.R175H"),
            ],
        )
        result = classify_sample(sample)
        assert result.primary_subtype == MolecularSubtype.MMRd
        assert result.multiple_classifier.is_multiple is True
        assert "p53abn" in result.multiple_classifier.secondary_features

    def test_biallelic_mmr_no_msi(self):
        """Biallelic MMR mutations without MSI data → moderate confidence."""
        sample = SampleInput(
            metadata=SampleMetadata(sample_id="BIALLELIC", tmb=30.0),
            variants=[
                make_variant(
                    hugo_symbol="MLH1",
                    hgvsp_short="p.R100X",
                    variant_classification="Nonsense_Mutation",
                ),
                make_variant(
                    hugo_symbol="MLH1",
                    hgvsp_short="p.W200fs",
                    variant_classification="Frame_Shift_Del",
                ),
            ],
        )
        result = classify_sample(sample)
        assert result.primary_subtype == MolecularSubtype.MMRd
        assert result.confidence in (ConfidenceLevel.MODERATE, ConfidenceLevel.HIGH)


class TestClassifyP53abn:
    def test_p53abn(self, p53abn_sample):
        result = classify_sample(p53abn_sample)
        assert result.primary_subtype == MolecularSubtype.p53abn

    def test_p53_hotspot(self):
        sample = SampleInput(
            metadata=SampleMetadata(
                sample_id="P53_HOTSPOT",
                tmb=4.0,
                msi_pct=0.0,
                fraction_genome_altered=0.5,
                signature_weights={"SBS1": 0.4, "SBS5": 0.3},
            ),
            variants=[make_variant(hugo_symbol="TP53", hgvsp_short="p.R248Q")],
        )
        result = classify_sample(sample)
        assert result.primary_subtype == MolecularSubtype.p53abn

    def test_p53_truncating(self):
        sample = SampleInput(
            metadata=SampleMetadata(sample_id="P53_TRUNC", msi_pct=0.0),
            variants=[
                make_variant(
                    hugo_symbol="TP53",
                    hgvsp_short="p.E180fs",
                    variant_classification="Frame_Shift_Del",
                ),
            ],
        )
        result = classify_sample(sample)
        assert result.primary_subtype == MolecularSubtype.p53abn


class TestClassifyNSMP:
    def test_nsmp(self, nsmp_sample):
        result = classify_sample(nsmp_sample)
        assert result.primary_subtype == MolecularSubtype.NSMP

    def test_nsmp_no_variants(self):
        sample = SampleInput(
            metadata=SampleMetadata(sample_id="EMPTY", msi_pct=0.0),
            variants=[],
        )
        result = classify_sample(sample)
        assert result.primary_subtype == MolecularSubtype.NSMP

    def test_nsmp_benign_tp53(self):
        """P72R should not trigger p53abn."""
        sample = SampleInput(
            metadata=SampleMetadata(sample_id="BENIGN_TP53", msi_pct=0.0),
            variants=[make_variant(hugo_symbol="TP53", hgvsp_short="p.P72R")],
        )
        result = classify_sample(sample)
        assert result.primary_subtype == MolecularSubtype.NSMP


class TestClassificationPath:
    def test_path_recorded(self, polemut_sample):
        result = classify_sample(polemut_sample)
        assert len(result.classification_path) >= 1
        assert result.classification_path[0].step == 1
        assert result.classification_path[0].test == "POLE_EDM"
        assert result.classification_path[0].result == "positive"

    def test_nsmp_has_all_steps(self, nsmp_sample):
        result = classify_sample(nsmp_sample)
        assert len(result.classification_path) == 4
        steps = [s.step for s in result.classification_path]
        assert steps == [1, 2, 3, 4]


class TestSecondaryEvidence:
    def test_tmb_concordance_polemut(self, polemut_sample):
        result = classify_sample(polemut_sample)
        assert result.secondary_evidence.tmb is not None
        assert result.secondary_evidence.tmb.concordant is True

    def test_cna_concordance_p53abn(self, p53abn_sample):
        result = classify_sample(p53abn_sample)
        assert result.secondary_evidence.cna_burden is not None
        assert result.secondary_evidence.cna_burden.concordant is True

    def test_signature_concordance_polemut(self, polemut_sample):
        result = classify_sample(polemut_sample)
        assert result.secondary_evidence.signatures is not None
        assert result.secondary_evidence.signatures.concordant is True


class TestClinicalNotes:
    def test_polemut_swedish_guidelines(self, polemut_sample):
        result = classify_sample(polemut_sample)
        assert any("Swedish guidelines" in n for n in result.clinical_notes)
        assert any("no adjuvant" in n for n in result.clinical_notes)

    def test_p53abn_swedish_guidelines(self, p53abn_sample):
        result = classify_sample(p53abn_sample)
        assert any("chemotherapy" in n for n in result.clinical_notes)

    def test_mmrd_ihc_recommendation(self, mmrd_sample):
        result = classify_sample(mmrd_sample)
        assert any("IHC" in n or "Lynch" in n for n in result.clinical_notes)
