"""Pydantic models for ec-molsubtype classification pipeline."""

from __future__ import annotations

from enum import Enum
from typing import Any

from pydantic import BaseModel, Field


class MolecularSubtype(str, Enum):
    """The four molecular subtypes of endometrial cancer (WHO 2020)."""

    POLEmut = "POLEmut"
    MMRd = "MMRd"
    p53abn = "p53abn"
    NSMP = "NSMP"


class ConfidenceLevel(str, Enum):
    """Confidence level for the classification."""

    HIGH = "high"
    MODERATE = "moderate"
    LOW = "low"
    DISCORDANT = "discordant"


class PoleTier(str, Enum):
    """POLE variant pathogenicity tier."""

    TIER1 = "tier1"
    TIER2 = "tier2"
    VUS = "vus"
    NOT_EDM = "not_edm"


class VariantClassification(str, Enum):
    """MAF variant classification types."""

    MISSENSE = "Missense_Mutation"
    NONSENSE = "Nonsense_Mutation"
    FRAME_SHIFT_DEL = "Frame_Shift_Del"
    FRAME_SHIFT_INS = "Frame_Shift_Ins"
    SPLICE_SITE = "Splice_Site"
    IN_FRAME_DEL = "In_Frame_Del"
    IN_FRAME_INS = "In_Frame_Ins"
    TRANSLATION_START_SITE = "Translation_Start_Site"
    NONSTOP_MUTATION = "Nonstop_Mutation"
    SILENT = "Silent"
    RNA = "RNA"
    INTRON = "Intron"
    THREE_PRIME_UTR = "3'UTR"
    FIVE_PRIME_UTR = "5'UTR"
    IGR = "IGR"
    FIVE_PRIME_FLANK = "5'Flank"
    THREE_PRIME_FLANK = "3'Flank"


TRUNCATING_CLASSIFICATIONS = {
    VariantClassification.NONSENSE,
    VariantClassification.FRAME_SHIFT_DEL,
    VariantClassification.FRAME_SHIFT_INS,
    VariantClassification.SPLICE_SITE,
    VariantClassification.NONSTOP_MUTATION,
}


class Variant(BaseModel):
    """A single somatic variant from MAF/VCF."""

    hugo_symbol: str
    chromosome: str = ""
    start_position: int = 0
    end_position: int = 0
    reference_allele: str = ""
    tumor_seq_allele2: str = ""
    variant_classification: str = ""
    hgvsp_short: str = ""
    t_alt_count: int = 0
    t_ref_count: int = 0
    clinvar_classification: str = ""
    sift: str = ""
    polyphen: str = ""

    @property
    def vaf(self) -> float | None:
        """Variant allele frequency."""
        total = self.t_alt_count + self.t_ref_count
        if total == 0:
            return None
        return self.t_alt_count / total

    @property
    def is_truncating(self) -> bool:
        """Whether this is a truncating variant."""
        try:
            vc = VariantClassification(self.variant_classification)
            return vc in TRUNCATING_CLASSIFICATIONS
        except ValueError:
            return False

    @property
    def protein_change(self) -> str:
        """Normalized protein change (strip leading 'p.')."""
        hgvs = self.hgvsp_short.strip()
        if hgvs.startswith("p."):
            return hgvs[2:]
        return hgvs


class SampleMetadata(BaseModel):
    """Per-sample supplementary metadata."""

    sample_id: str
    tmb: float | None = None
    panel_size_mb: float | None = None
    msi_pct: float | None = Field(None, ge=0.0, le=100.0)
    msi_status: str | None = None
    fraction_genome_altered: float | None = Field(None, ge=0.0, le=1.0)
    signature_weights: dict[str, float] | None = None


class SampleInput(BaseModel):
    """Complete input for a single sample: variants + metadata."""

    metadata: SampleMetadata
    variants: list[Variant] = Field(default_factory=list)


class ClassificationStep(BaseModel):
    """A single step in the classification algorithm."""

    step: int
    test: str
    result: str  # "positive" or "negative"
    variant: str | None = None
    tier: str | None = None
    msi_pct: float | None = None
    msi_threshold: float | None = None
    mmr_mutations: list[str] | None = None
    confidence: str | None = None
    details: str | None = None


class MultipleClassifier(BaseModel):
    """Information about multiple-classifier status."""

    is_multiple: bool = False
    secondary_features: list[str] = Field(default_factory=list)
    tp53_variant: str | None = None
    mmr_evidence: str | None = None
    pole_variant: str | None = None


class EvidenceItem(BaseModel):
    """A single piece of secondary evidence."""

    value: float | None = None
    expected: str = ""
    expected_range: str = ""
    concordant: bool | None = None
    details: dict[str, Any] = Field(default_factory=dict)


class SecondaryEvidence(BaseModel):
    """Aggregated secondary evidence for classification."""

    tmb: EvidenceItem | None = None
    cna_burden: EvidenceItem | None = None
    signatures: EvidenceItem | None = None
    substitution_profile: EvidenceItem | None = None


class ClassificationResult(BaseModel):
    """Complete classification result for a single sample."""

    sample_id: str
    primary_subtype: MolecularSubtype
    confidence: ConfidenceLevel
    classification_path: list[ClassificationStep] = Field(default_factory=list)
    multiple_classifier: MultipleClassifier = Field(default_factory=MultipleClassifier)
    secondary_evidence: SecondaryEvidence = Field(default_factory=SecondaryEvidence)
    figo_molecular_annotation: str = ""
    clinical_notes: list[str] = Field(default_factory=list)
    flags: list[str] = Field(default_factory=list)
    version: str = "0.1.0"
