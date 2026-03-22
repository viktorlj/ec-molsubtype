"""Tests for MAF/VCF parsing."""

import json
import pytest
from pathlib import Path

from ec_molsubtype.io import load_maf_variants, load_sample, load_sample_metadata, parse_maf


TEST_DATA_DIR = Path(__file__).parent / "data"


@pytest.fixture
def mini_maf(tmp_path):
    """Create a minimal MAF file for testing."""
    content = (
        "Hugo_Symbol\tChromosome\tStart_Position\tEnd_Position\t"
        "Reference_Allele\tTumor_Seq_Allele2\tVariant_Classification\t"
        "HGVSp_Short\tt_alt_count\tt_ref_count\tTumor_Sample_Barcode\n"
        "POLE\tchr12\t133250250\t133250250\tC\tG\tMissense_Mutation\t"
        "p.P286R\t80\t20\tSAMPLE_001\n"
        "TP53\tchr17\t7577120\t7577120\tC\tT\tMissense_Mutation\t"
        "p.R248W\t60\t40\tSAMPLE_001\n"
        "MSH6\tchr2\t47641560\t47641560\tG\t-\tFrame_Shift_Del\t"
        "p.F1088fs\t45\t55\tSAMPLE_001\n"
    )
    maf_path = tmp_path / "test.maf"
    maf_path.write_text(content)
    return maf_path


@pytest.fixture
def mini_metadata(tmp_path):
    """Create a minimal metadata JSON file."""
    data = {
        "sample_id": "SAMPLE_001",
        "tmb": 250.0,
        "panel_size_mb": 1.5,
        "msi_pct": 2.0,
        "fraction_genome_altered": 0.05,
        "signature_weights": {"SBS10a": 0.4, "SBS10b": 0.2},
    }
    meta_path = tmp_path / "metadata.json"
    meta_path.write_text(json.dumps(data))
    return meta_path


class TestParseMaf:
    def test_parse_basic(self, mini_maf):
        df = parse_maf(mini_maf)
        assert len(df) == 3
        assert "Hugo_Symbol" in df.columns
        assert "HGVSp_Short" in df.columns

    def test_parse_with_comments(self, tmp_path):
        content = (
            "# version 2.4\n"
            "# comment line\n"
            "Hugo_Symbol\tVariant_Classification\tHGVSp_Short\n"
            "POLE\tMissense_Mutation\tp.P286R\n"
        )
        path = tmp_path / "commented.maf"
        path.write_text(content)
        df = parse_maf(path)
        assert len(df) == 1

    def test_missing_columns_raises(self, tmp_path):
        content = "Hugo_Symbol\tSomeOther\nPOLE\tvalue\n"
        path = tmp_path / "bad.maf"
        path.write_text(content)
        with pytest.raises(ValueError, match="missing required columns"):
            parse_maf(path)


class TestLoadMafVariants:
    def test_loads_variants(self, mini_maf):
        variants = load_maf_variants(mini_maf)
        assert len(variants) == 3
        assert variants[0].hugo_symbol == "POLE"
        assert variants[0].hgvsp_short == "p.P286R"


class TestLoadSampleMetadata:
    def test_loads_metadata(self, mini_metadata):
        meta = load_sample_metadata(mini_metadata)
        assert meta.sample_id == "SAMPLE_001"
        assert meta.tmb == 250.0
        assert meta.msi_pct == 2.0


class TestLoadSample:
    def test_loads_with_metadata(self, mini_maf, mini_metadata):
        sample = load_sample(mini_maf, metadata_path=mini_metadata)
        assert sample.metadata.sample_id == "SAMPLE_001"
        assert len(sample.variants) == 3

    def test_loads_without_metadata(self, mini_maf):
        sample = load_sample(mini_maf)
        assert sample.metadata.sample_id == "SAMPLE_001"  # from Tumor_Sample_Barcode
        assert len(sample.variants) == 3
