"""
Tests for AnnoRefine Python bindings
"""

import pytest
import tempfile
import os
from pathlib import Path

try:
    import annorefine
    ANNOREFINE_AVAILABLE = True
except ImportError:
    ANNOREFINE_AVAILABLE = False


def test_basic_imports():
    """Test that basic Python modules work"""
    import sys
    import os
    assert sys.version_info.major == 3
    assert os.path.exists(".")


@pytest.mark.skipif(not ANNOREFINE_AVAILABLE, reason="AnnoRefine Python bindings not available")
class TestAnnoRefine:
    """Test suite for AnnoRefine Python bindings"""

    def test_version(self):
        """Test version function"""
        version = annorefine.version()
        assert isinstance(version, str)
        assert len(version) > 0
        # Should follow CalVer format YYYY.M.D
        parts = version.split('.')
        assert len(parts) == 3
        assert int(parts[0]) >= 2025  # Year should be 2025 or later

    def test_refinement_config_default(self):
        """Test RefinementConfig with default values"""
        config = annorefine.RefinementConfig()
        
        assert config.min_coverage == 5
        assert config.min_splice_support == 3
        assert config.max_utr_extension == 1000
        assert config.enable_novel_gene_detection == False
        assert config.min_novel_gene_coverage == 10
        assert config.min_novel_gene_length == 300
        assert config.min_exon_length == 50
        assert config.validate_splice_sites == True

    def test_refinement_config_custom(self):
        """Test RefinementConfig with custom values"""
        config = annorefine.RefinementConfig(
            min_coverage=10,
            min_splice_support=5,
            max_utr_extension=2000,
            enable_novel_gene_detection=True,
            min_novel_gene_coverage=15,
            min_novel_gene_length=500,
            min_exon_length=100,
            validate_splice_sites=False
        )
        
        assert config.min_coverage == 10
        assert config.min_splice_support == 5
        assert config.max_utr_extension == 2000
        assert config.enable_novel_gene_detection == True
        assert config.min_novel_gene_coverage == 15
        assert config.min_novel_gene_length == 500
        assert config.min_exon_length == 100
        assert config.validate_splice_sites == False

    def test_refinement_config_repr(self):
        """Test RefinementConfig string representation"""
        config = annorefine.RefinementConfig(min_coverage=10)
        repr_str = repr(config)
        assert "RefinementConfig" in repr_str
        assert "min_coverage=10" in repr_str

    def test_gene_model_repr(self):
        """Test GeneModel string representation"""
        # This would need actual gene model data from a real run
        # For now, just test that the class exists
        assert hasattr(annorefine, 'GeneModel')

    def test_refine_function_signature(self):
        """Test that refine function exists with correct signature"""
        assert hasattr(annorefine, 'refine')
        assert callable(annorefine.refine)

    def test_refine_annotations_function_signature(self):
        """Test that refine_annotations function exists with correct signature"""
        assert hasattr(annorefine, 'refine_annotations')
        assert callable(annorefine.refine_annotations)

    def test_module_metadata(self):
        """Test module metadata"""
        assert hasattr(annorefine, '__version__')
        assert hasattr(annorefine, '__author__')
        assert hasattr(annorefine, '__description__')
        
        assert annorefine.__author__ == "Jon Palmer"
        assert "annotation" in annorefine.__description__.lower()

    @pytest.mark.skipif(
        not os.path.exists("test_data") or os.getenv("CI"),
        reason="Test data directory not available or running in CI"
    )
    def test_refine_with_test_data(self):
        """Test refinement with actual test data (if available)"""
        test_data_dir = Path("test_data")
        
        # Look for test files
        fasta_files = list(test_data_dir.glob("*.fasta")) + list(test_data_dir.glob("*.fa"))
        gff3_files = list(test_data_dir.glob("*.gff3")) + list(test_data_dir.glob("*.gff"))
        bam_files = list(test_data_dir.glob("*.bam"))
        
        if not (fasta_files and gff3_files and bam_files):
            pytest.skip("Required test data files not found")
        
        with tempfile.NamedTemporaryFile(suffix=".gff3", delete=False) as tmp_output:
            try:
                result = annorefine.refine(
                    fasta_file=str(fasta_files[0]),
                    gff3_file=str(gff3_files[0]),
                    bam_file=str(bam_files[0]),
                    output_file=tmp_output.name,
                    min_coverage=1,  # Low threshold for test data
                    threads=1
                )
                
                # Check result structure
                assert isinstance(result, dict)
                assert 'genes_processed' in result
                assert 'genes_failed' in result
                assert 'output_file' in result
                assert result['output_file'] == tmp_output.name
                
                # Check that output file was created
                assert os.path.exists(tmp_output.name)
                assert os.path.getsize(tmp_output.name) > 0
                
            finally:
                # Clean up
                if os.path.exists(tmp_output.name):
                    os.unlink(tmp_output.name)

    def test_error_handling_missing_files(self):
        """Test error handling for missing input files"""
        with tempfile.NamedTemporaryFile(suffix=".gff3") as tmp_output:
            with pytest.raises(Exception):  # Should raise some kind of error
                annorefine.refine(
                    fasta_file="nonexistent.fasta",
                    gff3_file="nonexistent.gff3",
                    bam_file="nonexistent.bam",
                    output_file=tmp_output.name
                )

    def test_config_parameter_types(self):
        """Test that config parameters accept correct types"""
        # Test integer parameters
        config = annorefine.RefinementConfig(min_coverage=10)
        assert config.min_coverage == 10
        
        # Test boolean parameters
        config = annorefine.RefinementConfig(enable_novel_gene_detection=True)
        assert config.enable_novel_gene_detection == True
        
        config = annorefine.RefinementConfig(validate_splice_sites=False)
        assert config.validate_splice_sites == False


if __name__ == "__main__":
    pytest.main([__file__])
