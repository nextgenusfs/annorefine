"""
AnnoRefine: Genome annotation refinement using RNA-seq data

This package provides Python bindings for the AnnoRefine Rust library,
enabling genome annotation refinement directly from Python.

Example usage:
    import annorefine
    
    # Basic refinement
    result = annorefine.refine_annotations(
        fasta_file="genome.fasta",
        gff3_file="annotations.gff3", 
        bam_file="alignments.bam",
        output_file="refined.gff3"
    )
    
    # Advanced configuration
    config = annorefine.RefinementConfig(
        min_coverage=10,
        enable_novel_gene_detection=True,
        validate_splice_sites=True
    )
    
    result = annorefine.refine_annotations(
        fasta_file="genome.fasta",
        gff3_file="annotations.gff3",
        bam_file="alignments.bam", 
        output_file="refined.gff3",
        config=config,
        threads=8
    )
    
    print(f"Processed {result['genes_processed']} genes")
    print(f"Found {result['novel_genes_detected']} novel genes")
"""

from ._annorefine import (
    refine_annotations,
    version,
    current_num_threads,
    test_interruptible_operation,
    PyRefinementConfig as RefinementConfig,
    PyGeneModel as GeneModel,
)

__version__ = version()
__author__ = "Jon Palmer"
__email__ = "nextgenusfs@gmail.com"
__description__ = "Genome annotation refinement using RNA-seq data"

__all__ = [
    "refine",
    "refine_annotations",
    "version",
    "current_num_threads",
    "test_interruptible_operation",
    "RefinementConfig",
    "GeneModel",
]


def refine(
    fasta_file: str,
    gff3_file: str,
    bam_file: str,
    output_file: str,
    *,
    min_coverage: int = 5,
    min_splice_support: int = 3,
    max_utr_extension: int = 1000,
    enable_novel_gene_detection: bool = False,
    min_novel_gene_coverage: int = 10,
    min_novel_gene_length: int = 300,
    min_exon_length: int = 50,
    validate_splice_sites: bool = True,
    threads: int = None,
) -> dict:
    """
    Convenience function for annotation refinement with keyword arguments.
    
    Args:
        fasta_file: Path to genome FASTA file
        gff3_file: Path to input GFF3 annotations
        bam_file: Path to RNA-seq BAM alignments
        output_file: Path for refined GFF3 output
        min_coverage: Minimum coverage for UTR extensions
        min_splice_support: Minimum reads supporting splice junctions
        max_utr_extension: Maximum UTR extension length (bp)
        enable_novel_gene_detection: Enable novel gene discovery
        min_novel_gene_coverage: Minimum coverage for novel genes
        min_novel_gene_length: Minimum length for novel genes (bp)
        min_exon_length: Minimum exon length (bp)
        validate_splice_sites: Validate canonical splice sites
        threads: Number of threads (None for auto-detect, uses custom thread pool)
        
    Returns:
        Dictionary with refinement statistics and results
        
    Example:
        >>> result = annorefine.refine(
        ...     fasta_file="genome.fasta",
        ...     gff3_file="genes.gff3",
        ...     bam_file="rna_seq.bam", 
        ...     output_file="refined.gff3",
        ...     enable_novel_gene_detection=True,
        ...     threads=8
        ... )
        >>> print(f"Processed {result['genes_processed']} genes")
    """
    config = RefinementConfig(
        min_coverage=min_coverage,
        min_splice_support=min_splice_support,
        max_utr_extension=max_utr_extension,
        enable_novel_gene_detection=enable_novel_gene_detection,
        min_novel_gene_coverage=min_novel_gene_coverage,
        min_novel_gene_length=min_novel_gene_length,
        min_exon_length=min_exon_length,
        validate_splice_sites=validate_splice_sites,
    )
    
    return refine_annotations(
        fasta_file=fasta_file,
        gff3_file=gff3_file,
        bam_file=bam_file,
        output_file=output_file,
        config=config,
        threads=threads,
    )
