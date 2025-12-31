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
    bam2hints_convert,
    join_hints,
    filter_hints,
    version,
    current_num_threads,
    test_interruptible_operation,
    PyRefinementConfig as RefinementConfig,
    PyBam2HintsConfig as Bam2HintsConfig,
    PyGeneModel as GeneModel,
)

__version__ = version()
__author__ = "Jon Palmer"
__email__ = "nextgenusfs@gmail.com"
__description__ = "Genome annotation refinement using RNA-seq data"

__all__ = [
    "refine",
    "refine_annotations",
    "bam2hints",
    "bam2hints_convert",
    "join_hints",
    "filter_hints",
    "version",
    "current_num_threads",
    "test_interruptible_operation",
    "RefinementConfig",
    "Bam2HintsConfig",
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
    strand_bias_threshold: float = 0.65,
    max_reads_for_strand_detection: int = 10000,
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
        strand_bias_threshold: Threshold for detecting stranded RNA-seq (0.5-1.0)
        max_reads_for_strand_detection: Maximum reads to sample for strand detection
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
        strand_bias_threshold=strand_bias_threshold,
        max_reads_for_strand_detection=max_reads_for_strand_detection,
    )
    
    return refine_annotations(
        fasta_file=fasta_file,
        gff3_file=gff3_file,
        bam_file=bam_file,
        output_file=output_file,
        config=config,
        threads=threads,
    )


def bam2hints(
    bam_file: str,
    output_file: str,
    library_type: str,
    *,
    priority: int = 4,
    max_gap_len: int = 14,
    min_intron_len: int = 32,
    max_intron_len: int = 350000,
    min_end_block_len: int = 8,
    max_query_gap_len: int = 5,
    exonpart_cutoff: int = 10,
    source: str = "E",
    introns_only: bool = False,
    no_multiplicity: bool = False,
    remove_redundant: bool = False,
    max_coverage: int = 0,
    splice_sites_on: bool = False,
    truncated_splice_sites: bool = False,
    score: float = 0.0,
    max_gene_len: int = 400000,
    threads: int = None,
    contig: str = None,
    region: tuple = None,
    contig_map: dict = None,
) -> dict:
    """
    Convenience function for BAM to Augustus hints conversion with keyword arguments.

    Args:
        bam_file: Path to input BAM file (must be sorted and indexed if using contig/region filtering)
        output_file: Path for output GFF hints file
        library_type: Library strandedness specification. Options: "FR", "RF", "UU"
                     FR = paired-end forward/reverse (fr-secondstrand)
                     RF = paired-end reverse/forward (fr-firststrand)
                     UU = paired-end unstranded
        priority: Priority of hint group (default: 4)
        max_gap_len: Gaps at most this length are simply closed (default: 14)
        min_intron_len: Minimum intron length (default: 32)
        max_intron_len: Maximum intron length (default: 350000)
        min_end_block_len: Minimum length of a 'dangling' exon (default: 8)
        max_query_gap_len: Maximum length of gap in query sequence (default: 5)
        exonpart_cutoff: BP cut off of each exonpart hint at end of alignment (default: 10)
        source: Source identifier (default: "E")
        introns_only: Only retrieve intron hints (default: False)
        no_multiplicity: Do not summarize multiple identical intron hints (default: False)
        remove_redundant: Only keep the strongest hint for a region (default: False)
        max_coverage: Maximum number of hints at a given position, 0=unlimited (default: 0)
        splice_sites_on: Include splice site (dss, ass) hints (default: False)
        truncated_splice_sites: Include splice sites from truncated alignments (default: False)
        score: Score value to fill in score column (default: 0.0)
        max_gene_len: Alignments spanning more than this are ignored (default: 400000)
        threads: Number of threads to use for parallel processing (default: None, uses all available)
        contig: Filter to only process alignments on this contig (default: None, process all)
        region: Filter to only process alignments in this region as (contig, start, end) tuple (default: None)
        contig_map: Dictionary to rename contigs in output (default: None, empty dict)
                   Keys are input contig names, values are output contig names

    Returns:
        Dictionary with conversion statistics and results

    Example:
        >>> # Process entire BAM file
        >>> result = annorefine.bam2hints(
        ...     bam_file="alignments.bam",
        ...     output_file="hints.gff",
        ...     priority=5,
        ...     source="RNA",
        ...     splice_sites_on=True
        ... )
        >>> print(f"Generated {result['total_hints_generated']} hints")

        >>> # Process only a specific contig
        >>> result = annorefine.bam2hints(
        ...     bam_file="alignments.bam",
        ...     output_file="contig1_hints.gff",
        ...     contig="contig_1",
        ...     introns_only=True
        ... )

        >>> # Process only a specific region
        >>> result = annorefine.bam2hints(
        ...     bam_file="alignments.bam",
        ...     output_file="region_hints.gff",
        ...     region=("chr1", 1000, 50000)
        ... )
    """
    # Match CLI behavior: enable splice sites by default when not introns-only
    # unless explicitly disabled
    effective_splice_sites_on = splice_sites_on or not introns_only

    config = Bam2HintsConfig(
        library_type=library_type,
        priority=priority,
        max_gap_len=max_gap_len,
        min_intron_len=min_intron_len,
        max_intron_len=max_intron_len,
        min_end_block_len=min_end_block_len,
        max_query_gap_len=max_query_gap_len,
        exonpart_cutoff=exonpart_cutoff,
        source=source,
        introns_only=introns_only,
        no_multiplicity=no_multiplicity,
        remove_redundant=remove_redundant,
        max_coverage=max_coverage,
        splice_sites_on=effective_splice_sites_on,
        truncated_splice_sites=truncated_splice_sites,
        score=score,
        max_gene_len=max_gene_len,
        contig_map=contig_map if contig_map is not None else {},
    )

    return bam2hints_convert(
        bam_file=bam_file,
        output_file=output_file,
        library_type=library_type,
        config=config,
        threads=threads,
        contig=contig,
        region=region,
    )
