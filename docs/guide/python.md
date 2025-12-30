# Python API Guide

AnnoRefine provides a comprehensive Python API for integration into bioinformatics pipelines.

## Installation

```bash
pip install annorefine
```

## Quick Start

```python
import annorefine

# Check version
print(annorefine.version())

# Convert BAM to hints
result = annorefine.bam2hints(
    bam_file="alignments.bam",
    output_file="hints.gff",
    library_type="RF"
)

# Refine annotations
result = annorefine.refine(
    fasta_file="genome.fa",
    gff3_file="annotations.gff3",
    bam_file="alignments.bam",
    output_file="refined.gff3"
)
```

## Core Functions

### bam2hints()

Convert BAM alignments to Augustus/GeneMark hints.

```python
result = annorefine.bam2hints(
    bam_file: str,              # Input BAM file (required)
    output_file: str,           # Output GFF file (required)
    library_type: str,          # "FR", "RF", or "UU" (required)
    priority: int = 4,          # Hint priority
    source: str = "E",          # Source identifier
    introns_only: bool = False, # Only output intron hints
    threads: int = None,        # Number of threads (None = all)
    contig: str = None,         # Filter to specific contig
    region: tuple = None        # Filter to region (contig, start, end)
)
```

**Returns:** Dictionary with statistics

```python
{
    'total_hints_generated': 117988,
    'intron_hints': 20956,
    'exon_hints': 8309,
    'exonpart_hints': 48119,
    'dss_hints': 20302,
    'ass_hints': 20302
}
```

### join_hints()

Join and merge hints from multiple sources.

```python
result = annorefine.join_hints(
    input_files: List[str],     # List of input GFF files (required)
    output_file: str,           # Output GFF file (required)
    introns_only: bool = False  # Filter to introns only
)
```

**Returns:** Dictionary with statistics

```python
{
    'input_files': 2,
    'total_input_hints': 235976,
    'output_hints': 116772,
    'output_file': 'joined.gff'
}
```

### filter_hints()

Filter hints by type, multiplicity, or contig.

```python
result = annorefine.filter_hints(
    input_file: str,            # Input GFF file (required)
    output_file: str,           # Output GFF file (required)
    hint_types: List[str] = None,  # Filter by hint types
    min_mult: int = None,       # Minimum multiplicity
    contig: str = None          # Filter to specific contig
)
```

### refine()

Refine gene annotations using RNA-seq evidence.

```python
result = annorefine.refine(
    fasta_file: str,                        # Genome FASTA (required)
    gff3_file: str,                         # Input GFF3 (required)
    bam_file: str,                          # RNA-seq BAM (required)
    output_file: str,                       # Output GFF3 (required)
    min_coverage: int = 5,                  # Min coverage for UTR extension
    min_splice_support: int = 3,            # Min reads for splice junction
    max_utr_extension: int = 1000,          # Max UTR extension length
    enable_novel_gene_detection: bool = False,  # Detect novel genes
    min_novel_gene_coverage: int = 10,      # Min coverage for novel genes
    min_novel_gene_length: int = 300,       # Min length for novel genes
    min_exon_length: int = 50,              # Min exon length
    validate_splice_sites: bool = True,     # Validate splice sites
    strand_bias_threshold: float = 0.65,    # Strand detection threshold
    threads: int = None                     # Number of threads
)
```

## Funannotate2 Integration

AnnoRefine is designed for seamless integration with funannotate2 pipelines.

### Per-Contig Parallel Processing

```python
import annorefine
from multiprocessing import Pool

def process_contig(contig):
    """Process a single contig in parallel"""
    # Generate BAM hints
    annorefine.bam2hints(
        bam_file="alignments.bam",
        output_file=f"{contig}.bam_hints.gff",
        library_type="RF",
        contig=contig,
        threads=2
    )
    
    # Join with other evidence
    annorefine.join_hints(
        input_files=[
            f"{contig}.bam_hints.gff",
            f"{contig}.protein_hints.gff"
        ],
        output_file=f"{contig}.augustus_hints.gff"
    )
    
    # Create GeneMark hints (introns only)
    annorefine.join_hints(
        input_files=[
            f"{contig}.bam_hints.gff",
            f"{contig}.protein_hints.gff"
        ],
        output_file=f"{contig}.genemark_hints.gff",
        introns_only=True
    )
    
    return contig

# Process all contigs in parallel
contigs = ["chr1", "chr2", "chr3"]
with Pool(processes=8) as pool:
    results = pool.map(process_contig, contigs)
```

### Complete Workflow

```python
import annorefine

# Step 1: Generate hints from RNA-seq
bam_result = annorefine.bam2hints(
    bam_file="rnaseq.bam",
    output_file="bam_hints.gff",
    library_type="RF",
    threads=8
)

# Step 2: Join hints from multiple sources
augustus_result = annorefine.join_hints(
    input_files=[
        "bam_hints.gff",
        "protein_hints.gff",
        "transcript_hints.gff"
    ],
    output_file="augustus_hints.gff"
)

# Step 3: Create GeneMark hints
genemark_result = annorefine.join_hints(
    input_files=[
        "bam_hints.gff",
        "protein_hints.gff"
    ],
    output_file="genemark_hints.gff",
    introns_only=True
)

# Step 4: Refine existing annotations
refine_result = annorefine.refine(
    fasta_file="genome.fa",
    gff3_file="initial_annotations.gff3",
    bam_file="rnaseq.bam",
    output_file="refined_annotations.gff3",
    enable_novel_gene_detection=True,
    threads=8
)

print(f"Augustus hints: {augustus_result['output_hints']}")
print(f"GeneMark hints: {genemark_result['output_hints']}")
print(f"Genes refined: {refine_result['genes_processed']}")
print(f"Novel genes: {refine_result['novel_genes_detected']}")
```

## Configuration Objects

For advanced usage, use configuration objects:

```python
import annorefine

# Create custom configuration
config = annorefine.RefinementConfig(
    min_coverage=10,
    min_splice_support=5,
    max_utr_extension=2000,
    enable_novel_gene_detection=True,
    validate_splice_sites=True
)

# Use with refine_annotations
result = annorefine.refine_annotations(
    fasta_file="genome.fa",
    gff3_file="annotations.gff3",
    bam_file="alignments.bam",
    output_file="refined.gff3",
    config=config,
    threads=8
)
```

## Error Handling

```python
import annorefine

try:
    result = annorefine.bam2hints(
        bam_file="alignments.bam",
        output_file="hints.gff",
        library_type="RF"
    )
except Exception as e:
    print(f"Error: {e}")
```

## Performance Tips

1. **Use all available cores**: Set `threads=None` to use all CPU cores
2. **Process by contig**: For large genomes, process contigs in parallel
3. **Index BAM files**: Ensure BAM files have `.bai` index for region filtering
4. **Batch processing**: Process multiple samples in parallel using multiprocessing

## Next Steps

- [API Reference](../api/overview.md) - Detailed API documentation
- [BAM to Hints Guide](bam2hints.md) - Learn about hint generation
- [UTR Refinement Guide](utrs.md) - Learn about annotation refinement

