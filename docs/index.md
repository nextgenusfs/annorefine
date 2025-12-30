# AnnoRefine Documentation

**Genome annotation refinement toolkit using RNA-seq data**

AnnoRefine is a high-performance toolkit for refining genome annotations using RNA-seq evidence. It provides both command-line tools and Python bindings for integration into bioinformatics pipelines.

## Key Features

- **BAM to Hints Conversion**: Convert RNA-seq BAM alignments to Augustus/GeneMark hints format
- **UTR Refinement**: Extend and refine UTRs using RNA-seq evidence
- **Novel Gene Detection**: Identify novel genes from RNA-seq data
- **Python Bindings**: Full Python API for pipeline integration
- **High Performance**: Written in Rust with parallel processing support
- **Memory Efficient**: Optimized for large genomes and deep sequencing data

## Quick Start

### Installation

```bash
pip install annorefine
```

### Command Line Usage

Convert BAM alignments to Augustus hints:

```bash
annorefine bam2hints \
    --input alignments.bam \
    --output hints.gff \
    --stranded RF \
    --threads 4
```

Refine UTRs using RNA-seq evidence:

```bash
annorefine utrs \
    --fasta genome.fa \
    --gff3 annotations.gff3 \
    --bam alignments.bam \
    --output refined.gff3
```

### Python API Usage

```python
import annorefine

# Convert BAM to hints
result = annorefine.bam2hints(
    bam_file="alignments.bam",
    output_file="hints.gff",
    library_type="RF",
    threads=4
)

# Join hints from multiple sources
result = annorefine.join_hints(
    input_files=["bam_hints.gff", "protein_hints.gff"],
    output_file="joined_hints.gff"
)

# Refine annotations
result = annorefine.refine(
    fasta_file="genome.fa",
    gff3_file="annotations.gff3",
    bam_file="alignments.bam",
    output_file="refined.gff3",
    enable_novel_gene_detection=True
)
```

## Use Cases

### Augustus Gene Prediction

AnnoRefine generates hints in the format expected by Augustus gene predictor:

```bash
# Generate all hint types for Augustus
annorefine bam2hints -i alignments.bam -o augustus_hints.gff --stranded RF

# Run Augustus with hints
augustus --species=your_species --hintsfile=augustus_hints.gff genome.fa
```

### GeneMark-ETP Integration

Generate intron-only hints for GeneMark:

```python
import annorefine

# Generate intron hints for GeneMark
result = annorefine.join_hints(
    input_files=["bam_hints.gff", "protein_hints.gff"],
    output_file="genemark_hints.gff",
    introns_only=True
)
```

### Funannotate2 Pipeline

AnnoRefine is designed for integration with funannotate2:

```python
import annorefine

# Per-contig processing for parallel execution
for contig in contigs:
    # Generate hints for this contig
    annorefine.bam2hints(
        bam_file="alignments.bam",
        output_file=f"{contig}.hints.gff",
        library_type="RF",
        contig=contig
    )
    
    # Join with other evidence
    annorefine.join_hints(
        input_files=[
            f"{contig}.bam_hints.gff",
            f"{contig}.protein_hints.gff"
        ],
        output_file=f"{contig}.augustus_hints.gff"
    )
```

## Documentation

- **[Installation Guide](guide/installation.md)** - Installation instructions and requirements
- **[BAM to Hints](guide/bam2hints.md)** - Convert RNA-seq alignments to hints
- **[UTR Refinement](guide/utrs.md)** - Refine UTRs using RNA-seq evidence
- **[Python API](guide/python.md)** - Complete Python API documentation
- **[API Reference](api/overview.md)** - Detailed API reference

## Performance

AnnoRefine is optimized for performance:

- Written in Rust for maximum speed
- Parallel processing with configurable thread count
- Memory-efficient streaming of BAM files
- Optimized for large genomes (tested on multi-GB genomes)

## Links

- [GitHub Repository](https://github.com/nextgenusfs/annorefine)
- [PyPI Package](https://pypi.org/project/annorefine/)
- [Issue Tracker](https://github.com/nextgenusfs/annorefine/issues)
- [License](https://github.com/nextgenusfs/annorefine/blob/main/LICENSE)

## Citation

If you use AnnoRefine in your research, please cite:

```
Palmer, J. (2025). AnnoRefine: Genome annotation refinement toolkit using RNA-seq data.
https://github.com/nextgenusfs/annorefine
```

