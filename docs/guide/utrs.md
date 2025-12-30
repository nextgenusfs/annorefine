# UTR Refinement

AnnoRefine can extend and refine UTRs (Untranslated Regions) in existing gene annotations using RNA-seq evidence.

## Overview

The `utrs` command refines gene annotations by:

- Extending 5' and 3' UTRs based on RNA-seq coverage
- Detecting novel genes from RNA-seq data
- Validating splice sites
- Improving gene model accuracy

## Command Line Usage

### Basic Usage

```bash
annorefine utrs \
    --fasta genome.fa \
    --gff3 annotations.gff3 \
    --bam alignments.bam \
    --output refined.gff3
```

### With Novel Gene Detection

```bash
annorefine utrs \
    --fasta genome.fa \
    --gff3 annotations.gff3 \
    --bam alignments.bam \
    --output refined.gff3 \
    --enable-novel-gene-detection \
    --min-novel-gene-coverage 10
```

### Advanced Options

```bash
annorefine utrs \
    --fasta genome.fa \
    --gff3 annotations.gff3 \
    --bam alignments.bam \
    --output refined.gff3 \
    --min-coverage 5 \
    --min-splice-support 3 \
    --max-utr-extension 1000 \
    --validate-splice-sites \
    --threads 8
```

## Python API Usage

### Basic Refinement

```python
import annorefine

result = annorefine.refine(
    fasta_file="genome.fa",
    gff3_file="annotations.gff3",
    bam_file="alignments.bam",
    output_file="refined.gff3"
)

print(f"Processed {result['genes_processed']} genes")
print(f"UTRs extended: {result['utrs_extended']}")
```

### With Novel Gene Detection

```python
result = annorefine.refine(
    fasta_file="genome.fa",
    gff3_file="annotations.gff3",
    bam_file="alignments.bam",
    output_file="refined.gff3",
    enable_novel_gene_detection=True,
    min_novel_gene_coverage=10,
    min_novel_gene_length=300
)

print(f"Novel genes detected: {result['novel_genes_detected']}")
```

### Advanced Configuration

```python
result = annorefine.refine(
    fasta_file="genome.fa",
    gff3_file="annotations.gff3",
    bam_file="alignments.bam",
    output_file="refined.gff3",
    min_coverage=5,
    min_splice_support=3,
    max_utr_extension=1000,
    enable_novel_gene_detection=True,
    min_novel_gene_coverage=10,
    min_novel_gene_length=300,
    min_exon_length=50,
    validate_splice_sites=True,
    strand_bias_threshold=0.65,
    threads=8
)
```

### Using Configuration Objects

```python
import annorefine

# Create a custom configuration
config = annorefine.RefinementConfig(
    min_coverage=10,
    min_splice_support=5,
    max_utr_extension=2000,
    enable_novel_gene_detection=True,
    min_novel_gene_coverage=15,
    min_novel_gene_length=500,
    validate_splice_sites=True
)

# Run refinement with custom config
result = annorefine.refine_annotations(
    fasta_file="genome.fa",
    gff3_file="annotations.gff3",
    bam_file="alignments.bam",
    output_file="refined.gff3",
    config=config,
    threads=8
)
```

## Parameters

### Coverage Parameters

- **min_coverage**: Minimum RNA-seq coverage to extend UTRs (default: 5)
- **min_splice_support**: Minimum reads supporting a splice junction (default: 3)
- **min_novel_gene_coverage**: Minimum coverage for novel gene detection (default: 10)

### Extension Parameters

- **max_utr_extension**: Maximum bases to extend UTRs (default: 1000)
- **min_exon_length**: Minimum exon length for novel genes (default: 50)
- **min_novel_gene_length**: Minimum total length for novel genes (default: 300)

### Validation Parameters

- **validate_splice_sites**: Check for canonical splice sites (GT-AG) (default: True)
- **strand_bias_threshold**: Threshold for strand detection (default: 0.65)

### Performance Parameters

- **threads**: Number of threads for parallel processing (default: all available)

## Output

The output is a refined GFF3 file with:

- Extended UTRs based on RNA-seq evidence
- Novel genes (if enabled)
- Improved gene models
- Validated splice sites

## Example Workflow

```python
import annorefine

# Step 1: Refine existing annotations
result = annorefine.refine(
    fasta_file="genome.fa",
    gff3_file="initial_annotations.gff3",
    bam_file="rnaseq.bam",
    output_file="refined_annotations.gff3",
    enable_novel_gene_detection=True,
    threads=8
)

print(f"Refinement complete:")
print(f"  Genes processed: {result['genes_processed']}")
print(f"  UTRs extended: {result['utrs_extended']}")
print(f"  Novel genes: {result['novel_genes_detected']}")

# Step 2: Generate hints for Augustus
annorefine.bam2hints(
    bam_file="rnaseq.bam",
    output_file="hints.gff",
    library_type="RF",
    threads=8
)

# Step 3: Run Augustus with refined annotations and hints
# (external command)
```

## Performance Tips

1. **Use multiple threads**: Set `threads` parameter to match CPU cores
2. **Adjust coverage thresholds**: Lower for low-coverage data, higher for deep sequencing
3. **Validate splice sites**: Enable for better quality, disable for speed
4. **Novel gene detection**: Disable if not needed to improve performance

## Next Steps

- [Python API Guide](python.md) - Complete Python API documentation
- [BAM to Hints](bam2hints.md) - Generate hints for gene prediction
- [API Reference](../api/overview.md) - Detailed API reference

