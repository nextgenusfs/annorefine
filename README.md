# Annorefine

[![CI](https://github.com/nextgenusfs/annorefine/workflows/CI/badge.svg)](https://github.com/nextgenusfs/annorefine/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Genome annotation refinement toolkit using RNA-seq data**

Annorefine is a high-performance toolkit for refining genome annotations using RNA-seq evidence. Built in Rust for maximum performance and reliability.

## Table of Contents

- [Features](#features)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Subcommands](#subcommands)
- [Input Requirements](#input-requirements)
- [Algorithm Overview](#algorithm-overview)
- [Performance & Scalability](#performance--scalability)
- [Integration Examples](#integration-examples)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)
- [License](#license)

## Features

### 🔧 **annorefine utrs** - RNA-seq-based Annotation Refinement
- 🧬 **UTR Extension & Trimming**: Intelligently extends and trims UTRs based on RNA-seq coverage
- 🔀 **Splice Site Refinement**: Adjusts splice sites using junction evidence from alignments
- 🆕 **Novel Gene Detection**: Optionally identifies new genes from RNA-seq evidence
- 🧭 **Comprehensive Strand Support**: Handles all RNA-seq library types (F, R, U, RF, FR, UU)
- ⚡ **High Performance**: Multi-threaded processing with Rust backend
- ✅ **Robust Validation**: Ensures all refined models remain structurally valid

### 🎯 **annorefine bam2hints** - BAM to Augustus Hints Conversion
- 🔍 **Splice Junction Detection**: Extracts intron hints from RNA-seq alignments
- 📊 **Exon/Exonpart Hints**: Generates exon and exonpart hints for gene prediction
- 🧬 **Splice Site Hints**: Identifies donor and acceptor splice sites
- 📈 **Multiplicity Support**: Counts supporting reads for each hint
- 🎛️ **Configurable Filtering**: Adjustable parameters for gap lengths, coverage, and quality
- 📝 **Augustus Compatible**: Outputs hints in standard Augustus GFF format



## Installation

### Python Package (Recommended)

Install from PyPI with pip:

```bash
pip install annorefine
```

**Supported Python versions:** 3.9, 3.10, 3.11, 3.12, 3.13
**Supported platforms:** Linux (x86_64), macOS (Intel & Apple Silicon)

### Standalone Binary

Download pre-built binaries from [GitHub Releases](https://github.com/nextgenusfs/annorefine/releases):

```bash
# Linux
wget https://github.com/nextgenusfs/annorefine/releases/download/v2025.9.18/annorefine-linux-x86_64
chmod +x annorefine-linux-x86_64

# macOS
wget https://github.com/nextgenusfs/annorefine/releases/download/v2025.9.18/annorefine-macos-arm64
chmod +x annorefine-macos-arm64
```

### Build from Source

```bash
git clone https://github.com/nextgenusfs/annorefine.git
cd annorefine
cargo build --release
# Binary will be at target/release/annorefine
```

## Quick Start

### Python API

```python
import annorefine

# Simple refinement with default settings
result = annorefine.refine(
    fasta_file="genome.fasta",
    gff3_file="annotations.gff3",
    bam_file="alignments.bam",
    output_file="refined.gff3"
)

print(f"Processed {result['genes_processed']} genes")
print(f"Novel genes detected: {result['novel_genes_detected']}")
```

### Command Line

AnnoRefine provides annotation refinement using RNA-seq evidence:

#### 🔧 Refine Existing Annotations

```bash
# Basic refinement with auto-detected library type
annorefine utrs \
    --fasta genome.fa \
    --gff3 annotations.gff3 \
    --bam alignments.bam \
    --output refined_annotations.gff3

# Specify library type for better accuracy
annorefine utrs \
    --fasta genome.fa \
    --gff3 annotations.gff3 \
    --bam alignments.bam \
    --output refined_annotations.gff3 \
    --stranded RF \
    --detect-novel-genes \
    --verbose

#### 🎯 Generate Augustus Hints

```bash
# Convert BAM alignments to Augustus hints
annorefine bam2hints \
    --in alignments.bam \
    --out hints.gff \
    --priority 4 \
    --source E

# Generate comprehensive hints for gene prediction
annorefine bam2hints \
    --in alignments.bam \
    --out hints.gff \
    --exonhints \
    --ssOn \
    --maxcoverage 1000 \
    --verbose
```



## Subcommands

### annorefine utrs

The `utrs` subcommand improves existing gene annotations using RNA-seq evidence:

**Key Features:**
- **UTR Extension & Trimming**: Extends or trims UTRs based on RNA-seq coverage patterns
- **Splice Site Refinement**: Adjusts intron boundaries using splice junction evidence
- **Novel Gene Detection**: Optionally discovers new genes from RNA-seq data
- **Strand-Aware Processing**: Supports all RNA-seq library types (F, R, U, RF, FR, UU)

**Library Type Support:**
- `auto` - Auto-detect library type (default)
- `F` - Single-end forward stranded
- `R` - Single-end reverse stranded
- `U` - Single-end unstranded
- `RF` - Paired-end reverse/forward
- `FR` - Paired-end forward/reverse
- `UU` - Paired-end unstranded

### annorefine bam2hints

The `bam2hints` subcommand converts BAM alignments into hints for Augustus gene prediction:

**Key Features:**
- **Intron Hint Generation**: Extracts splice junctions from RNA-seq alignments
- **Exon/Exonpart Hints**: Generates hints for exon boundaries and internal exons
- **Splice Site Detection**: Identifies donor and acceptor splice sites
- **Multiplicity Counting**: Aggregates identical hints with support counts
- **Coverage Filtering**: Limits hints in high-coverage regions to prevent overload
- **Augustus Compatibility**: Outputs standard Augustus GFF format

**Basic Usage:**
```bash
# Generate intron hints only (default)
annorefine bam2hints --in alignments.bam --out hints.gff

# Generate all hint types
annorefine bam2hints --in alignments.bam --out hints.gff --exonhints --ssOn

# Custom parameters
annorefine bam2hints \
    --in alignments.bam \
    --out hints.gff \
    --priority 5 \
    --source RNA \
    --minintronlen 50 \
    --maxintronlen 500000 \
    --maxcoverage 1000
```

**Important Parameters:**
- `--priority`: Priority level for hints (default: 4)
- `--source`: Source identifier for hints (default: E)
- `--minintronlen`: Minimum intron length (default: 32)
- `--maxintronlen`: Maximum intron length (default: 350000)
- `--maxcoverage`: Maximum hints per position (default: 0 = unlimited)
- `--exonhints`: Enable exon and exonpart hints
- `--ssOn`: Enable splice site hints



## Python API Reference

### Main Functions

#### `annorefine.refine()`

Convenience function with keyword arguments for all parameters:

```python
result = annorefine.refine(
    fasta_file="genome.fasta",
    gff3_file="annotations.gff3",
    bam_file="alignments.bam",
    output_file="refined.gff3",
    # Optional parameters with defaults:
    min_coverage=5,                      # Minimum coverage for UTR extension
    min_splice_support=3,                # Minimum reads supporting splice junctions
    max_utr_extension=1000,              # Maximum UTR extension length (bp)
    enable_novel_gene_detection=False,   # Enable novel gene discovery
    min_novel_gene_coverage=10,          # Minimum coverage for novel genes
    min_novel_gene_length=300,           # Minimum length for novel genes (bp)
    min_exon_length=50,                  # Minimum exon length (bp)
    validate_splice_sites=True,          # Validate canonical splice sites
    threads=None                         # Number of threads (None = auto-detect)
)
```

#### `annorefine.refine_annotations()`

Lower-level function using configuration object:

```python
# Create configuration
config = annorefine.RefinementConfig(
    min_coverage=10,
    enable_novel_gene_detection=True,
    validate_splice_sites=True
)

# Run refinement
result = annorefine.refine_annotations(
    fasta_file="genome.fasta",
    gff3_file="annotations.gff3",
    bam_file="alignments.bam",
    output_file="refined.gff3",
    config=config,
    threads=8
)
```

### Configuration Options

#### `RefinementConfig` Parameters

```python
config = annorefine.RefinementConfig(
    min_coverage=5,                      # int: Minimum coverage threshold for UTR extension
    min_splice_support=3,                # int: Minimum supporting reads for splice junctions
    max_utr_extension=1000,              # int: Maximum UTR extension length in base pairs
    enable_novel_gene_detection=False,   # bool: Enable discovery of novel genes
    min_novel_gene_coverage=10,          # int: Minimum coverage for novel gene detection
    min_novel_gene_length=300,           # int: Minimum length for novel genes (bp)
    min_exon_length=50,                  # int: Minimum exon length (bp)
    validate_splice_sites=True           # bool: Validate canonical splice sites (GT-AG, GC-AG, AT-AC)
)
```

### Return Values

Both functions return a dictionary with refinement statistics:

```python
{
    'genes_processed': 1250,              # Number of genes processed
    'genes_failed': 3,                    # Number of genes that failed processing
    'transcripts_with_structure_changes': 45,  # Transcripts with modified exon/intron structure
    'transcripts_with_5utr_extension': 234,     # Transcripts with 5' UTR extensions
    'transcripts_with_3utr_extension': 456,     # Transcripts with 3' UTR extensions
    'novel_genes_detected': 12,           # Number of novel genes discovered
    'output_file': 'refined.gff3'        # Path to output file
}
```

### Utility Functions

```python
# Get version information
version = annorefine.version()
print(f"AnnoRefine version: {version}")

# Get current thread count
threads = annorefine.current_num_threads()
print(f"Using {threads} threads")

# Test interrupt handling (for development)
annorefine.test_interruptible_operation(duration_seconds=5)
```

### Error Handling

```python
try:
    result = annorefine.refine(
        fasta_file="genome.fasta",
        gff3_file="annotations.gff3",
        bam_file="alignments.bam",
        output_file="refined.gff3"
    )
except FileNotFoundError as e:
    print(f"Input file not found: {e}")
except Exception as e:
    print(f"Refinement failed: {e}")
```

## Command Line Reference

### Basic Usage

```bash
annorefine --fasta genome.fasta --gff3 annotations.gff3 --bam alignments.bam --output refined.gff3
```

### All Options

```bash
annorefine [OPTIONS] --fasta <FILE> --gff3 <FILE> --bam <FILE> --output <FILE>
```

#### Required Arguments
- `-f, --fasta <FILE>`: Input genome FASTA file
- `-g, --gff3 <FILE>`: Input GFF3 annotation file
- `-b, --bam <FILE>`: Input BAM alignment file (RNA-seq, must be indexed)
- `-o, --output <FILE>`: Output refined GFF3 file

#### Coverage & Extension Options
- `--min-coverage <N>`: Minimum coverage threshold for UTR extension (default: 5)
- `--min-splice-support <N>`: Minimum supporting reads for splice junctions (default: 3)
- `--max-utr-extension <N>`: Maximum UTR extension length in bp (default: 1000)

#### Novel Gene Detection
- `--detect-novel-genes`: Enable novel gene detection from RNA-seq evidence
- `--min-novel-coverage <N>`: Minimum coverage for novel genes (default: 10)
- `--min-novel-length <N>`: Minimum length for novel genes in bp (default: 300)

#### Quality Control
- `--min-exon-length <N>`: Minimum exon length in bp (default: 50)
- `--no-splice-validation`: Disable canonical splice site validation

#### Performance & Output
- `-t, --threads <N>`: Number of threads for parallel processing (default: auto-detect)
- `-v, --verbose`: Enable verbose output (shows warnings and debug info)
- `--log-file <FILE>`: Write detailed log to file

#### Information
- `-h, --help`: Print help information
- `-V, --version`: Print version information

### Examples

```bash
# Basic refinement
annorefine -f genome.fa -g genes.gff3 -b rna_seq.bam -o refined.gff3

# High-sensitivity novel gene detection
annorefine -f genome.fa -g genes.gff3 -b rna_seq.bam -o refined.gff3 \
    --detect-novel-genes --min-novel-coverage 5 --min-coverage 3

# Performance optimization
annorefine -f genome.fa -g genes.gff3 -b rna_seq.bam -o refined.gff3 \
    --threads 16 --log-file refinement.log

# Conservative refinement (higher thresholds)
annorefine -f genome.fa -g genes.gff3 -b rna_seq.bam -o refined.gff3 \
    --min-coverage 10 --min-splice-support 5 --max-utr-extension 500
```

## Input Requirements

### Required Files

1. **Genome FASTA file** (`.fasta`, `.fa`, `.fna`)
   - Reference genome sequences
   - Can be compressed (`.gz`)

2. **GFF3 annotation file** (`.gff3`, `.gff`)
   - Must contain `gene`, `mRNA`, `exon`, and `CDS` features
   - Follows GFF3 specification
   - Can be compressed (`.gz`)

3. **BAM alignment file** (`.bam`)
   - RNA-seq alignments to the reference genome
   - **Must be sorted and indexed** (`.bam.bai` file required)
   - Splice-aware alignment recommended (STAR, HISAT2, etc.)

### File Preparation

#### Creating BAM Files

```bash
# Option 1: Using STAR (recommended for RNA-seq)
STAR --runMode genomeGenerate --genomeDir genome_index --genomeFastaFiles genome.fasta
STAR --genomeDir genome_index --readFilesIn reads_R1.fastq reads_R2.fastq \
     --outSAMtype BAM SortedByCoordinate --outFileNamePrefix sample_

# Index the BAM file
samtools index sample_Aligned.sortedByCoord.out.bam

# Option 2: Using HISAT2
hisat2-build genome.fasta genome_index
hisat2 -x genome_index -1 reads_R1.fastq -2 reads_R2.fastq | \
    samtools sort -o alignments.bam
samtools index alignments.bam
```

#### Validating Input Files

```python
import annorefine

# Check if files are accessible and properly formatted
try:
    result = annorefine.refine(
        fasta_file="genome.fasta",
        gff3_file="annotations.gff3",
        bam_file="alignments.bam",
        output_file="test_output.gff3"
    )
    print("✅ All input files are valid")
except Exception as e:
    print(f"❌ Input validation failed: {e}")
```

## Algorithm Overview

AnnoRefine uses a multi-step approach to refine genome annotations:

### 1. Input Parsing & Validation
- Load genome sequences from FASTA
- Parse gene models from GFF3 (genes → transcripts → exons → CDS)
- Index BAM file for efficient region-based queries
- Validate file formats and cross-references

### 2. Evidence Extraction
For each gene region:
- Extract RNA-seq coverage profiles
- Identify splice junctions with read support
- Calculate coverage statistics and junction confidence

### 3. Gene Model Refinement
- **UTR Extension**: Extend 5' and 3' UTRs based on continuous coverage
- **Splice Site Adjustment**: Refine intron/exon boundaries using well-supported junctions
- **Structure Validation**: Ensure all changes maintain valid gene model structure
- **CDS Preservation**: Maintain coding sequence integrity

### 4. Novel Gene Detection (Optional)
- Identify regions with RNA-seq coverage but no existing annotations
- Predict gene structures from splice junction patterns
- Filter candidates by coverage, length, and splice site quality

### 5. Output Generation
- Write refined annotations in GFF3 format
- Preserve original feature attributes and metadata
- Add refinement statistics and provenance information

## Performance & Scalability

- **Multi-threaded**: Parallel processing of gene regions
- **Memory efficient**: Streaming BAM processing, minimal memory footprint
- **Fast I/O**: Optimized file parsing and writing
- **Scalable**: Handles mammalian-sized genomes efficiently

**Typical performance:**
- Human genome (~20K genes): 10-30 minutes on 8 cores
- Plant genome (~30K genes): 15-45 minutes on 8 cores
- Memory usage: 2-8 GB depending on genome size

## Integration Examples

### Nextflow Pipeline

```nextflow
process ANNOREFINE {
    conda 'pip::annorefine'

    input:
    path genome_fasta
    path annotations_gff3
    path alignments_bam
    path alignments_bai

    output:
    path "refined.gff3"

    script:
    """
    annorefine \\
        --fasta ${genome_fasta} \\
        --gff3 ${annotations_gff3} \\
        --bam ${alignments_bam} \\
        --output refined.gff3 \\
        --threads ${task.cpus} \\
        --detect-novel-genes
    """
}
```

### Snakemake Rule

```python
rule annorefine:
    input:
        fasta="genome.fasta",
        gff3="annotations.gff3",
        bam="alignments.bam",
        bai="alignments.bam.bai"
    output:
        "refined_annotations.gff3"
    conda:
        "envs/annorefine.yaml"  # pip: annorefine
    threads: 8
    shell:
        """
        annorefine --fasta {input.fasta} --gff3 {input.gff3} \\
                   --bam {input.bam} --output {output} \\
                   --threads {threads} --detect-novel-genes
        """
```

## Troubleshooting

### Common Issues

**"BAM file not indexed"**
```bash
samtools index alignments.bam
```

**"No splice junctions found"**
- Ensure BAM contains splice-aware alignments (use STAR, HISAT2, not BWA)
- Check that RNA-seq reads span introns

**"Low refinement rate"**
- Increase RNA-seq depth (>50M reads recommended)
- Lower `--min-coverage` threshold
- Check RNA-seq quality and mapping rate

**"Memory usage too high"**
- Reduce `--threads` parameter
- Process smaller genomic regions separately

### Getting Help

- 📖 **Documentation**: [GitHub Wiki](https://github.com/nextgenusfs/annorefine/wiki)
- 🐛 **Bug Reports**: [GitHub Issues](https://github.com/nextgenusfs/annorefine/issues)
- 💬 **Discussions**: [GitHub Discussions](https://github.com/nextgenusfs/annorefine/discussions)

## Citation

If you use AnnoRefine in your research, please cite:

```
Palmer, J. (2025). AnnoRefine: High-performance genome annotation refinement using RNA-seq data.
GitHub: https://github.com/nextgenusfs/annorefine
```

## License

MIT License - see [LICENSE](LICENSE) file for details.

