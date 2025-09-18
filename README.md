# AnnoRefine

[![CI](https://github.com/nextgenusfs/annorefine/workflows/CI/badge.svg)](https://github.com/nextgenusfs/annorefine/actions)
[![PyPI version](https://badge.fury.io/py/annorefine.svg)](https://badge.fury.io/py/annorefine)
[![Python versions](https://img.shields.io/pypi/pyversions/annorefine.svg)](https://pypi.org/project/annorefine/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Genome annotation refinement using RNA-seq data**

AnnoRefine is a high-performance tool for refining genome annotations using RNA-seq evidence. Available as both a Rust command-line tool and Python library.

## Installation

### Python Package (Recommended)

```bash
pip install annorefine
```

### Rust Binary

Download pre-built binaries from [GitHub Releases](https://github.com/nextgenusfs/annorefine/releases) or build from source:

```bash
cargo install annorefine
```

## Quick Start

### Python API

```python
import annorefine

# Basic usage
result = annorefine.refine(
    fasta_file="genome.fasta",
    gff3_file="annotations.gff3",
    bam_file="alignments.bam",
    output_file="refined.gff3"
)

print(f"Processed {result['genes_processed']} genes")
print(f"Found {result['novel_genes_detected']} novel genes")

# Advanced configuration with config object
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

# Or use keyword arguments directly
result = annorefine.refine(
    fasta_file="genome.fasta",
    gff3_file="annotations.gff3",
    bam_file="alignments.bam",
    output_file="refined.gff3",
    min_coverage=10,
    enable_novel_gene_detection=True,
    validate_splice_sites=True,
    threads=8
)
```

### Command Line

```bash
# Basic refinement
annorefine -f genome.fasta -g annotations.gff3 -b alignments.bam -o refined.gff3

# With novel gene detection
annorefine -f genome.fasta -g annotations.gff3 -b alignments.bam -o refined.gff3 \
    --detect-novel-genes --min-novel-coverage 10
```

## Overview

AnnoRefine takes as input:
- A genome FASTA file
- A genome annotation file in GFF3 format  
- A genome BAM alignment file of RNA-seq data

The tool parses the GFF3 gene models and uses the transcriptomic alignments in the BAM file to refine gene model predictions. Refinements can include:
- Adding 5' UTR extensions
- Adding 3' UTR extensions  
- Changing intron/exon structure based on alignments (while maintaining valid gene models)

## Installation

```bash
git clone <repository>
cd annorefine
cargo build --release
```

## Usage

```bash
annorefine --fasta genome.fasta --gff3 annotations.gff3 --bam alignments.bam --output refined_annotations.gff3
```

### Options

- `-f, --fasta <FILE>`: Input genome FASTA file
- `-g, --gff3 <FILE>`: Input GFF3 annotation file
- `-b, --bam <FILE>`: Input BAM alignment file (RNA-seq)
- `-o, --output <FILE>`: Output refined GFF3 file
- `--min-coverage <N>`: Minimum coverage threshold for UTR extension (default: 5)
- `--min-splice-support <N>`: Minimum number of supporting reads for splice junction (default: 3)
- `--max-utr-extension <N>`: Maximum UTR extension length in bp (default: 1000)
- `-v, --verbose`: Verbose output
- `-h, --help`: Print help
- `-V, --version`: Print version

## Testing

### Test Data

Sample test files are provided in the `test_data/` directory:
- `test_genome.fasta`: A simple test genome sequence
- `test_annotations.gff3`: Sample gene annotations

### Creating Test BAM Files

To test the full pipeline, you'll need a BAM file. You can create one using tools like:

1. **Using minimap2 and samtools:**
   ```bash
   # Align RNA-seq reads to genome
   minimap2 -ax splice test_data/test_genome.fasta reads.fastq > alignments.sam
   
   # Convert to BAM and index
   samtools view -bS alignments.sam > alignments.bam
   samtools sort alignments.bam > sorted_alignments.bam
   samtools index sorted_alignments.bam
   ```

2. **Using STAR:**
   ```bash
   # Build genome index
   STAR --runMode genomeGenerate --genomeDir genome_index --genomeFastaFiles test_data/test_genome.fasta
   
   # Align reads
   STAR --genomeDir genome_index --readFilesIn reads.fastq --outFileNamePrefix aligned_
   ```

### Running Tests

```bash
# Test with sample data (requires a BAM file)
cargo run -- -f test_data/test_genome.fasta -g test_data/test_annotations.gff3 -b your_alignments.bam -o refined_output.gff3 -v

# Test FASTA parsing only
cargo test test_fasta

# Test GFF3 parsing only  
cargo test test_gff3

# Run all tests
cargo test
```

## Architecture

The application is structured into several modules:

- `types.rs`: Core data structures for genome sequences, annotations, and alignments
- `fasta.rs`: FASTA file parsing and genome sequence handling
- `gff3.rs`: GFF3 file parsing and gene model construction
- `bam.rs`: BAM file reading and RNA-seq alignment processing
- `refinement.rs`: Gene model refinement engine using RNA-seq evidence
- `output.rs`: Output generation for refined gene models in GFF3 format

## Algorithm

1. **Parse Input Files**: Load genome sequences, gene models, and RNA-seq alignments
2. **Extract Evidence**: For each gene region, extract splice junctions and coverage from RNA-seq data
3. **Refine Structure**: Adjust intron/exon boundaries based on well-supported splice junctions
4. **Extend UTRs**: Extend 5' and 3' UTRs based on continuous coverage evidence
5. **Validate Models**: Ensure all refined gene models remain structurally valid
6. **Output Results**: Write refined annotations in GFF3 format

## Requirements

- Rust 1.70 or later
- Input files:
  - Genome FASTA file
  - GFF3 annotation file with gene/mRNA/exon/CDS features
  - Indexed BAM file with RNA-seq alignments

## License

