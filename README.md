# AnnoRefine

[![CI](https://github.com/nextgenusfs/annorefine/workflows/CI/badge.svg)](https://github.com/nextgenusfs/annorefine/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![PyPI](https://img.shields.io/pypi/v/annorefine.svg)](https://pypi.org/project/annorefine/)
[![Documentation](https://img.shields.io/badge/docs-mkdocs-blue.svg)](https://nextgenusfs.github.io/annorefine/)

**High-performance genome annotation refinement toolkit using RNA-seq data**

AnnoRefine is a Rust-based toolkit for refining genome annotations and generating gene prediction hints from RNA-seq evidence. It provides both command-line tools and Python bindings for seamless integration into bioinformatics pipelines.

📖 **[Full Documentation](https://nextgenusfs.github.io/annorefine/)**

## Features

- 🔧 **UTR Refinement** - Extend and trim UTRs based on RNA-seq coverage
- 🔀 **Splice Site Refinement** - Adjust intron boundaries using junction evidence
- 🆕 **Novel Gene Detection** - Discover new genes from RNA-seq data
- 🎯 **Hint Generation** - Convert BAM alignments to Augustus/GeneMark hints
- 📊 **Hint Processing** - Join and filter hints from multiple sources
- ⚡ **High Performance** - Multi-threaded Rust implementation
- 🐍 **Python Bindings** - Easy integration into Python workflows
- 🧭 **Strand-Aware** - Supports all RNA-seq library types (FR, RF, UU)



## Installation

**Python Package (Recommended):**
```bash
pip install annorefine
```

**Standalone Binary:**
Download from [GitHub Releases](https://github.com/nextgenusfs/annorefine/releases)

**Build from Source:**
```bash
git clone https://github.com/nextgenusfs/annorefine.git
cd annorefine
cargo build --release
```

See the [Installation Guide](https://nextgenusfs.github.io/annorefine/guide/installation/) for detailed instructions.

## Quick Start

**Python API:**
```python
import annorefine

# Refine annotations
result = annorefine.refine(
    fasta_file="genome.fa",
    gff3_file="annotations.gff3",
    bam_file="alignments.bam",
    output_file="refined.gff3"
)

# Generate hints for gene prediction
result = annorefine.bam2hints(
    bam_file="alignments.bam",
    output_file="hints.gff",
    library_type="RF"
)

# Join hints from multiple sources
result = annorefine.join_hints(
    input_files=["bam_hints.gff", "protein_hints.gff"],
    output_file="joined_hints.gff"
)
```

**Command Line:**
```bash
# Refine annotations
annorefine utrs \
    --fasta genome.fa \
    --gff3 annotations.gff3 \
    --bam alignments.bam \
    --output refined.gff3

# Generate hints
annorefine bam2hints \
    --in alignments.bam \
    --out hints.gff \
    --stranded RF

# Join hints
annorefine join-hints \
    --input bam_hints.gff protein_hints.gff \
    --output joined_hints.gff
```

See the [User Guide](https://nextgenusfs.github.io/annorefine/guide/bam2hints/) for more examples.



## Use Cases

- **Annotation Refinement** - Improve existing gene models with RNA-seq evidence
- **Augustus Gene Prediction** - Generate hints for ab initio gene prediction
- **GeneMark-ETP** - Create intron-only hints for GeneMark
- **funannotate2 Integration** - Seamless integration with gene prediction pipelines

## Documentation

- 📖 [User Guide](https://nextgenusfs.github.io/annorefine/guide/installation/)
- 🐍 [Python API Reference](https://nextgenusfs.github.io/annorefine/api/functions/)
- 💻 [Command Line Reference](https://nextgenusfs.github.io/annorefine/api/overview/)
- 🚀 [funannotate2 Integration](https://nextgenusfs.github.io/annorefine/guide/python/)



## Performance

- **Multi-threaded** - Parallel processing with Rust backend
- **Memory efficient** - Streaming BAM processing
- **Scalable** - Handles mammalian-sized genomes efficiently

**Typical performance:**
- Human genome (~20K genes): 10-30 minutes on 8 cores
- Memory usage: 2-8 GB depending on genome size

## Support

- 📖 [Documentation](https://nextgenusfs.github.io/annorefine/)
- 🐛 [Bug Reports](https://github.com/nextgenusfs/annorefine/issues)
- 💬 [Discussions](https://github.com/nextgenusfs/annorefine/discussions)

## Citation

```
Palmer, J. (2025). AnnoRefine: High-performance genome annotation refinement using RNA-seq data.
GitHub: https://github.com/nextgenusfs/annorefine
```

## License

MIT License - see [LICENSE](LICENSE) file for details.

---

**Built with ❤️ in Rust** | [Documentation](https://nextgenusfs.github.io/annorefine/) | [PyPI](https://pypi.org/project/annorefine/)

