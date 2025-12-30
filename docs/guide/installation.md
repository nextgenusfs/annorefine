# Installation

## Requirements

- Python 3.8 or later
- pip package manager

## Install from PyPI

The easiest way to install AnnoRefine is from PyPI:

```bash
pip install annorefine
```

This will install both the command-line tools and Python bindings.

## Install from Source

To install from source (for development or latest features):

```bash
# Clone the repository
git clone https://github.com/nextgenusfs/annorefine.git
cd annorefine

# Install Rust (if not already installed)
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

# Install maturin for building Python extensions
pip install maturin

# Build and install in development mode
maturin develop --release

# Or build a wheel
maturin build --release
pip install target/wheels/annorefine-*.whl
```

## Verify Installation

Check that AnnoRefine is installed correctly:

```bash
# Check command-line tool
annorefine --version

# Check Python bindings
python -c "import annorefine; print(annorefine.version())"
```

## Dependencies

AnnoRefine has minimal dependencies:

- **Rust dependencies** (automatically handled during build):
  - rust-htslib for BAM file processing
  - rayon for parallel processing
  - clap for command-line interface

- **Python dependencies** (automatically installed):
  - No additional Python dependencies required!

## System Requirements

- **Memory**: Depends on genome size and BAM file size. Typically 2-8 GB for most genomes.
- **Disk Space**: Minimal, but ensure sufficient space for output files.
- **CPU**: Multi-core recommended for parallel processing (use `--threads` option).

## Platform Support

AnnoRefine is tested on:

- Linux (Ubuntu 20.04+, CentOS 7+)
- macOS (10.15+)
- Windows (via WSL2)

## Troubleshooting

### Build Errors

If you encounter build errors when installing from source:

1. Ensure Rust is up to date:
   ```bash
   rustup update
   ```

2. Ensure you have required system libraries:
   ```bash
   # Ubuntu/Debian
   sudo apt-get install build-essential libssl-dev pkg-config
   
   # macOS
   brew install openssl pkg-config
   ```

### Import Errors

If you get import errors in Python:

1. Ensure the package is installed:
   ```bash
   pip list | grep annorefine
   ```

2. Try reinstalling:
   ```bash
   pip uninstall annorefine
   pip install annorefine
   ```

## Next Steps

- [BAM to Hints Guide](bam2hints.md) - Learn how to convert BAM files to hints
- [UTR Refinement Guide](utrs.md) - Learn how to refine UTRs
- [Python API Guide](python.md) - Learn how to use the Python API

