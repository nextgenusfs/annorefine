# Python Functions Reference

Complete reference for all Python functions in AnnoRefine.

## Hint Generation

### bam2hints()

Convert BAM alignments to Augustus/GeneMark hints format.

**Signature:**
```python
annorefine.bam2hints(
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
    source: str = 'E',
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
    contig_map: dict = None
) -> dict
```

**Parameters:**

- `bam_file` (str): Path to input BAM file (must be sorted and indexed if using contig/region filtering)
- `output_file` (str): Path for output GFF hints file
- `library_type` (str): Library strandedness specification. Options: "FR", "RF", "UU"
- `priority` (int): Priority of hint group (default: 4)
- `source` (str): Source identifier (default: "E")
- `introns_only` (bool): Only retrieve intron hints (default: False)
- `threads` (int): Number of threads to use (default: None, uses all available)
- `contig` (str): Filter to only process alignments on this contig (default: None)
- `region` (tuple): Filter to region as (contig, start, end) tuple (default: None)
- `contig_map` (dict): Dictionary to rename contigs in output. Keys are input contig names, values are output contig names (default: None, no mapping)

**Returns:** Dictionary with conversion statistics

**Example:**

```python
import annorefine

# Basic usage
result = annorefine.bam2hints(
    bam_file="alignments.bam",
    output_file="hints.gff",
    library_type="RF"
)

# Advanced usage
result = annorefine.bam2hints(
    bam_file="alignments.bam",
    output_file="hints.gff",
    library_type="RF",
    priority=4,
    source="E",
    threads=8,
    min_intron_len=32,
    max_intron_len=350000,
    contig="chr1"
)

# With contig name mapping
contig_map = {
    'NC_000001.11': 'chr1',
    'NC_000002.12': 'chr2',
    'NC_000003.12': 'chr3'
}

result = annorefine.bam2hints(
    bam_file="alignments.bam",
    output_file="hints.gff",
    library_type="RF",
    contig_map=contig_map  # Rename contigs in output
)

print(f"Generated {result['total_hints_generated']} hints")
```

## Hint Processing

### join_hints()

Join and merge hints from multiple GFF files.

**Signature:**
```python
annorefine.join_hints(
    input_files: List[str],
    output_file: str,
    introns_only: bool = False
) -> dict
```

**Parameters:**

- `input_files` (List[str]): List of input GFF hint files to join
- `output_file` (str): Path for output joined hints file
- `introns_only` (bool): If True, only output intron hints (useful for GeneMark) (default: False)

**Returns:** Dictionary with joining statistics

**Example:**

```python
import annorefine

# Join hints from multiple sources
result = annorefine.join_hints(
    input_files=[
        "bam_hints.gff",
        "protein_hints.gff",
        "transcript_hints.gff"
    ],
    output_file="joined_hints.gff"
)

# Join and filter to introns only
result = annorefine.join_hints(
    input_files=["bam_hints.gff", "protein_hints.gff"],
    output_file="genemark_hints.gff",
    introns_only=True
)

print(f"Merged {result['total_input_hints']} into {result['output_hints']} hints")
```

### filter_hints()

Filter hints by type, multiplicity, or contig.

**Signature:**
```python
annorefine.filter_hints(
    input_file: str,
    output_file: str,
    hint_types: List[str] = None,
    min_mult: int = None,
    contig: str = None
) -> dict
```

**Parameters:**

- `input_file` (str): Input GFF file
- `output_file` (str): Output GFF file
- `hint_types` (List[str]): Filter by hint types (e.g., ["intron", "exon"]) (default: None)
- `min_mult` (int): Minimum multiplicity (default: None)
- `contig` (str): Filter to specific contig (default: None)

**Returns:** Dictionary with filter statistics

**Example:**

```python
import annorefine

# Filter by hint type
result = annorefine.filter_hints(
    input_file="all_hints.gff",
    output_file="intron_hints.gff",
    hint_types=["intron"]
)

# Filter by multiplicity
result = annorefine.filter_hints(
    input_file="all_hints.gff",
    output_file="high_confidence_hints.gff",
    min_mult=10
)

# Filter by contig
result = annorefine.filter_hints(
    input_file="all_hints.gff",
    output_file="chr1_hints.gff",
    contig="chr1"
)
```

## Annotation Refinement

### refine()

Convenience function for annotation refinement with keyword arguments.

**Signature:**
```python
annorefine.refine(
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
    threads: int = None
) -> dict
```

**Parameters:**

- `fasta_file` (str): Path to genome FASTA file
- `gff3_file` (str): Path to input GFF3 annotation file
- `bam_file` (str): Path to RNA-seq BAM file
- `output_file` (str): Path for output refined GFF3 file
- `min_coverage` (int): Minimum coverage for UTR extension (default: 5)
- `min_splice_support` (int): Minimum reads supporting splice junction (default: 3)
- `max_utr_extension` (int): Maximum UTR extension length (default: 1000)
- `enable_novel_gene_detection` (bool): Enable novel gene detection (default: False)
- `min_novel_gene_coverage` (int): Minimum coverage for novel genes (default: 10)
- `min_novel_gene_length` (int): Minimum length for novel genes (default: 300)
- `validate_splice_sites` (bool): Validate splice sites (default: True)
- `threads` (int): Number of threads (default: None, uses all available)

**Returns:** Dictionary with refinement statistics

**Example:**

```python
import annorefine

# Basic refinement
result = annorefine.refine(
    fasta_file="genome.fa",
    gff3_file="annotations.gff3",
    bam_file="alignments.bam",
    output_file="refined.gff3"
)

# With novel gene detection
result = annorefine.refine(
    fasta_file="genome.fa",
    gff3_file="annotations.gff3",
    bam_file="alignments.bam",
    output_file="refined.gff3",
    enable_novel_gene_detection=True,
    min_novel_gene_coverage=10,
    min_novel_gene_length=300,
    threads=8
)

print(f"Processed {result['genes_processed']} genes")
print(f"Novel genes: {result['novel_genes_detected']}")
```

### refine_annotations()

Refine annotations using a configuration object.

**Signature:**
```python
annorefine.refine_annotations(
    fasta_file: str,
    gff3_file: str,
    bam_file: str,
    output_file: str,
    config: RefinementConfig = None,
    threads: int = None
) -> dict
```

**Parameters:**

- `fasta_file` (str): Path to genome FASTA file
- `gff3_file` (str): Path to input GFF3 annotation file
- `bam_file` (str): Path to RNA-seq BAM file
- `output_file` (str): Path for output refined GFF3 file
- `config` (RefinementConfig): Configuration object (default: None, uses defaults)
- `threads` (int): Number of threads (default: None, uses all available)

**Returns:** Dictionary with refinement statistics

**Example:**

```python
import annorefine

# Using configuration object
config = annorefine.RefinementConfig(
    min_coverage=10,
    min_splice_support=5,
    max_utr_extension=2000,
    enable_novel_gene_detection=True
)

result = annorefine.refine_annotations(
    fasta_file="genome.fa",
    gff3_file="annotations.gff3",
    bam_file="alignments.bam",
    output_file="refined.gff3",
    config=config,
    threads=8
)
```

## Configuration Classes

### RefinementConfig

Configuration object for annotation refinement.

**Constructor:**
```python
annorefine.RefinementConfig(
    min_coverage: int = 5,
    min_splice_support: int = 3,
    max_utr_extension: int = 1000,
    enable_novel_gene_detection: bool = False,
    min_novel_gene_coverage: int = 10,
    min_novel_gene_length: int = 300,
    min_exon_length: int = 50,
    validate_splice_sites: bool = True,
    strand_bias_threshold: float = 0.65,
    max_reads_for_strand_detection: int = 10000
)
```

**Example:**

```python
import annorefine

config = annorefine.RefinementConfig(
    min_coverage=10,
    min_splice_support=5,
    max_utr_extension=2000,
    enable_novel_gene_detection=True,
    min_novel_gene_coverage=15,
    min_novel_gene_length=500,
    min_exon_length=50,
    validate_splice_sites=True,
    strand_bias_threshold=0.7
)
```

### Bam2HintsConfig

Configuration object for BAM to hints conversion.

**Constructor:**
```python
annorefine.Bam2HintsConfig(
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
    contig_map: dict = None
)
```

**Example:**

```python
import annorefine

config = annorefine.Bam2HintsConfig(
    priority=4,
    max_gap_len=14,
    min_intron_len=32,
    max_intron_len=350000,
    min_end_block_len=8,
    source="E",
    introns_only=False
)

# With contig mapping
contig_map = {'NC_000001.11': 'chr1', 'NC_000002.12': 'chr2'}
config = annorefine.Bam2HintsConfig(
    priority=4,
    source="E",
    contig_map=contig_map
)
```

## Utility Functions

### version()

Get the AnnoRefine version string.

```python
import annorefine

version = annorefine.version()
print(f"AnnoRefine version: {version}")
```

**Returns:** `str` - Version string (e.g., "2025.9.18")

### current_num_threads()

Get the current number of threads configured for parallel processing.

```python
import annorefine

num_threads = annorefine.current_num_threads()
print(f"Using {num_threads} threads")
```

**Returns:** `int` - Number of threads

## Return Values

All main functions return dictionaries with statistics:

### bam2hints() Returns

```python
{
    'total_hints_generated': int,  # Total hints generated
    'intron_hints': int,           # Number of intron hints
    'exon_hints': int,             # Number of exon hints
    'exonpart_hints': int,         # Number of exonpart hints
    'dss_hints': int,              # Number of donor splice site hints
    'ass_hints': int               # Number of acceptor splice site hints
}
```

### join_hints() Returns

```python
{
    'input_files': int,            # Number of input files
    'total_input_hints': int,      # Total hints from all inputs
    'output_hints': int,           # Hints after merging
    'output_file': str             # Output file path
}
```

### refine() Returns

```python
{
    'genes_processed': int,        # Number of genes processed
    'utrs_extended': int,          # Number of UTRs extended
    'novel_genes_detected': int,   # Number of novel genes found
    'splice_junctions_refined': int  # Number of splice junctions refined
}
```

