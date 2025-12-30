//! Core data structures for genome annotation refinement

use std::collections::HashMap;
use thiserror::Error;

/// Errors that can occur during annotation refinement
#[derive(Error, Debug)]
pub enum AnnoRefineError {
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),

    #[error("FASTA parsing error: {0}")]
    FastaParse(String),

    #[error("GFF3 parsing error: {0}")]
    Gff3Parse(String),

    #[error("BAM parsing error: {0}")]
    BamParse(String),

    #[error("Invalid input: {0}")]
    InvalidInput(String),

    #[error("Invalid gene model: {0}")]
    InvalidGeneModel(String),

    #[error("Refinement error: {0}")]
    Refinement(String),
}

pub type Result<T> = std::result::Result<T, AnnoRefineError>;

/// Represents a genomic coordinate
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct GenomicPosition {
    pub chromosome: u32, // Index into chromosome names
    pub position: u64,
}

/// Represents a genomic interval
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GenomicInterval {
    pub chromosome: String,
    pub start: u64, // 1-based, inclusive
    pub end: u64,   // 1-based, inclusive
    pub strand: Strand,
}

/// DNA strand orientation
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Strand {
    Forward,
    Reverse,
    Unknown,
}

impl std::fmt::Display for Strand {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Strand::Forward => write!(f, "+"),
            Strand::Reverse => write!(f, "-"),
            Strand::Unknown => write!(f, "."),
        }
    }
}

impl std::str::FromStr for Strand {
    type Err = AnnoRefineError;

    fn from_str(s: &str) -> Result<Self> {
        match s {
            "+" => Ok(Strand::Forward),
            "-" => Ok(Strand::Reverse),
            "." => Ok(Strand::Unknown),
            _ => Err(AnnoRefineError::Gff3Parse(format!("Invalid strand: {}", s))),
        }
    }
}

/// Represents a genome sequence
#[derive(Debug, Clone)]
pub struct GenomeSequence {
    pub id: String,
    pub description: Option<String>,
    pub sequence: Vec<u8>,
}

/// Collection of genome sequences
#[derive(Debug)]
pub struct Genome {
    pub sequences: HashMap<String, GenomeSequence>,
    pub sequence_order: Vec<String>,
}

impl Genome {
    pub fn new() -> Self {
        Self {
            sequences: HashMap::new(),
            sequence_order: Vec::new(),
        }
    }

    pub fn add_sequence(&mut self, sequence: GenomeSequence) {
        self.sequence_order.push(sequence.id.clone());
        self.sequences.insert(sequence.id.clone(), sequence);
    }

    pub fn get_sequence(&self, id: &str) -> Option<&GenomeSequence> {
        self.sequences.get(id)
    }

    pub fn get_subsequence(&self, interval: &GenomicInterval) -> Result<Vec<u8>> {
        let seq = self.get_sequence(&interval.chromosome).ok_or_else(|| {
            AnnoRefineError::FastaParse(format!("Chromosome not found: {}", interval.chromosome))
        })?;

        let start = (interval.start - 1) as usize; // Convert to 0-based
        let end = interval.end as usize;

        if start >= seq.sequence.len() || end > seq.sequence.len() {
            return Err(AnnoRefineError::FastaParse(format!(
                "Interval out of bounds: {}:{}-{}",
                interval.chromosome, interval.start, interval.end
            )));
        }

        let mut subseq = seq.sequence[start..end].to_vec();

        // Reverse complement if on reverse strand
        if interval.strand == Strand::Reverse {
            subseq = reverse_complement(&subseq);
        }

        Ok(subseq)
    }
}

/// Reverse complement a DNA sequence
fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&base| match base.to_ascii_uppercase() {
            b'A' => b'T',
            b'T' => b'A',
            b'G' => b'C',
            b'C' => b'G',
            b'N' => b'N',
            _ => base,
        })
        .collect()
}

/// GFF3 feature types
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FeatureType {
    Gene,
    Mrna,
    Exon,
    Cds,
    FivePrimeUtr,
    ThreePrimeUtr,
    Other(String),
}

impl std::fmt::Display for FeatureType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            FeatureType::Gene => write!(f, "gene"),
            FeatureType::Mrna => write!(f, "mRNA"),
            FeatureType::Exon => write!(f, "exon"),
            FeatureType::Cds => write!(f, "CDS"),
            FeatureType::FivePrimeUtr => write!(f, "five_prime_UTR"),
            FeatureType::ThreePrimeUtr => write!(f, "three_prime_UTR"),
            FeatureType::Other(s) => write!(f, "{}", s),
        }
    }
}

impl std::str::FromStr for FeatureType {
    type Err = AnnoRefineError;

    fn from_str(s: &str) -> Result<Self> {
        match s {
            "gene" => Ok(FeatureType::Gene),
            "mRNA" => Ok(FeatureType::Mrna),
            "exon" => Ok(FeatureType::Exon),
            "CDS" => Ok(FeatureType::Cds),
            "five_prime_UTR" => Ok(FeatureType::FivePrimeUtr),
            "three_prime_UTR" => Ok(FeatureType::ThreePrimeUtr),
            _ => Ok(FeatureType::Other(s.to_string())),
        }
    }
}

/// Represents a GFF3 feature
#[derive(Debug, Clone)]
pub struct Gff3Feature {
    pub seqid: String,
    pub source: String,
    pub feature_type: FeatureType,
    pub start: u64,
    pub end: u64,
    pub score: Option<f64>,
    pub strand: Strand,
    pub phase: Option<u8>,
    pub attributes: HashMap<String, Vec<String>>,
}

impl Gff3Feature {
    pub fn get_id(&self) -> Option<&str> {
        self.attributes.get("ID")?.first().map(|s| s.as_str())
    }

    pub fn get_parent(&self) -> Option<&str> {
        self.attributes.get("Parent")?.first().map(|s| s.as_str())
    }

    pub fn get_name(&self) -> Option<&str> {
        self.attributes.get("Name")?.first().map(|s| s.as_str())
    }

    pub fn interval(&self) -> GenomicInterval {
        GenomicInterval {
            chromosome: self.seqid.clone(),
            start: self.start,
            end: self.end,
            strand: self.strand,
        }
    }
}

/// Represents an RNA-seq alignment
#[derive(Debug, Clone)]
pub struct RnaSeqAlignment {
    pub read_name: String,
    pub chromosome: String,
    pub start: u64,
    pub end: u64,
    pub strand: Strand,
    pub cigar: String,
    pub mapping_quality: u8,
    pub splice_junctions: Vec<SpliceJunction>,
    pub is_paired: bool,
    pub is_first_in_pair: bool,
}

/// Represents a splice junction from RNA-seq data
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SpliceJunction {
    pub chromosome: String,
    pub donor_pos: u64,    // Last position of upstream exon
    pub acceptor_pos: u64, // First position of downstream exon
    pub strand: Strand,
    pub support_count: u32,
}

/// Represents a gene model with hierarchical structure
#[derive(Debug, Clone)]
pub struct GeneModel {
    pub id: String,
    pub name: Option<String>,
    pub chromosome: String,
    pub start: u64,
    pub end: u64,
    pub strand: Strand,
    pub transcripts: Vec<Transcript>,
    pub original_attributes: HashMap<String, Vec<String>>,
    pub has_structural_changes: bool,
    pub original_cds_lengths: HashMap<String, u64>, // transcript_id -> original CDS length
}

impl GeneModel {
    pub fn interval(&self) -> GenomicInterval {
        GenomicInterval {
            chromosome: self.chromosome.clone(),
            start: self.start,
            end: self.end,
            strand: self.strand,
        }
    }

    pub fn update_boundaries(&mut self) {
        if self.transcripts.is_empty() {
            return;
        }

        self.start = self.transcripts.iter().map(|t| t.start).min().unwrap();
        self.end = self.transcripts.iter().map(|t| t.end).max().unwrap();
    }
}

/// Represents a transcript (mRNA) within a gene
#[derive(Debug, Clone)]
pub struct Transcript {
    pub id: String,
    pub name: Option<String>,
    pub start: u64,
    pub end: u64,
    pub exons: Vec<Exon>,
    pub cds_regions: Vec<CdsRegion>,
    pub five_prime_utr: Option<Vec<GenomicInterval>>,
    pub three_prime_utr: Option<Vec<GenomicInterval>>,
    pub original_attributes: HashMap<String, Vec<String>>,
}

impl Transcript {
    pub fn update_boundaries(&mut self) {
        if self.exons.is_empty() {
            return;
        }

        self.start = self.exons.iter().map(|e| e.start).min().unwrap();
        self.end = self.exons.iter().map(|e| e.end).max().unwrap();
    }

    pub fn get_splice_junctions(&self, strand: Strand) -> Vec<SpliceJunction> {
        let mut junctions = Vec::new();

        if self.exons.len() < 2 {
            return junctions;
        }

        for i in 0..self.exons.len() - 1 {
            let donor_pos = self.exons[i].end;
            let acceptor_pos = self.exons[i + 1].start;

            junctions.push(SpliceJunction {
                chromosome: "".to_string(), // Will be filled by caller
                donor_pos,
                acceptor_pos,
                strand,
                support_count: 0, // Will be updated based on RNA-seq evidence
            });
        }

        junctions
    }
}

/// Represents an exon
#[derive(Debug, Clone)]
pub struct Exon {
    pub id: Option<String>,
    pub start: u64,
    pub end: u64,
}

/// Represents a CDS region
#[derive(Debug, Clone)]
pub struct CdsRegion {
    pub id: Option<String>,
    pub start: u64,
    pub end: u64,
    pub phase: Option<u8>,
}

/// Strand bias detection result
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum StrandBias {
    /// Unstranded RNA-seq data (roughly equal forward/reverse reads)
    Unstranded,
    /// Forward-stranded RNA-seq data (bias towards forward reads)
    ForwardStranded,
    /// Reverse-stranded RNA-seq data (bias towards reverse reads)
    ReverseStranded,
}

impl std::fmt::Display for StrandBias {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            StrandBias::Unstranded => write!(f, "unstranded"),
            StrandBias::ForwardStranded => write!(f, "forward-stranded"),
            StrandBias::ReverseStranded => write!(f, "reverse-stranded"),
        }
    }
}

/// Comprehensive RNA-seq library specification
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LibraryType {
    /// Auto-detect library type
    Auto,
    /// Single-end forward stranded (F)
    SingleForward,
    /// Single-end reverse stranded (R)
    SingleReverse,
    /// Single-end unstranded (U)
    SingleUnstranded,
    /// Paired-end first-in-pair reverse (RF)
    PairedRF,
    /// Paired-end first-in-pair forward (FR)
    PairedFR,
    /// Paired-end unstranded (UU)
    PairedUnstranded,
}

impl std::fmt::Display for LibraryType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            LibraryType::Auto => write!(f, "auto"),
            LibraryType::SingleForward => write!(f, "F (single-end forward)"),
            LibraryType::SingleReverse => write!(f, "R (single-end reverse)"),
            LibraryType::SingleUnstranded => write!(f, "U (single-end unstranded)"),
            LibraryType::PairedRF => write!(f, "RF (paired-end reverse/forward)"),
            LibraryType::PairedFR => write!(f, "FR (paired-end forward/reverse)"),
            LibraryType::PairedUnstranded => write!(f, "UU (paired-end unstranded)"),
        }
    }
}

impl std::str::FromStr for LibraryType {
    type Err = AnnoRefineError;

    fn from_str(s: &str) -> Result<Self> {
        match s.to_uppercase().as_str() {
            "AUTO" => Ok(LibraryType::Auto),
            "F" => Ok(LibraryType::SingleForward),
            "R" => Ok(LibraryType::SingleReverse),
            "U" => Ok(LibraryType::SingleUnstranded),
            "RF" => Ok(LibraryType::PairedRF),
            "FR" => Ok(LibraryType::PairedFR),
            "UU" => Ok(LibraryType::PairedUnstranded),
            _ => Err(AnnoRefineError::InvalidInput(format!(
                "Invalid library type: {}. Valid options: auto, F, R, U, RF, FR, UU", s
            ))),
        }
    }
}

/// Configuration parameters for refinement
#[derive(Debug, Clone)]
pub struct RefinementConfig {
    pub min_coverage: u32,
    pub min_splice_support: u32,
    pub max_utr_extension: u32,
    pub enable_novel_gene_detection: bool,
    pub min_novel_gene_coverage: u32,
    pub min_novel_gene_length: u64,
    pub min_exon_length: u64,
    pub validate_splice_sites: bool,
    pub strand_bias_threshold: f64,
    pub max_reads_for_strand_detection: u32,
    pub library_type: LibraryType,
}

impl Default for RefinementConfig {
    fn default() -> Self {
        Self {
            min_coverage: 5,
            min_splice_support: 3,
            max_utr_extension: 1000,
            enable_novel_gene_detection: false,
            min_novel_gene_coverage: 10,
            min_novel_gene_length: 300,
            min_exon_length: 50,
            validate_splice_sites: true,
            strand_bias_threshold: 0.65, // 65% threshold for strand bias detection
            max_reads_for_strand_detection: 10000, // Sample up to 10k reads for detection
            library_type: LibraryType::Auto,
        }
    }
}

/// Represents a potential novel gene detected from RNA-seq evidence
#[derive(Debug, Clone)]
pub struct NovelGeneCandidate {
    pub chromosome: String,
    pub start: u64,
    pub end: u64,
    pub strand: Strand,
    pub exons: Vec<Exon>,
    pub splice_junctions: Vec<SpliceJunction>,
    pub coverage: f64,
    pub confidence_score: f64,
}

/// Represents a cluster of alignments that might form a novel gene
#[derive(Debug, Clone)]
pub struct AlignmentCluster {
    pub chromosome: String,
    pub start: u64,
    pub end: u64,
    pub strand: Strand,
    pub splice_junctions: Vec<SpliceJunction>,
    pub coverage_profile: Vec<u32>,
    pub alignment_count: u32,
}

/// Types of hints that can be generated for Augustus
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum HintType {
    /// Intron hint - indicates a likely intron region
    Intron,
    /// Exon hint - indicates a likely exon region (internal exons)
    Exon,
    /// Exonpart hint - indicates part of an exon (for terminal exons)
    Exonpart,
    /// Donor splice site hint
    DSS,
    /// Acceptor splice site hint
    ASS,
}

impl std::fmt::Display for HintType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            HintType::Intron => write!(f, "intron"),
            HintType::Exon => write!(f, "exon"),
            HintType::Exonpart => write!(f, "exonpart"),
            HintType::DSS => write!(f, "dss"),
            HintType::ASS => write!(f, "ass"),
        }
    }
}

/// Represents a hint for Augustus gene prediction
#[derive(Debug, Clone)]
pub struct AugustusHint {
    pub chromosome: String,
    pub start: u64,        // 1-based, inclusive
    pub end: u64,          // 1-based, inclusive
    pub hint_type: HintType,
    pub strand: Strand,
    pub score: f64,
    pub source: String,    // Source identifier (e.g., "E" for EST, "b2h" for bam2hints)
    pub priority: u32,     // Priority level for the hint
    pub multiplicity: u32, // Number of supporting alignments (for intron hints)
}

impl AugustusHint {
    /// Create a new hint
    pub fn new(
        chromosome: String,
        start: u64,
        end: u64,
        hint_type: HintType,
        strand: Strand,
        source: String,
        priority: u32,
    ) -> Self {
        Self {
            chromosome,
            start,
            end,
            hint_type,
            strand,
            score: 0.0,
            source,
            priority,
            multiplicity: 1,
        }
    }

    /// Create an intron hint with multiplicity
    pub fn intron_hint(
        chromosome: String,
        start: u64,
        end: u64,
        strand: Strand,
        source: String,
        priority: u32,
        multiplicity: u32,
    ) -> Self {
        Self {
            chromosome,
            start,
            end,
            hint_type: HintType::Intron,
            strand,
            score: 0.0,
            source,
            priority,
            multiplicity,
        }
    }

    /// Format as GFF3 line for Augustus
    pub fn to_gff3_line(&self) -> String {
        let mut attributes = Vec::new();

        if self.multiplicity > 1 {
            attributes.push(format!("mult={}", self.multiplicity));
        }
        attributes.push(format!("pri={}", self.priority));
        attributes.push(format!("src={}", self.source));

        let attributes_str = attributes.join(";");

        format!(
            "{}\tb2h\t{}\t{}\t{}\t{}\t{}\t.\t{}",
            self.chromosome,
            self.hint_type,
            self.start,
            self.end,
            self.score,
            self.strand,
            attributes_str
        )
    }
}

/// Configuration for bam2hints conversion
#[derive(Debug, Clone)]
pub struct Bam2HintsConfig {
    /// Priority of hint group (default: 4)
    pub priority: u32,
    /// Gaps at most this length are simply closed (default: 14)
    pub max_gap_len: u32,
    /// Alignments with gaps shorter than this and longer than max_gap_len are discarded (default: 32)
    pub min_intron_len: u32,
    /// Alignments with longer gaps are discarded (default: 350000)
    pub max_intron_len: u32,
    /// Minimum length of a 'dangling' exon (default: 8)
    pub min_end_block_len: u32,
    /// Maximum length of gap in query (cDNA) sequence (default: 5)
    pub max_query_gap_len: u32,
    /// This many bp are cut off of each exonpart hint at end of alignment (default: 10)
    pub exonpart_cutoff: u32,
    /// Source identifier (default: "E")
    pub source: String,
    /// Only retrieve intron hints (default: false)
    pub introns_only: bool,
    /// Do not summarize multiple identical intron hints (default: false)
    pub no_multiplicity: bool,
    /// Only keep the strongest hint for a region (default: false)
    pub remove_redundant: bool,
    /// Maximal number of hints at a given position (0: filtering deactivated, default: 0)
    pub max_coverage: u32,
    /// Include splice site (dss, ass) hints in output (default: false)
    pub splice_sites_on: bool,
    /// Include splice sites hints from the ends of a truncated alignment (default: false)
    pub truncated_splice_sites: bool,
    /// Fill this number in the score column (default: 0.0)
    pub score: f64,
    /// Alignments that span more than this are ignored (default: 400000)
    pub max_gene_len: u32,
    /// Library type for strand-specific processing
    pub library_type: LibraryType,
    /// Detected strand bias from BAM file
    pub strand_bias: StrandBias,
}

impl Default for Bam2HintsConfig {
    fn default() -> Self {
        Self {
            priority: 4,
            max_gap_len: 14,
            min_intron_len: 32,
            max_intron_len: 350000,
            min_end_block_len: 8,
            max_query_gap_len: 5,
            exonpart_cutoff: 10,
            source: "E".to_string(),
            introns_only: true,
            no_multiplicity: false,
            remove_redundant: false,
            max_coverage: 0,
            splice_sites_on: false,
            truncated_splice_sites: false,
            score: 0.0,
            max_gene_len: 400000,
            library_type: LibraryType::Auto,
            strand_bias: StrandBias::Unstranded,
        }
    }
}

/// Represents alignment blocks extracted from CIGAR parsing (similar to PSL format)
#[derive(Debug, Clone)]
pub struct AlignmentBlock {
    /// Length of the block
    pub length: u32,
    /// 1-based start position in query sequence
    pub query_start: u64,
    /// 1-based start position in target sequence
    pub target_start: u64,
}

/// Represents filtered alignment blocks for hint generation
#[derive(Debug, Clone)]
pub struct FilteredBlock {
    /// 1-based start coordinate of the block
    pub start: u64,
    /// 1-based end coordinate of the block
    pub end: u64,
    /// Whether the gap following this block is considered an intron
    pub following_intron_ok: bool,
}
