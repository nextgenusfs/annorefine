//! BAM to hints conversion functionality for Augustus gene prediction
//! 
//! This module implements the core algorithms to convert BAM alignments into
//! hints for Augustus gene prediction, similar to the original bam2hints tool.

use crate::types::{
    AlignmentBlock, AugustusHint, Bam2HintsConfig, FilteredBlock, HintType, Result, 
    RnaSeqAlignment, Strand, AnnoRefineError
};
use log::{debug, warn};

use std::collections::HashMap;

/// Convert BAM alignments to Augustus hints
pub struct Bam2HintsConverter {
    config: Bam2HintsConfig,
    hint_lists: HintLists,
}

/// Container for different types of hints
#[derive(Debug, Default)]
struct HintLists {
    intron_hints: Vec<AugustusHint>,
    exon_hints: Vec<AugustusHint>,
    exonpart_hints: Vec<AugustusHint>,
    dss_hints: Vec<AugustusHint>,
    ass_hints: Vec<AugustusHint>,
}

impl Bam2HintsConverter {
    /// Create a new converter with the given configuration
    pub fn new(config: Bam2HintsConfig) -> Self {
        Self {
            config,
            hint_lists: HintLists::default(),
        }
    }

    /// Process a single RNA-seq alignment and generate hints
    pub fn process_alignment(&mut self, alignment: &RnaSeqAlignment) -> Result<()> {
        // Filter by mapping quality
        if alignment.mapping_quality < self.config.min_mapping_quality {
            debug!("Dropping alignment {} - mapping quality {} below threshold {}",
                   alignment.read_name, alignment.mapping_quality, self.config.min_mapping_quality);
            return Ok(());
        }

        // Filter multi-mapping reads if enabled
        if self.config.filter_multimappers {
            if let Some(num_hits) = alignment.num_hits {
                if num_hits > 1 {
                    debug!("Dropping alignment {} - multi-mapper with {} hits",
                           alignment.read_name, num_hits);
                    return Ok(());
                }
            }
        }

        // Filter alignment based on library type and strand
        if !self.should_process_alignment(alignment) {
            return Ok(()); // Skip this alignment
        }

        // Parse CIGAR to extract alignment blocks (similar to PSL format)
        let blocks = self.parse_cigar_to_blocks(alignment)?;
        
        if blocks.is_empty() {
            return Ok(());
        }

        // Check if alignment exceeds maximum gene length
        let alignment_span = blocks.last().unwrap().target_start + blocks.last().unwrap().length as u64 
                            - blocks.first().unwrap().target_start;
        
        if alignment_span > self.config.max_gene_len as u64 {
            debug!("Dropping alignment {} - exceeds max gene length", alignment.read_name);
            return Ok(());
        }

        // Filter blocks based on gap lengths and quality criteria
        let filtered_blocks = self.filter_blocks(&blocks)?;
        
        if filtered_blocks.is_empty() {
            return Ok(());
        }

        // Determine the correct strand for hints based on library protocol
        let hint_strand = self.determine_hint_strand(alignment);

        // Generate hints from filtered blocks
        self.generate_hints_from_blocks(&filtered_blocks, &alignment.chromosome, hint_strand)?;

        Ok(())
    }

    /// Determine if an alignment should be processed based on library type and strand
    fn should_process_alignment(&self, alignment: &RnaSeqAlignment) -> bool {
        use crate::types::StrandBias;

        match self.config.strand_bias {
            StrandBias::Unstranded => {
                // For unstranded libraries, process all reads
                true
            }
            StrandBias::ForwardStranded => {
                // For FR libraries, process all reads (both mates contribute)
                true
            }
            StrandBias::ReverseStranded => {
                // For RF libraries, only process first-in-pair reads for strand specificity
                // This is the standard approach for RF stranded RNA-seq
                if alignment.is_paired {
                    alignment.is_first_in_pair
                } else {
                    // Process single-end reads
                    true
                }
            }
        }
    }

    /// Determine the correct strand for hints based on library protocol and alignment strand
    fn determine_hint_strand(&self, alignment: &RnaSeqAlignment) -> Strand {
        use crate::types::StrandBias;
        use log::debug;

        let hint_strand = match self.config.strand_bias {
            StrandBias::Unstranded => {
                // For unstranded libraries, hints should be strand-neutral
                Strand::Unknown
            }
            StrandBias::ForwardStranded => {
                // For FR libraries (fr-secondstrand):
                // - First-in-pair reads align same direction as gene
                // - Second-in-pair reads align opposite direction to gene
                if alignment.is_paired {
                    if alignment.is_first_in_pair {
                        // First-in-pair: alignment strand = gene strand
                        alignment.strand
                    } else {
                        // Second-in-pair: alignment strand is opposite to gene strand
                        match alignment.strand {
                            Strand::Forward => Strand::Reverse,
                            Strand::Reverse => Strand::Forward,
                            Strand::Unknown => Strand::Unknown,
                        }
                    }
                } else {
                    // Single-end: alignment strand = gene strand
                    alignment.strand
                }
            }
            StrandBias::ReverseStranded => {
                // For RF libraries (fr-firststrand):
                // We only process first-in-pair reads, and their alignment strand is opposite to gene strand
                if alignment.is_paired {
                    // Should only be first-in-pair reads due to filtering in should_process_alignment
                    debug_assert!(alignment.is_first_in_pair, "RF library should only process first-in-pair reads");
                    // First-in-pair: alignment strand is opposite to gene strand
                    match alignment.strand {
                        Strand::Forward => Strand::Reverse,
                        Strand::Reverse => Strand::Forward,
                        Strand::Unknown => Strand::Unknown,
                    }
                } else {
                    // Single-end: alignment strand is opposite to gene strand
                    match alignment.strand {
                        Strand::Forward => Strand::Reverse,
                        Strand::Reverse => Strand::Forward,
                        Strand::Unknown => Strand::Unknown,
                    }
                }
            }
        };

        debug!("Strand determination: read={}, alignment={:?}, is_paired={}, is_first_in_pair={}, library={:?}, hint={:?}",
               alignment.read_name, alignment.strand, alignment.is_paired, alignment.is_first_in_pair, self.config.strand_bias, hint_strand);
        hint_strand
    }

    /// Parse CIGAR string to extract alignment blocks
    fn parse_cigar_to_blocks(&self, alignment: &RnaSeqAlignment) -> Result<Vec<AlignmentBlock>> {
        let mut blocks: Vec<AlignmentBlock> = Vec::new();
        let mut query_offset = 1u64; // 1-based position in query
        let mut target_offset = alignment.start; // 1-based position in target

        // Parse the CIGAR string manually since we have it as a string
        let cigar_ops = self.parse_cigar_string(&alignment.cigar)?;

        for (op_type, length) in cigar_ops {
            match op_type {
                'M' | 'X' | '=' => {
                    // Match/mismatch - this is an aligned block
                    if let Some(last_block) = blocks.last_mut() {
                        // Check if we can extend the previous block
                        if last_block.target_start + last_block.length as u64 == target_offset &&
                           last_block.query_start + last_block.length as u64 == query_offset {
                            // Extend existing block
                            last_block.length += length;
                        } else {
                            // Create new block
                            blocks.push(AlignmentBlock {
                                length,
                                query_start: query_offset,
                                target_start: target_offset,
                            });
                        }
                    } else {
                        // First block
                        blocks.push(AlignmentBlock {
                            length,
                            query_start: query_offset,
                            target_start: target_offset,
                        });
                    }
                    
                    query_offset += length as u64;
                    target_offset += length as u64;
                }
                'D' | 'N' => {
                    // Deletion or RefSkip (intron) - advance target position only
                    target_offset += length as u64;
                }
                'I' => {
                    // Insertion - advance query position only
                    query_offset += length as u64;
                }
                'H' | 'S' | 'P' => {
                    // Hard clip, soft clip, padding - ignore
                }
                _ => {
                    warn!("Unknown CIGAR operation: {}", op_type);
                }
            }
        }

        Ok(blocks)
    }

    /// Parse CIGAR string into operations
    fn parse_cigar_string(&self, cigar: &str) -> Result<Vec<(char, u32)>> {
        let mut operations = Vec::new();
        let mut current_number = String::new();

        for ch in cigar.chars() {
            if ch.is_ascii_digit() {
                current_number.push(ch);
            } else {
                if !current_number.is_empty() {
                    let length = current_number.parse::<u32>().map_err(|e| {
                        AnnoRefineError::BamParse(format!("Invalid CIGAR length: {}", e))
                    })?;
                    operations.push((ch, length));
                    current_number.clear();
                }
            }
        }

        Ok(operations)
    }

    /// Filter alignment blocks based on gap lengths and quality criteria
    fn filter_blocks(&self, blocks: &[AlignmentBlock]) -> Result<Vec<FilteredBlock>> {
        let mut filtered = Vec::new();

        for (i, block) in blocks.iter().enumerate() {
            let gap_len = if i == 0 {
                self.config.min_intron_len as u64 // First block, use minimum intron length
            } else {
                // Calculate gap to previous block
                let prev_block = &blocks[i - 1];
                let prev_end = prev_block.target_start + prev_block.length as u64 - 1;
                block.target_start - prev_end - 1
            };

            let block_start = block.target_start;
            let block_end = block.target_start + block.length as u64 - 1;

            // Decide based on gap length
            if gap_len >= self.config.min_intron_len as u64 && gap_len <= self.config.max_intron_len as u64 {
                // Gap represents an intron, add new block
                let following_intron_ok = if i < blocks.len() - 1 {
                    // Check query gap length for next block
                    let next_block = &blocks[i + 1];
                    let query_gap = next_block.query_start - (block.query_start + block.length as u64);
                    query_gap <= self.config.max_query_gap_len as u64
                } else {
                    false
                };

                filtered.push(FilteredBlock {
                    start: block_start,
                    end: block_end,
                    following_intron_ok,
                });
            } else if gap_len <= self.config.max_gap_len as u64 && !filtered.is_empty() {
                // Gap represents a deletion, expand previous block
                if let Some(last_block) = filtered.last_mut() {
                    last_block.end = block_end;
                    
                    // Update following_intron_ok
                    if i < blocks.len() - 1 {
                        let next_block = &blocks[i + 1];
                        let query_gap = next_block.query_start - (block.query_start + block.length as u64);
                        last_block.following_intron_ok = query_gap <= self.config.max_query_gap_len as u64;
                    } else {
                        last_block.following_intron_ok = false;
                    }
                }
            } else if gap_len > self.config.max_intron_len as u64 {
                // Gap is too long, drop alignment
                return Ok(Vec::new());
            }
        }

        Ok(filtered)
    }

    /// Generate hints from filtered blocks
    fn generate_hints_from_blocks(
        &mut self,
        blocks: &[FilteredBlock],
        chromosome: &str,
        strand: Strand,
    ) -> Result<()> {
        for (i, block) in blocks.iter().enumerate() {
            if i == 0 {
                // First block of alignment
                if blocks.len() == 1 && !self.config.introns_only {
                    // Single-block alignment - only exonpart hint
                    if block.end - block.start + 1 >= 2 * self.config.exonpart_cutoff as u64 {
                        let hint = AugustusHint::new(
                            chromosome.to_string(),
                            block.start + self.config.exonpart_cutoff as u64,
                            block.end - self.config.exonpart_cutoff as u64,
                            HintType::Exonpart,
                            strand,
                            self.config.source.clone(),
                            self.config.priority,
                        );
                        self.hint_lists.exonpart_hints.push(hint);
                    }
                } else if block.end - block.start + 1 >= self.config.min_end_block_len as u64 {
                    // First block of multi-block alignment with minimum length
                    if !self.config.introns_only &&
                       block.end - block.start + 1 >= 2 * self.config.exonpart_cutoff as u64 {
                        // Exonpart hint
                        let hint_start = block.start + self.config.exonpart_cutoff as u64;
                        if hint_start < block.end {
                            let hint = AugustusHint::new(
                                chromosome.to_string(),
                                hint_start,
                                block.end,
                                HintType::Exonpart,
                                strand,
                                self.config.source.clone(),
                                self.config.priority,
                            );
                            self.hint_lists.exonpart_hints.push(hint);
                        }
                    }

                    // Generate splice site hints if enabled
                    if self.config.splice_sites_on && !self.config.introns_only && i < blocks.len() - 1 {
                        self.add_splice_site_hints(chromosome, block.end, blocks[i + 1].start, strand);
                    }

                    // Generate intron hint if conditions are met
                    if block.following_intron_ok && i + 1 < blocks.len() &&
                       (blocks.len() > 2 && i < blocks.len() - 2 ||
                        (blocks[i + 1].end >= blocks[i + 1].start &&
                         blocks[i + 1].end - blocks[i + 1].start + 1 >= self.config.min_end_block_len as u64)) &&
                       blocks[i + 1].start > 0 {
                        let hint = AugustusHint::new(
                            chromosome.to_string(),
                            block.end + 1,
                            blocks[i + 1].start - 1,
                            HintType::Intron,
                            strand,
                            self.config.source.clone(),
                            self.config.priority,
                        );
                        self.hint_lists.intron_hints.push(hint);
                    }
                }
            } else if i == blocks.len() - 1 && !self.config.introns_only {
                // Last block of multi-block alignment
                if block.end - block.start + 1 >= self.config.min_end_block_len as u64 &&
                   block.end - block.start + 1 >= 2 * self.config.exonpart_cutoff as u64 {
                    let hint_end = block.end - self.config.exonpart_cutoff as u64;
                    if hint_end > block.start {
                        let hint = AugustusHint::new(
                            chromosome.to_string(),
                            block.start,
                            hint_end,
                            HintType::Exonpart,
                            strand,
                            self.config.source.clone(),
                            self.config.priority,
                        );
                        self.hint_lists.exonpart_hints.push(hint);
                    }
                }
            } else {
                // Inner block of alignment
                if !self.config.introns_only {
                    let hint = AugustusHint::new(
                        chromosome.to_string(),
                        block.start,
                        block.end,
                        HintType::Exon,
                        strand,
                        self.config.source.clone(),
                        self.config.priority,
                    );
                    self.hint_lists.exon_hints.push(hint);
                }

                // Generate intron hint if conditions are met
                if block.following_intron_ok && i + 1 < blocks.len() &&
                   (blocks.len() > 2 && i < blocks.len() - 2 ||
                    (blocks[i + 1].end >= blocks[i + 1].start &&
                     blocks[i + 1].end - blocks[i + 1].start + 1 >= self.config.min_end_block_len as u64)) &&
                   blocks[i + 1].start > 0 {
                    let hint = AugustusHint::new(
                        chromosome.to_string(),
                        block.end + 1,
                        blocks[i + 1].start - 1,
                        HintType::Intron,
                        strand,
                        self.config.source.clone(),
                        self.config.priority,
                    );
                    self.hint_lists.intron_hints.push(hint);

                    // Generate splice site hints if enabled
                    if self.config.splice_sites_on && !self.config.introns_only {
                        self.add_splice_site_hints(chromosome, block.end, blocks[i + 1].start, strand);
                    }
                }
            }
        }

        Ok(())
    }

    /// Add splice site hints for donor and acceptor sites
    fn add_splice_site_hints(&mut self, chromosome: &str, donor_end: u64, acceptor_start: u64, strand: Strand) {
        // Donor splice site (end of exon + 1)
        let dss_hint = AugustusHint::new(
            chromosome.to_string(),
            donor_end + 1,
            donor_end + 1,
            HintType::DSS,
            strand,
            self.config.source.clone(),
            self.config.priority,
        );
        self.hint_lists.dss_hints.push(dss_hint);

        // Acceptor splice site (start of exon - 1)
        let ass_hint = AugustusHint::new(
            chromosome.to_string(),
            acceptor_start - 1,
            acceptor_start - 1,
            HintType::ASS,
            strand,
            self.config.source.clone(),
            self.config.priority,
        );
        self.hint_lists.ass_hints.push(ass_hint);
    }

    /// Get all generated hints, applying multiplicity compression for all hint types
    pub fn get_hints(&mut self) -> Vec<AugustusHint> {
        let mut all_hints = Vec::new();

        // Process intron hints with multiplicity compression first
        let intron_hints = if !self.config.no_multiplicity {
            self.compress_intron_hints()
        } else {
            self.hint_lists.intron_hints.drain(..).collect()
        };

        // Process other hint types with deduplication and reclassification
        let mut other_hints = if !self.config.no_multiplicity {
            self.compress_other_hints()
        } else {
            let mut hints = Vec::new();
            hints.extend(self.hint_lists.exonpart_hints.drain(..));
            hints.extend(self.hint_lists.exon_hints.drain(..));
            hints.extend(self.hint_lists.dss_hints.drain(..));
            hints.extend(self.hint_lists.ass_hints.drain(..));
            hints
        };

        // Reclassify exonpart hints that are flanked by introns as exon hints
        self.reclassify_internal_exons(&mut other_hints, &intron_hints);

        // Combine all hints
        all_hints.extend(intron_hints);
        all_hints.extend(other_hints);

        all_hints
    }

    /// Reclassify exonpart hints that are flanked by introns as exon hints
    fn reclassify_internal_exons(&self, other_hints: &mut Vec<AugustusHint>, intron_hints: &[AugustusHint]) {
        // Group introns by chromosome and strand for efficient lookup
        use std::collections::HashMap;
        let mut intron_map: HashMap<(String, Strand), Vec<(u64, u64)>> = HashMap::new();

        for intron in intron_hints {
            let key = (intron.chromosome.clone(), intron.strand);
            intron_map.entry(key).or_insert_with(Vec::new).push((intron.start, intron.end));
        }

        // Sort intron coordinates for each chromosome/strand
        for intervals in intron_map.values_mut() {
            intervals.sort_by_key(|&(start, _)| start);
        }

        // Reclassify exonpart hints
        for hint in other_hints.iter_mut() {
            if hint.hint_type == HintType::Exonpart {
                let key = (hint.chromosome.clone(), hint.strand);
                if let Some(introns) = intron_map.get(&key) {
                    // Check if this exonpart is flanked by introns
                    let is_internal = self.is_flanked_by_introns(hint.start, hint.end, introns);
                    if is_internal {
                        hint.hint_type = HintType::Exon;
                    }
                }
            }
        }
    }

    /// Check if a region is flanked by introns (indicating it's an internal exon)
    fn is_flanked_by_introns(&self, exon_start: u64, exon_end: u64, introns: &[(u64, u64)]) -> bool {
        let mut has_upstream_intron = false;
        let mut has_downstream_intron = false;

        for &(intron_start, intron_end) in introns {
            // Check for upstream intron (intron ends just before exon starts)
            if intron_end + 1 == exon_start {
                has_upstream_intron = true;
            }
            // Check for downstream intron (intron starts just after exon ends)
            if intron_start == exon_end + 1 {
                has_downstream_intron = true;
            }

            // Early exit if both found
            if has_upstream_intron && has_downstream_intron {
                return true;
            }
        }

        has_upstream_intron && has_downstream_intron
    }

    /// Compress identical intron hints by counting multiplicity
    /// For unstranded libraries, merge introns with same coordinates but different strands
    fn compress_intron_hints(&mut self) -> Vec<AugustusHint> {
        use crate::types::StrandBias;

        match self.config.strand_bias {
            StrandBias::Unstranded => {
                // For unstranded libraries, merge strand-specific introns into strand-neutral hints
                self.compress_intron_hints_strand_neutral()
            }
            StrandBias::ForwardStranded | StrandBias::ReverseStranded => {
                // For stranded libraries, keep strand information and filter appropriately
                self.compress_intron_hints_strand_specific()
            }
        }
    }

    /// Compress intron hints keeping strand information (original behavior)
    fn compress_intron_hints_strand_specific(&mut self) -> Vec<AugustusHint> {
        let mut intron_counts: HashMap<(String, u64, u64, Strand), u32> = HashMap::new();

        // Count occurrences of identical intron hints
        for hint in &self.hint_lists.intron_hints {
            let key = (
                hint.chromosome.clone(),
                hint.start,
                hint.end,
                hint.strand,
            );
            *intron_counts.entry(key).or_insert(0) += 1;
        }

        // Create compressed hints with deterministic ordering
        let mut compressed = Vec::new();
        let mut sorted_keys: Vec<_> = intron_counts.keys().collect();
        sorted_keys.sort_by(|a, b| {
            // Sort by chromosome, start, end, then strand (Forward before Reverse)
            a.0.cmp(&b.0)
                .then(a.1.cmp(&b.1))
                .then(a.2.cmp(&b.2))
                .then(a.3.cmp(&b.3))
        });

        for key in sorted_keys {
            let count = intron_counts[key];
            let hint = AugustusHint::intron_hint(
                key.0.clone(),
                key.1,
                key.2,
                key.3,
                self.config.source.clone(),
                self.config.priority,
                count,
            );
            compressed.push(hint);
        }

        compressed
    }

    /// Compress other hint types (exonpart, exon, dss, ass) with multiplicity counting
    fn compress_other_hints(&mut self) -> Vec<AugustusHint> {
        let mut all_other_hints = Vec::new();

        // Process exonpart hints with consolidation and multiplicity
        let consolidated_exonparts = self.consolidate_exonpart_hints_with_multiplicity();
        all_other_hints.extend(consolidated_exonparts);

        // Process exon hints with multiplicity counting
        let compressed_exons = self.compress_exon_hints_with_multiplicity();
        all_other_hints.extend(compressed_exons);

        // Process splice site hints with multiplicity counting
        let compressed_dss = Self::compress_splice_site_hints_static(&mut self.hint_lists.dss_hints, &self.config);
        all_other_hints.extend(compressed_dss);

        let compressed_ass = Self::compress_splice_site_hints_static(&mut self.hint_lists.ass_hints, &self.config);
        all_other_hints.extend(compressed_ass);

        // Sort for deterministic output
        all_other_hints.sort_by(|a, b| {
            a.chromosome.cmp(&b.chromosome)
                .then(a.start.cmp(&b.start))
                .then(a.end.cmp(&b.end))
                .then(a.hint_type.to_string().cmp(&b.hint_type.to_string()))
                .then(a.strand.cmp(&b.strand))
        });

        all_other_hints
    }

    /// Compress exon hints with multiplicity counting
    fn compress_exon_hints_with_multiplicity(&mut self) -> Vec<AugustusHint> {
        use std::collections::HashMap;

        let mut exon_counts: HashMap<(String, u64, u64, Strand), u32> = HashMap::new();

        // Count occurrences of each exon
        for hint in self.hint_lists.exon_hints.drain(..) {
            let key = (hint.chromosome, hint.start, hint.end, hint.strand);
            *exon_counts.entry(key).or_insert(0) += 1;
        }

        // Create compressed hints with multiplicity
        let mut compressed = Vec::new();
        for ((chromosome, start, end, strand), count) in exon_counts {
            let mut hint = AugustusHint::new(
                chromosome,
                start,
                end,
                HintType::Exon,
                strand,
                self.config.source.clone(),
                self.config.priority,
            );
            hint.multiplicity = count;
            compressed.push(hint);
        }

        compressed
    }

    /// Compress splice site hints with multiplicity counting
    fn compress_splice_site_hints_static(hints: &mut Vec<AugustusHint>, config: &Bam2HintsConfig) -> Vec<AugustusHint> {
        use std::collections::HashMap;

        let mut site_counts: HashMap<(String, u64, u64, HintType, Strand), u32> = HashMap::new();

        // Count occurrences of each splice site
        for hint in hints.drain(..) {
            let key = (hint.chromosome, hint.start, hint.end, hint.hint_type, hint.strand);
            *site_counts.entry(key).or_insert(0) += 1;
        }

        // Create compressed hints with multiplicity
        let mut compressed = Vec::new();
        for ((chromosome, start, end, hint_type, strand), count) in site_counts {
            let mut hint = AugustusHint::new(
                chromosome,
                start,
                end,
                hint_type,
                strand,
                config.source.clone(),
                config.priority,
            );
            hint.multiplicity = count;
            compressed.push(hint);
        }

        compressed
    }

    /// Consolidate exonpart hints with multiplicity - merge overlapping intervals and count support
    fn consolidate_exonpart_hints_with_multiplicity(&mut self) -> Vec<AugustusHint> {
        let exonpart_hints = self.hint_lists.exonpart_hints.drain(..).collect::<Vec<_>>();

        if exonpart_hints.is_empty() {
            return Vec::new();
        }

        // Group by chromosome and strand
        use std::collections::HashMap;
        let mut groups: HashMap<(String, Strand), Vec<AugustusHint>> = HashMap::new();

        for hint in exonpart_hints {
            let key = (hint.chromosome.clone(), hint.strand);
            groups.entry(key).or_insert_with(Vec::new).push(hint);
        }

        let mut consolidated = Vec::new();

        for ((chromosome, strand), mut hints) in groups {
            // Sort by start position
            hints.sort_by_key(|h| h.start);

            let mut merged = Vec::new();

            for hint in hints {
                if merged.is_empty() {
                    merged.push(hint);
                } else {
                    let last_idx = merged.len() - 1;
                    let last = &merged[last_idx];

                    // Check if current hint overlaps with the last merged hint
                    if hint.start <= last.end + 1 {
                        // Merge by extending the end coordinate and adding multiplicity
                        let new_end = std::cmp::max(last.end, hint.end);
                        let mut merged_hint = AugustusHint::new(
                            chromosome.clone(),
                            last.start,
                            new_end,
                            HintType::Exonpart,
                            strand,
                            last.source.clone(),
                            last.priority,
                        );
                        merged_hint.multiplicity = last.multiplicity + hint.multiplicity;
                        merged[last_idx] = merged_hint;
                    } else {
                        // No overlap, add as new hint
                        merged.push(hint);
                    }
                }
            }

            consolidated.extend(merged);
        }

        consolidated
    }



    /// Compress intron hints merging strand-specific introns (for unstranded libraries)
    fn compress_intron_hints_strand_neutral(&mut self) -> Vec<AugustusHint> {
        let mut intron_counts: HashMap<(String, u64, u64), u32> = HashMap::new();

        // Count occurrences of intron hints, ignoring strand
        for hint in &self.hint_lists.intron_hints {
            let key = (
                hint.chromosome.clone(),
                hint.start,
                hint.end,
            );
            *intron_counts.entry(key).or_insert(0) += 1;
        }

        // Create compressed hints with strand set to Unknown (strand-neutral)
        let mut compressed = Vec::new();
        let mut sorted_keys: Vec<_> = intron_counts.keys().collect();
        sorted_keys.sort_by(|a, b| {
            // Sort by chromosome, start, end
            a.0.cmp(&b.0)
                .then(a.1.cmp(&b.1))
                .then(a.2.cmp(&b.2))
        });

        for key in sorted_keys {
            let count = intron_counts[key];
            let hint = AugustusHint::intron_hint(
                key.0.clone(),
                key.1,
                key.2,
                Strand::Unknown, // Strand-neutral for merged introns
                self.config.source.clone(),
                self.config.priority,
                count,
            );
            compressed.push(hint);
        }

        compressed
    }

    /// Clear all hint lists
    pub fn clear(&mut self) {
        self.hint_lists = HintLists::default();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::{RnaSeqAlignment, SpliceJunction, Strand};

    fn create_test_alignment() -> RnaSeqAlignment {
        RnaSeqAlignment {
            read_name: "test_read".to_string(),
            chromosome: "chr1".to_string(),
            start: 1000,
            end: 2000,
            strand: Strand::Forward,
            cigar: "100M900N100M".to_string(), // 100bp match, 900bp intron, 100bp match
            mapping_quality: 60,
            splice_junctions: vec![SpliceJunction {
                chromosome: "chr1".to_string(),
                donor_pos: 1099,
                acceptor_pos: 2000,
                strand: Strand::Forward,
                support_count: 1,
            }],
            is_paired: false,
            is_first_in_pair: false,
        }
    }

    #[test]
    fn test_cigar_parsing() {
        let config = Bam2HintsConfig::default();
        let converter = Bam2HintsConverter::new(config);

        let cigar = "100M900N100M";
        let ops = converter.parse_cigar_string(cigar).unwrap();

        assert_eq!(ops.len(), 3);
        assert_eq!(ops[0], ('M', 100));
        assert_eq!(ops[1], ('N', 900));
        assert_eq!(ops[2], ('M', 100));
    }

    #[test]
    fn test_alignment_to_blocks() {
        let config = Bam2HintsConfig::default();
        let converter = Bam2HintsConverter::new(config);
        let alignment = create_test_alignment();

        let blocks = converter.parse_cigar_to_blocks(&alignment).unwrap();

        assert_eq!(blocks.len(), 2);
        assert_eq!(blocks[0].target_start, 1000);
        assert_eq!(blocks[0].length, 100);
        assert_eq!(blocks[1].target_start, 2000);
        assert_eq!(blocks[1].length, 100);
    }

    #[test]
    fn test_block_filtering() {
        let config = Bam2HintsConfig::default();
        let converter = Bam2HintsConverter::new(config);

        let blocks = vec![
            AlignmentBlock {
                length: 100,
                query_start: 1,
                target_start: 1000,
            },
            AlignmentBlock {
                length: 100,
                query_start: 101,
                target_start: 2000,
            },
        ];

        let filtered = converter.filter_blocks(&blocks).unwrap();

        assert_eq!(filtered.len(), 2);
        assert_eq!(filtered[0].start, 1000);
        assert_eq!(filtered[0].end, 1099);
        assert_eq!(filtered[1].start, 2000);
        assert_eq!(filtered[1].end, 2099);
    }

    #[test]
    fn test_hint_generation() {
        let config = Bam2HintsConfig::default();
        let mut converter = Bam2HintsConverter::new(config);
        let alignment = create_test_alignment();

        converter.process_alignment(&alignment).unwrap();
        let hints = converter.get_hints();

        // Should generate at least one intron hint
        assert!(!hints.is_empty());
        let intron_hints: Vec<_> = hints.iter()
            .filter(|h| h.hint_type == HintType::Intron)
            .collect();
        assert!(!intron_hints.is_empty());
    }

    #[test]
    fn test_multiplicity_compression() {
        let config = Bam2HintsConfig::default();
        let mut converter = Bam2HintsConverter::new(config);
        let alignment = create_test_alignment();

        // Process the same alignment multiple times
        for _ in 0..3 {
            converter.process_alignment(&alignment).unwrap();
        }

        let hints = converter.get_hints();
        let intron_hints: Vec<_> = hints.iter()
            .filter(|h| h.hint_type == HintType::Intron)
            .collect();

        // Should have compressed identical intron hints
        assert_eq!(intron_hints.len(), 1);
        assert_eq!(intron_hints[0].multiplicity, 3);
    }

    #[test]
    fn test_augustus_hint_gff_format() {
        let hint = AugustusHint::intron_hint(
            "chr1".to_string(),
            1000,
            2000,
            Strand::Forward,
            "E".to_string(),
            4,
            5,
        );

        let gff_line = hint.to_gff3_line();
        assert_eq!(gff_line, "chr1\tb2h\tintron\t1000\t2000\t0\t+\t.\tmult=5;pri=4;src=E");
    }

    #[test]
    fn test_config_defaults() {
        let config = Bam2HintsConfig::default();

        assert_eq!(config.priority, 4);
        assert_eq!(config.max_gap_len, 14);
        assert_eq!(config.min_intron_len, 32);
        assert_eq!(config.max_intron_len, 350000);
        assert_eq!(config.source, "E");
        assert!(config.introns_only); // Default is true for Rust (CLI default)
        assert!(!config.no_multiplicity);
    }
}
