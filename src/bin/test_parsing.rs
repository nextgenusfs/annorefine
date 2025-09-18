//! Simple test program to verify FASTA and GFF3 parsing

use annorefine::{
    fasta::{get_sequence_stats, parse_fasta_file},
    gff3::parse_gff3_file,
    output::{validate_gene_models, Gff3Writer},
};
use std::env;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    env_logger::init();

    let args: Vec<String> = env::args().collect();
    if args.len() != 3 {
        eprintln!("Usage: {} <fasta_file> <gff3_file>", args[0]);
        std::process::exit(1);
    }

    let fasta_file = &args[1];
    let gff3_file = &args[2];

    println!("Testing FASTA parsing...");
    let genome = parse_fasta_file(fasta_file)?;
    let stats = get_sequence_stats(&genome);
    println!("Genome stats: {}", stats);

    println!("\nTesting GFF3 parsing...");
    let gene_models = parse_gff3_file(gff3_file)?;
    println!("Loaded {} gene models", gene_models.len());

    for gene in &gene_models {
        println!(
            "Gene: {} ({}:{}-{}) with {} transcripts",
            gene.id,
            gene.chromosome,
            gene.start,
            gene.end,
            gene.transcripts.len()
        );

        for transcript in &gene.transcripts {
            println!(
                "  Transcript: {} with {} exons, {} CDS regions",
                transcript.id,
                transcript.exons.len(),
                transcript.cds_regions.len()
            );
        }
    }

    println!("\nValidating gene models...");
    validate_gene_models(&gene_models)?;
    println!("All gene models are valid!");

    println!("\nTesting GFF3 output...");
    let mut writer = Gff3Writer::new("test_output.gff3")?;
    writer.write_gene_models(&gene_models)?;
    println!("Successfully wrote output to test_output.gff3");

    Ok(())
}
