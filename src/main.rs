use anyhow::Result;
use clap::{Parser, Subcommand};
use log::info;

mod commands;

use commands::bam2hints::Bam2HintsCommand;
use commands::utrs::UtrsCommand;

/// AnnoRefine: Genome annotation refinement toolkit
#[derive(Parser)]
#[command(name = "annorefine")]
#[command(about = "Genome annotation refinement toolkit using RNA-seq data")]
#[command(version)]
struct Args {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Extend and refine UTRs using RNA-seq data
    Utrs(UtrsCommand),
    /// Convert BAM alignments to Augustus hints
    Bam2hints(Bam2HintsCommand),
}

fn main() -> Result<()> {
    let args = Args::parse();

    match args.command {
        Commands::Utrs(cmd) => {
            info!("Starting UTR extension and refinement");
            cmd.run()
        }
        Commands::Bam2hints(cmd) => {
            info!("Starting BAM to hints conversion");
            cmd.run()
        }
    }
}
