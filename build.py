#!/usr/bin/env python3
"""
Build script to compile the Rust binary and include it in the wheel.
This script is called by maturin during the build process.
"""
import os
import subprocess
import sys
import shutil
from pathlib import Path


def build_binary():
    """Build the Rust binary and copy it to the data/scripts directory."""
    print("Building Rust binary for inclusion in wheel...")
    
    # Get the project root
    project_root = Path(__file__).parent
    
    # Create the data/scripts directory
    scripts_dir = project_root / "annorefine.data" / "scripts"
    scripts_dir.mkdir(parents=True, exist_ok=True)
    
    # Build the binary in release mode
    result = subprocess.run(
        ["cargo", "build", "--release", "--bin", "annorefine"],
        cwd=project_root,
        capture_output=True,
        text=True
    )
    
    if result.returncode != 0:
        print(f"Error building binary: {result.stderr}", file=sys.stderr)
        sys.exit(1)
    
    # Determine the binary name based on platform
    binary_name = "annorefine.exe" if sys.platform == "win32" else "annorefine"
    
    # Copy the binary to the scripts directory
    source_binary = project_root / "target" / "release" / binary_name
    dest_binary = scripts_dir / binary_name
    
    if not source_binary.exists():
        print(f"Error: Binary not found at {source_binary}", file=sys.stderr)
        sys.exit(1)
    
    shutil.copy2(source_binary, dest_binary)
    
    # Make it executable on Unix-like systems
    if sys.platform != "win32":
        os.chmod(dest_binary, 0o755)
    
    print(f"Binary copied to {dest_binary}")


if __name__ == "__main__":
    build_binary()

