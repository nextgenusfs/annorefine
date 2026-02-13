#!/usr/bin/env python3
"""
Sync version between Cargo.toml and pyproject.toml.

Usage:
    python sync_version.py                  # Show current versions
    python sync_version.py 2026.2.13        # Set version in both files
    python sync_version.py --check          # Check if versions match (exit 1 if not)
"""
import re
import sys
from pathlib import Path


def get_cargo_version():
    """Extract version from Cargo.toml."""
    cargo_path = Path("Cargo.toml")
    content = cargo_path.read_text()
    match = re.search(r'^version\s*=\s*"([^"]+)"', content, re.MULTILINE)
    if match:
        return match.group(1)
    return None


def get_pyproject_version():
    """Extract version from pyproject.toml."""
    pyproject_path = Path("pyproject.toml")
    content = pyproject_path.read_text()
    match = re.search(r'^version\s*=\s*"([^"]+)"', content, re.MULTILINE)
    if match:
        return match.group(1)
    return None


def set_cargo_version(version):
    """Set version in Cargo.toml."""
    cargo_path = Path("Cargo.toml")
    content = cargo_path.read_text()
    new_content = re.sub(
        r'^version\s*=\s*"[^"]+"',
        f'version = "{version}"',
        content,
        count=1,
        flags=re.MULTILINE
    )
    cargo_path.write_text(new_content)
    print(f"✓ Updated Cargo.toml to version {version}")


def set_pyproject_version(version):
    """Set version in pyproject.toml."""
    pyproject_path = Path("pyproject.toml")
    content = pyproject_path.read_text()
    new_content = re.sub(
        r'^version\s*=\s*"[^"]+"',
        f'version = "{version}"',
        content,
        count=1,
        flags=re.MULTILINE
    )
    pyproject_path.write_text(new_content)
    print(f"✓ Updated pyproject.toml to version {version}")


def validate_version(version):
    """Validate version format (CalVer: YYYY.M.D or YYYY.M.D.PATCH)."""
    pattern = r'^\d{4}\.\d{1,2}\.\d{1,2}(\.\d+)?$'
    if not re.match(pattern, version):
        print(f"Error: Invalid version format '{version}'", file=sys.stderr)
        print("Expected format: YYYY.M.D (e.g., 2026.2.13)", file=sys.stderr)
        sys.exit(1)


def main():
    if len(sys.argv) == 1:
        # Show current versions
        cargo_ver = get_cargo_version()
        pyproject_ver = get_pyproject_version()
        
        print("Current versions:")
        print(f"  Cargo.toml:     {cargo_ver}")
        print(f"  pyproject.toml: {pyproject_ver}")
        
        if cargo_ver == pyproject_ver:
            print("\n✓ Versions are in sync")
        else:
            print("\n✗ Versions are OUT OF SYNC!")
            print("\nTo sync versions, run:")
            print(f"  python sync_version.py {cargo_ver or pyproject_ver}")
            sys.exit(1)
    
    elif sys.argv[1] == "--check":
        # Check if versions match
        cargo_ver = get_cargo_version()
        pyproject_ver = get_pyproject_version()
        
        if cargo_ver == pyproject_ver:
            print(f"✓ Versions match: {cargo_ver}")
            sys.exit(0)
        else:
            print(f"✗ Version mismatch!", file=sys.stderr)
            print(f"  Cargo.toml:     {cargo_ver}", file=sys.stderr)
            print(f"  pyproject.toml: {pyproject_ver}", file=sys.stderr)
            sys.exit(1)
    
    else:
        # Set version in both files
        new_version = sys.argv[1]
        validate_version(new_version)
        
        set_cargo_version(new_version)
        set_pyproject_version(new_version)
        
        print(f"\n✓ Successfully synced version to {new_version}")
        print("\nNext steps:")
        print(f"  git add Cargo.toml pyproject.toml")
        print(f"  git commit -m 'Bump version to {new_version}'")
        print(f"  git tag -a v{new_version} -m 'Release v{new_version}'")
        print(f"  git push origin main v{new_version}")


if __name__ == "__main__":
    main()

