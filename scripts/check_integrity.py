#!/usr/bin/env python3
"""Check the versioned publication assets and exact unpacked source bytes.

Copyright 2026 Rohit Kumar Jha
SPDX-License-Identifier: Apache-2.0
"""
import hashlib
import json
from pathlib import Path, PurePosixPath
import sys
import zipfile

ROOT = Path(__file__).resolve().parent.parent


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def check_inventory(files, inventory, prefix=""):
    listed = set()
    for line in files[inventory].decode("utf-8").splitlines():
        expected, relative = line.split("  ", 1)
        target = prefix + relative
        require(target in files and target != inventory, f"Invalid entry: {target}")
        require(target not in listed, f"Duplicate entry: {target}")
        require(expected == sha256(files[target]), f"Hash mismatch: {target}")
        listed.add(target)
    require(
        listed == {n for n in files if n.startswith(prefix) and n != inventory},
        f"Incomplete checksum coverage: {inventory}",
    )


def main():
    manifest = json.loads((ROOT / "release-manifest.json").read_text())
    slug = manifest["paper"]
    release = ROOT / "release"
    files = {p.name: p.read_bytes() for p in release.iterdir() if p.is_file()}
    require(set(files) == set(manifest["files"]), "Release file list changed")
    for name, expected in manifest["files"].items():
        require(sha256(files[name]) == expected, f"Release changed: {name}")
    check_inventory(files, "SHA256SUMS.txt")
    archive_path = release / f"{slug}-v1-sources.zip"
    with zipfile.ZipFile(archive_path) as archive:
        require(archive.testzip() is None, "Source archive CRC failure")
        names = archive.namelist()
        require(len(names) == len(set(names)), "Duplicate archive members")
        for name in names:
            path = PurePosixPath(name)
            require(not path.is_absolute() and ".." not in path.parts and "\\" not in name,
                    f"Unsafe archive path: {name}")
            require(not any(part in {".git", ".lake", "__MACOSX", ".DS_Store"}
                            for part in path.parts), f"Unwanted archive member: {name}")
        sources = {name: archive.read(name) for name in names}
    check_inventory(sources, "SHA256SUMS.txt")
    check_inventory(sources, "supplement/SHA256SUMS.txt", "supplement/")
    for name, data in sources.items():
        path = ROOT / "sources" / name
        require(path.is_file() and not path.is_symlink(), f"Missing regular source file: {name}")
        require(path.read_bytes() == data, f"Unpacked source differs: {name}")
    require(len(sources) == manifest["source_members"], "Source member count changed")
    require((ROOT / "LICENSES/Apache-2.0.txt").read_bytes() ==
            sources["LICENSES/Apache-2.0.txt"], "Apache license copy changed")
    print(f"PASS: {slug}: {len(files)} release files; "
          f"{len(sources)} source members; complete hashes and exact source identity.")
    print("This is an integrity check, not a run of the exact audits or a verification of the analytic proof.")


if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError, KeyError, zipfile.BadZipFile) as error:
        print(f"FAIL: {error}", file=sys.stderr)
        sys.exit(1)
