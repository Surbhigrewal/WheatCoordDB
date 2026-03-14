#!/usr/bin/env python3
"""
rename_chromosomes.py
=====================
Renames chromosomes to Chr1A..Chr7D convention in FASTA and GFF3 files.

Handles four input naming styles:
  Ensembl style    :  1A, 1B ... 7D         (CS FASTA + GFF, 10-wheat panel)
  ENA colon style  :  chromosome: 1A        (CS_v1, Kariega, Attraktion,
                                             Renan_v2, Paragon_v3, Cadenza_v2,
                                             Fielder)
  ENA comma style  :  chromosome 1A,        (Aikang58, CS_IAAS, Sumai3,
                                             JIN50, MOV — no colon)
  CS_CAU style     :  chromosome A1 / D6    (subgenome letter before number)

Output naming: Chr1A, Chr1B, Chr1D, Chr2A ... Chr7D

Usage:
  # CS assembly + GFF (Ensembl style)
  python3 rename_chromosomes.py \
      --fasta Triticum_aestivum_refseqv2.IWGSC_RefSeq_v2.1.dna.toplevel.fa.gz \
      --gff   Triticum_aestivum_refseqv2.IWGSC_RefSeq_v2.1.62.chr.gff3.gz \
      --outdir renamed/CS/ \
      --name CS

  # 10-wheat target (ENA accession style, chr name in header description)
  python3 rename_chromosomes.py \
      --fasta GCA_903993795.1_10wheat_assembly_jagger_genomic.fasta \
      --outdir renamed/Jagger/ \
      --name Jagger
"""

import argparse
import gzip
import re
import subprocess
import sys
from pathlib import Path

# ── Target chromosome names ───────────────────────────────────────────────────
WHEAT_CHRS = [f"{g}{s}" for g in range(1, 8) for s in ["A", "B", "D"]]  # 1A..7D


def canonical(raw):
    """
    Given any chromosome name token, return 'Chr1A' style or None.
    Covers:
      1A / 1a / chr1A / Chr1A / CHR1A
      LR* accessions resolved separately via FASTA header description
    """
    m = re.search(r"([1-7][ABDabd])\b", raw)
    if m:
        return f"Chr{m.group(1).upper()}"
    return None


def parse_chr_from_description(description):
    """
    Try all known description styles to extract a canonical Chr name.
    Returns 'Chr1A' string or None.

    Styles handled (tried in order):
      A) 'chromosome: 1A'   — colon style (CS_v1, Kariega, Attraktion,
                               Renan_v2, Paragon_v3, Cadenza_v2, Fielder)
      B) 'chromosome 1A,'   — no colon, comma/space-terminated
                               (Aikang58, CS_IAAS, Sumai3, JIN50, MOV)
      C) 'chromosome A1'    — CS_CAU reversed style (subgenome before number)
    """
    # Style A: colon before standard identifier e.g. 'chromosome: 1A'
    m = re.search(r"chromosome:\s*([1-7][ABDabd])\b", description, re.I)
    if m:
        return f"Chr{m.group(1).upper()}"

    # Style B: 'chromosome' then space then standard identifier (no colon)
    # e.g. 'chromosome 1A,' or 'chromosome 1A '
    m = re.search(r"chromosome\s+([1-7][ABDabd])\b", description, re.I)
    if m:
        return f"Chr{m.group(1).upper()}"

    # Style C: CS_CAU reversed — subgenome letter then number
    # e.g. 'chromosome A1', 'chromosome D6', 'chromosome B3'
    # Must come AFTER Style B to avoid mis-matching 'chromosome 1A' as 'A' + '1'
    m = re.search(r"chromosome\s+([ABDabd])([1-7])\b", description, re.I)
    if m:
        subgenome = m.group(1).upper()
        number = m.group(2)
        return f"Chr{number}{subgenome}"

    return None


def build_mapping_from_fasta_headers(fasta_path):
    """
    Stream FASTA headers and build {old_name: Chr1A} mapping.

    Cases handled:
      A) Ensembl / short style: '>1A dna:chromosome ...'
         seq_id itself matches e.g. '1A' -> 'Chr1A'
      B) ENA colon style: '>OU963872.1 ... chromosome: 1A'
      C) ENA no-colon style: '>CM077322.1 ... chromosome 1A,'
      D) CS_CAU reversed: '>CM100877.1 ... chromosome A1'

    Returns:
      mapping     dict  old_name -> Chr1A
      unplaced    list  (old_name, description) for sequences not mapped
    """
    mapping = {}
    unplaced = []

    opener = gzip.open if str(fasta_path).endswith(".gz") else open

    print(f"  Scanning headers: {Path(fasta_path).name}", file=sys.stderr)
    with opener(fasta_path, "rt") as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            line = line.rstrip()
            parts = line[1:].split(" ", 1)
            seq_id = parts[0]
            description = parts[1] if len(parts) > 1 else ""

            # Case A: seq_id itself is a chromosome name (e.g. "1A")
            chr_name = canonical(seq_id)

            # Cases B/C/D: seq_id is an accession — look in description
            if chr_name is None:
                chr_name = parse_chr_from_description(description)

            if chr_name:
                mapping[seq_id] = chr_name
            else:
                unplaced.append((seq_id, description[:80]))

    return mapping, unplaced


def rename_fasta(fasta_path, out_fasta, mapping):
    """Stream FASTA, rename headers. Unplaced sequences written unchanged."""
    opener_in = gzip.open if str(fasta_path).endswith(".gz") else open

    renamed = 0
    kept = 0

    print(f"  Writing FASTA: {Path(out_fasta).name}", file=sys.stderr)
    with opener_in(fasta_path, "rt") as fh, open(out_fasta, "w") as out:
        for line in fh:
            if line.startswith(">"):
                seq_id = line[1:].split()[0]
                if seq_id in mapping:
                    out.write(f">{mapping[seq_id]}\n")
                    renamed += 1
                else:
                    out.write(line)
                    kept += 1
            else:
                out.write(line)

    print(f"  Renamed: {renamed} chromosomes | kept unchanged: {kept} (unplaced/organelle)",
          file=sys.stderr)
    return renamed


def rename_gff(gff_path, out_gff, mapping):
    """
    Rename seqid (col 1) in GFF3.
    Also renames sequence-region pragma lines:
      ##sequence-region   1A 1 598660471  ->  ##sequence-region   Chr1A 1 598660471
    """
    opener_in = gzip.open if str(gff_path).endswith(".gz") else open

    renamed_features = 0
    renamed_pragmas = 0
    skipped = 0

    print(f"  Writing GFF3:  {Path(out_gff).name}", file=sys.stderr)
    with opener_in(gff_path, "rt") as fh, open(out_gff, "w") as out:
        for line in fh:
            # Handle sequence-region pragma
            if line.startswith("##sequence-region"):
                parts = line.split()
                if len(parts) >= 2:
                    old = parts[1]
                    new = mapping.get(old) or canonical(old)
                    if new:
                        line = line.replace(old, new, 1)
                        renamed_pragmas += 1
                out.write(line)
                continue

            # Skip other comment lines
            if line.startswith("#"):
                out.write(line)
                continue

            # Feature lines: rename col 1
            parts = line.split("\t", 1)
            if len(parts) < 2:
                out.write(line)
                continue

            seqid = parts[0]
            new_seqid = mapping.get(seqid) or canonical(seqid)
            if new_seqid:
                out.write(new_seqid + "\t" + parts[1])
                renamed_features += 1
            else:
                out.write(line)
                skipped += 1

    print(f"  GFF features renamed: {renamed_features} | "
          f"pragmas renamed: {renamed_pragmas} | "
          f"unplaced kept: {skipped}", file=sys.stderr)


def run_samtools_faidx(fasta_path):
    r = subprocess.run(["samtools", "faidx", str(fasta_path)],
                       capture_output=True, text=True)
    if r.returncode != 0:
        print(f"  WARN: samtools faidx failed: {r.stderr.strip()}", file=sys.stderr)
        return False
    return True


def write_chrom_sizes(fasta_path):
    fai = Path(str(fasta_path) + ".fai")
    sizes_out = fasta_path.with_suffix("").with_suffix(".chrom.sizes")
    if not fai.exists():
        print(f"  WARN: .fai not found, skipping chrom.sizes", file=sys.stderr)
        return
    with open(fai) as fh, open(sizes_out, "w") as out:
        for line in fh:
            p = line.split("\t")
            out.write(f"{p[0]}\t{p[1]}\n")
    print(f"  Chrom sizes:   {sizes_out.name}", file=sys.stderr)


def validate(mapping, unplaced, name):
    found = set(mapping.values())
    expected = {f"Chr{c}" for c in WHEAT_CHRS}
    missing = expected - found
    print(f"\n  {name}: {len(found)}/21 chromosomes mapped", file=sys.stderr)
    if missing:
        print(f"  !! MISSING: {sorted(missing)}", file=sys.stderr)
    if unplaced:
        print(f"  Unplaced/organelle sequences: {len(unplaced)}", file=sys.stderr)


def save_mapping_tsv(mapping, unplaced, outdir, name):
    tsv = outdir / f"{name}_chr_mapping.tsv"
    with open(tsv, "w") as fh:
        fh.write("original_id\trenamed_to\n")
        for old, new in sorted(mapping.items(), key=lambda x: x[1]):
            fh.write(f"{old}\t{new}\n")
        if unplaced:
            fh.write("# Unplaced / organelle (not renamed)\n")
            for seq_id, desc in unplaced[:50]:
                fh.write(f"{seq_id}\t# {desc}\n")
    print(f"  Mapping TSV:   {tsv.name}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description="Rename wheat chromosomes to Chr1A..Chr7D in FASTA and GFF3"
    )
    parser.add_argument("--fasta", required=True,
                        help="Input FASTA (.fa/.fna/.fasta, optionally .gz)")
    parser.add_argument("--gff", default=None,
                        help="Input GFF3 (.gff3/.gff, optionally .gz) — optional")
    parser.add_argument("--outdir", required=True,
                        help="Output directory")
    parser.add_argument("--name", required=True,
                        help="Short name used in output filenames (e.g. CS, Jagger)")
    parser.add_argument("--no-index", action="store_true",
                        help="Skip samtools faidx indexing")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"\n{'='*60}", file=sys.stderr)
    print(f"  Assembly: {args.name}", file=sys.stderr)
    print(f"  FASTA in: {args.fasta}", file=sys.stderr)
    if args.gff:
        print(f"  GFF3  in: {args.gff}", file=sys.stderr)
    print(f"  Out dir:  {outdir}", file=sys.stderr)

    # Build mapping from FASTA headers
    mapping, unplaced = build_mapping_from_fasta_headers(args.fasta)
    validate(mapping, unplaced, args.name)
    save_mapping_tsv(mapping, unplaced, outdir, args.name)

    # Rename FASTA
    out_fasta = outdir / f"{args.name}.fa"
    rename_fasta(args.fasta, out_fasta, mapping)

    # Index + chrom sizes
    if not args.no_index:
        print(f"  Running samtools faidx...", file=sys.stderr)
        if run_samtools_faidx(out_fasta):
            write_chrom_sizes(out_fasta)

    # Rename GFF3
    if args.gff:
        out_gff = outdir / f"{args.name}.gff3"
        rename_gff(args.gff, out_gff, mapping)

    print(f"\n  DONE: {args.name}", file=sys.stderr)
    print(f"  Output FASTA: {out_fasta}", file=sys.stderr)
    if args.gff:
        print(f"  Output GFF3:  {outdir / (args.name + '.gff3')}", file=sys.stderr)
    print(f"{'='*60}", file=sys.stderr)


if __name__ == "__main__":
    main()
