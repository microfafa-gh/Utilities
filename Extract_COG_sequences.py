#!/usr/bin/env python3

import os
from pathlib import Path
from collections import defaultdict
from Bio import SeqIO
import re
import sys
import csv

# === PATHS ===
prokka_root = Path("./prokka")
eggnog_root = Path("./eggnog")
fna_root = Path("./genomes_selected")
output_root_nt = Path("./cog_sequences_nt")
output_root_aa = Path("./cog_sequences_aa")
log_file = Path("./cog_extraction_log.txt")
cog_list_file = Path("./cog_gene_list.txt")  # expects tab-separated file with 'cog' and optional 'gene' columns
genome_list_file = Path("./genome_list.txt")  # optional: list of genome IDs to include

# === LOGGING ===
log_file.write_text("=== COG Sequence Extraction Log ===\n")
def log(msg):
    print(msg)
    log_file.write_text(log_file.read_text() + msg + "\n")

# === OVERWRITE POLICY ===
overwrite_policy = "ask"  # Options: 'ask', 'always', 'never', 'dry-run'
if len(sys.argv) > 1:
    arg = sys.argv[1].lower()
    if arg in ["ask", "always", "never", "dry-run"]:
        overwrite_policy = arg
    else:
        log("[ERROR] Invalid overwrite policy. Use: ask, always, never, or dry-run.")
        sys.exit(1)

is_dry_run = overwrite_policy == "dry-run"

# === LOAD COG LIST FROM FILE ===
cog_to_gene = {}
with open(cog_list_file) as f:
    reader = csv.DictReader(f, delimiter="\t")
    if "cog" not in reader.fieldnames:
        log("[ERROR] 'cog' column missing in COG list file")
        sys.exit(1)
    for row in reader:
        cog = row["cog"].strip()
        gene = row["gene"].strip() if "gene" in row and row["gene"].strip() else cog
        if cog:
            cog_to_gene[cog] = gene

if not cog_to_gene:
    log("[ERROR] No valid COGs loaded from list")
    sys.exit(1)

# === GENOME LIST ===
if genome_list_file.exists():
    with open(genome_list_file) as f:
        genome_ids = [line.strip() for line in f if line.strip()]
else:
    genome_ids = [folder.name for folder in prokka_root.iterdir() if folder.is_dir()]

total_genomes = len(genome_ids)

# Map genome_id -> {locus_tag -> COG}
genome_cog_hits = defaultdict(dict)
parsed_count = 0

# === STEP 1: Parse eggNOG .annotations files ===
log("\n[Step 1] Parsing eggNOG annotations...")
for genome_id in genome_ids:
    pattern = re.compile(rf"{re.escape(genome_id)}(\.\d+)?\.emapper\.annotations")
    matched_files = [f for f in eggnog_root.glob("*.annotations") if pattern.match(f.name)]
    if not matched_files:
        log(f"[WARN] No annotations found for genome {genome_id}")
        continue
    eggnog_file = matched_files[0]
    with open(eggnog_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 6:
                continue
            locus_tag = parts[0]
            cog_field = parts[4]
            if cog_field == "-":
                continue
            cogs = [c.split("@")[0] for c in cog_field.split(",")]
            for cog in cogs:
                if cog in cog_to_gene:
                    genome_cog_hits[genome_id][locus_tag] = cog

# === STEP 2: Extract sequences ===
log("\n[Step 2] Extracting gene sequences...")
for i, genome_id in enumerate(genome_ids, 1):
    log(f"[INFO] Processing genome {i}/{total_genomes}: {genome_id}")

    gff_path = prokka_root / genome_id / f"{genome_id}.gff"
    faa_path = prokka_root / genome_id / f"{genome_id}.faa"
    fna_candidates = list(fna_root.glob(f"{genome_id}*.fna"))

    if not gff_path.exists() or not faa_path.exists() or not fna_candidates:
        log(f"[SKIP] Missing files for {genome_id} (gff/faa/fna)")
        continue

    if genome_id not in genome_cog_hits:
        log(f"[SKIP] No COG hits found for {genome_id} in annotations")
        continue

    fna_path = fna_candidates[0]
    parsed_count += 1

    genome_seqs_nt = SeqIO.to_dict(SeqIO.parse(fna_path, "fasta"))
    genome_seqs_aa = SeqIO.to_dict(SeqIO.parse(faa_path, "fasta"))

    gff_coords = {}
    with open(gff_path) as gff:
        for line in gff:
            if line.startswith("#") or "\tCDS\t" not in line:
                continue
            parts = line.strip().split("\t")
            seq_id, _, _, start, end, _, strand, _, attr = parts
            locus_tag = None
            for item in attr.split(";"):
                if item.startswith("locus_tag="):
                    locus_tag = item.split("=")[1]
                    break
            if locus_tag:
                gff_coords[locus_tag] = (seq_id, int(start), int(end), strand)

    for locus_tag, cog in genome_cog_hits[genome_id].items():
        if locus_tag not in gff_coords:
            continue

        seq_id, start, end, strand = gff_coords[locus_tag]
        seq_nt = genome_seqs_nt[seq_id].seq[start - 1:end]
        if strand == "-":
            seq_nt = seq_nt.reverse_complement()

        header = f"{genome_id}|{locus_tag}|{cog}|{cog_to_gene[cog]}"
        gene_name = cog_to_gene[cog].replace(" ", "_")

        out_nt_dir = output_root_nt / f"{cog}_{gene_name}"
        out_nt_dir.mkdir(parents=True, exist_ok=True)
        out_nt_file = out_nt_dir / f"{genome_id}_{cog}.fna"

        if out_nt_file.exists():
            if overwrite_policy == "never":
                log(f"[SKIP] Existing file (nt) {out_nt_file}, not overwritten.")
                continue
            elif overwrite_policy == "ask":
                response = input(f"[ASK] Overwrite {out_nt_file}? (y/n): ")
                if response.lower() != "y":
                    continue

        if is_dry_run:
            log(f"[DRY-RUN] Would write {out_nt_file}")
        else:
            with open(out_nt_file, "w") as f:
                f.write(f">{header}\n{seq_nt}\n")

        matched_aa = None
        for record in genome_seqs_aa.values():
            if locus_tag in record.id:
                matched_aa = record.seq
                break

        out_aa_dir = output_root_aa / f"{cog}_{gene_name}"
        out_aa_dir.mkdir(parents=True, exist_ok=True)
        out_aa_file = out_aa_dir / f"{genome_id}_{cog}.faa"

        if matched_aa:
            if out_aa_file.exists():
                if overwrite_policy == "never":
                    log(f"[SKIP] Existing file (aa) {out_aa_file}, not overwritten.")
                    continue
                elif overwrite_policy == "ask":
                    response = input(f"[ASK] Overwrite {out_aa_file}? (y/n): ")
                    if response.lower() != "y":
                        continue

            if is_dry_run:
                log(f"[DRY-RUN] Would write {out_aa_file}")
            else:
                with open(out_aa_file, "w") as f:
                    f.write(f">{header}\n{matched_aa}\n")

log(f"\n Processed {parsed_count} genomes with valid input and COG annotations.")
