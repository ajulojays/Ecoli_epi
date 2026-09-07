#!/usr/bin/env python3
"""
script_13_build_master.py

Builds the master AMRFinderPlus report across 704 E. coli target assemblies.
Includes high-precision nanosecond timing, NaN data cleaning, and 
metadata/platform integration.
"""

import os
import sys
import time
import pandas as pd

# ==========================================
# 1. NANOSECOND TIMING INITIALIZATION
# ==========================================
start_time_ns = time.time_ns()
timestamp_str = time.strftime('%Y-%m-%d %H:%M:%S')
print(f"[{timestamp_str}] [INFO] Starting AMRFinderPlus Master Report Build...")

# ==========================================
# 2. PATH CONFIGURATION
# ==========================================
AMR_DIR = "/home/samuelajulo/epi/Ecoli_epi/results/13_amrfinderplus"
MANIFEST = os.path.join(AMR_DIR, "assemblies_manifest.tsv")
OUTPUT_FILE = os.path.join(AMR_DIR, "master_AMRFinder_report.tsv")
TARGET_METADATA = os.path.join(AMR_DIR, "target_metadata_beef_cattle_unspecified_preharvest.csv")

# ==========================================
# 3. VERIFY MANIFEST AND SAMPLES
# ==========================================
if not os.path.exists(MANIFEST):
    print(f"[ERROR] Manifest file not found: {MANIFEST}")
    sys.exit(1)

target_samples = []
with open(MANIFEST, "r") as f:
    for line in f:
        if line.strip():
            sample_id = line.strip().split("\t")[0]
            target_samples.append(sample_id)

total_samples = len(target_samples)
print(f"[INFO] Successfully loaded {total_samples} target sample IDs from manifest.")

# ==========================================
# 4. STREAM CONCATENATION (PREVENTS DEADLOCKS)
# ==========================================
written_rows = 0
found_files = 0
header_written = False

with open(OUTPUT_FILE, "w") as out_h:
    for sample_id in target_samples:
        fpath = os.path.join(AMR_DIR, f"{sample_id}.tsv")
        
        if not os.path.exists(fpath):
            fpath_alt = os.path.join(AMR_DIR, f"{sample_id}_amrfinder.tsv")
            if os.path.exists(fpath_alt):
                fpath = fpath_alt
            else:
                print(f"[WARN] Output TSV missing for sample: {sample_id}")
                continue

        found_files += 1
        with open(fpath, "r") as in_h:
            lines = in_h.readlines()
            if not lines:
                continue

            # Write header once with 'Sample' prepended
            if not header_written:
                header = lines[0].strip()
                out_h.write(f"Sample\t{header}\n")
                header_written = True

            # Write data rows
            for line in lines[1:]:
                line_clean = line.strip()
                if line_clean:
                    out_h.write(f"{sample_id}\t{line_clean}\n")
                    written_rows += 1

print(f"[SUCCESS] Aggregated {found_files} / {total_samples} individual sample files.")
print(f"[INFO] Raw hit records concatenated: {written_rows}")

# ==========================================
# 5. DATA SANITIZATION (NaN HANDLING & METADATA)
# ==========================================
if os.path.exists(OUTPUT_FILE) and written_rows > 0:
    print("[INFO] Sanitizing NaN fields and merging target metadata...")
    df = pd.read_csv(OUTPUT_FILE, sep="\t")
    
    # Replace empty/NaN values with 'NA' for standardized downstream R/Python parsing
    df.fillna("NA", inplace=True)

    # Attach cohort metadata (e.g. year, sequencing platform/Nano tags) if present
    if os.path.exists(TARGET_METADATA):
        meta = pd.read_csv(TARGET_METADATA)
        sample_col = df.columns[0]
        meta_id_col = meta.columns[0]
        
        df[sample_col] = df[sample_col].astype(str).str.strip()
        meta[meta_id_col] = meta[meta_id_col].astype(str).str.strip()
        
        df = df.merge(meta, left_on=sample_col, right_on=meta_id_col, how="left")
        print(f"[INFO] Metadata successfully merged ({len(meta.columns)} attributes attached).")

    # Overwrite master report with clean formatted table
    df.to_csv(OUTPUT_FILE, sep="\t", index=False)

# ==========================================
# 6. NANOSECOND EXECUTION TIMING SUMMARY
# ==========================================
end_time_ns = time.time_ns()
elapsed_ns = end_time_ns - start_time_ns
elapsed_ms = elapsed_ns / 1_000_000
elapsed_sec = elapsed_ns / 1_000_000_000

print("\n----------------------------------------------------------")
print(" ✅ Master Report Build Completed Successfully")
print(f" Master File Output : {OUTPUT_FILE}")
print(f" Total Samples      : {found_files} / {total_samples}")
print(f" Total Hits         : {written_rows}")
print(f" Execution Time     : {elapsed_ns:,} ns ({elapsed_ms:.2f} ms / {elapsed_sec:.3f} s)")
print("----------------------------------------------------------")
