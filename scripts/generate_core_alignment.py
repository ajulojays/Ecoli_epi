#!/usr/bin/env python3
import os
import sys
import csv
import subprocess
from multiprocessing import Pool
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm

# ==============================================================================
# CONFIGURATION
# ==============================================================================
PANAROO_DIR = "/home/samuelajulo/epi/Ecoli_epi/results/panaroo"
PA_FILE = os.path.join(PANAROO_DIR, "gene_presence_absence.csv")
FASTA_FILE = os.path.join(PANAROO_DIR, "combined_DNA_CDS.fasta")
OUT_FILE = os.path.join(PANAROO_DIR, "core_gene_alignment.aln")
THREADS = 64
THRESHOLD = 0.95

# 1. Parse Gene Presence/Absence Matrix
print("[INFO] Parsing gene_presence_absence.csv...")
if not os.path.exists(PA_FILE):
    print(f"[ERROR] Could not find matrix file at {PA_FILE}"); sys.exit(1)

with open(PA_FILE, 'r') as f:
    reader = csv.reader(f)
    header = next(reader)
    
    # FIXED: Panaroo metadata occupies columns 0-13. Samples start strictly at 14.
    samples = header[14:]  
    num_samples = len(samples)
    
    core_genes = []
    for row in reader:
        gene_name = row[0]
        sample_entries = row[14:]
        present_count = sum(1 for x in sample_entries if x.strip())
        if (present_count / num_samples) >= THRESHOLD:
            core_genes.append((gene_name, sample_entries))

print(f"[INFO] Detected {len(core_genes)} core genes across {num_samples} samples (Threshold: {THRESHOLD*100}%).")

# 2. Load Global FASTA Sequences into RAM
print("[INFO] Loading combined_DNA_CDS.fasta into memory...")
if not os.path.exists(FASTA_FILE):
    print(f"[ERROR] Could not find FASTA file at {FASTA_FILE}"); sys.exit(1)
seq_dict = SeqIO.to_dict(SeqIO.parse(FASTA_FILE, "fasta"))
print("[INFO] Sequence database successfully loaded.")

# 3. Worker Function for Threaded MAFFT Processing
def align_gene(gene_data):
    gene_name, sample_entries = gene_data
    unaligned_records = []
    
    for idx, entry in enumerate(sample_entries):
        if entry.strip():
            gene_id = entry.split(';')[0].split('\t')[0].split()[0].strip()
            if gene_id in seq_dict:
                # Safe object duplication to prevent parent memory corruption
                orig_record = seq_dict[gene_id]
                record = SeqRecord(orig_record.seq, id=samples[idx], description="")
                unaligned_records.append(record)
    
    if not unaligned_records:
        return gene_name, {}

    pid = os.getpid()
    tmp_in = f"/dev/shm/tmp_in_{gene_name}_{pid}.fasta"
    tmp_out = f"/dev/shm/tmp_out_{gene_name}_{pid}.fasta"
    
    SeqIO.write(unaligned_records, tmp_in, "fasta")
    
    # Run MAFFT alignment locally using shared memory storage
    subprocess.run(f"mafft --auto --quiet --threads 1 {tmp_in} > {tmp_out}", shell=True)
    
    aligned_dict = {}
    alignment_len = 0
    if os.path.exists(tmp_out) and os.path.getsize(tmp_out) > 0:
        for record in SeqIO.parse(tmp_out, "fasta"):
            aligned_dict[record.id] = str(record.seq)
        if aligned_dict:
            alignment_len = len(next(iter(aligned_dict.values())))
            
    # Clean up temporary RAM disk files
    for temp_file in [tmp_in, tmp_out]:
        if os.path.exists(temp_file):
            os.remove(temp_file)
            
    # Map back alignments and fill missing data slots with gap blocks
    final_gene_alignment = {}
    for s in samples:
        if s in aligned_dict:
            final_gene_alignment[s] = aligned_dict[s]
        else:
            final_gene_alignment[s] = "-" * alignment_len
            
    return gene_name, final_gene_alignment

# 4. Process Multi-Core Executions
print(f"[INFO] Aligning sequence blocks using {THREADS} processing cores...")
sample_alignments = {s: [] for s in samples}

with Pool(THREADS) as pool:
    results = list(tqdm(pool.imap(align_gene, core_genes), total=len(core_genes)))
    
for gene_name, alignment in results:
    if not alignment:
        continue
    for s in samples:
        sample_alignments[s].append(alignment[s])

# 5. Concatenate and Write Outputs
print(f"[INFO] Compiling alignments into final matrix layout: {OUT_FILE}...")
with open(OUT_FILE, 'w') as out_f:
    for s in samples:
        full_seq = "".join(sample_alignments[s])
        out_f.write(f">{s}\n{full_seq}\n")
        
print("[DONE] Core genome alignment file is ready for downstream tree inference.")
