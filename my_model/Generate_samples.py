print("Script started", flush=True)
import os
import gzip
import random
import numpy as np
from glob import glob
import argparse

# === PARSE INPUT ===
print("Start parsing arguments etc", flush=True)
parser = argparse.ArgumentParser()
parser.add_argument('--cfDNA_dir', required=True, help='Path to healthy cfDNA background samples')
parser.add_argument('--tissue_dir', required=True, help='Path to tissue samples (cirrhosis + tumour)')
parser.add_argument('--output_dir', required=True, help='Path to preferred output directory')
args = parser.parse_args()

cfDNA_dir = args.cfDNA_dir
tissue_dir = args.tissue_dir
output_dir = args.output_dir

# Define subdirectories and metadata path
OUTPUT_DIR = os.path.join(output_dir, "synthetic_samples")
META_FILE = os.path.join(output_dir, "synthetic_sample_metadata.tsv")

# Create output directory if it doesn't exist
os.makedirs(OUTPUT_DIR, exist_ok=True)

# === CONFIG ===
READS_PER_SAMPLE = 70_000_000
N_HEALTHY = 10
N_LIVER_CONTROL = 5
N_CIRRHOSIS = 100
N_TUMOUR = 100

HEALTHY_LIVER_DIST = lambda: np.random.uniform(0.001, 0.01) # uniform 0.1%-1%
CIRRHOSIS_DIST = lambda: np.random.uniform(0.01, 0.1)  # uniform 1%–10%
TUMOUR_DIST = lambda: np.random.beta(0.3, 6) * 0.15     # skewed < 15%

# === COLLECT FILES ===
print("Start collecting files", flush=True)
cfdna_background_files = glob(os.path.join(cfDNA_dir, "*.bed.gz"))
tissue_files = glob(os.path.join(tissue_dir, "*.bed.gz"))
healthy_liver_files = [f for f in tissue_files if "Cirrhosis" not in f and "Tumour" not in f]
cirrhosis_files = [f for f in tissue_files if "Cirrhosis" in f]
tumour_files = [f for f in tissue_files if "Tumour" in f]

assert cfdna_background_files, "No cfdna background files found"
assert healthy_liver_files, "No healthyliver files found"
assert cirrhosis_files, "No cirrhosis files found"
assert tumour_files, "No tumour files found"

# === DEFINE FUNCTIONS ===
def load_all_reads(file_list):
    """Load all reads from all files in the list, return a list of lines"""
    all_reads = []
    for f in file_list:
        with gzip.open(f, 'rt') as infile:
            all_reads.extend(infile.readlines())
    return all_reads

def sample_reads_from_pool(pool, n):
    """Randomly sample n reads from a pre-loaded pool of reads"""
    return random.sample(pool, n)

def write_sample(reads, out_path):
    with gzip.open(out_path, 'wt') as out:
        for line in reads:
            out.write(line)

# === CREATE SAMPLE POOLS ===
print("Loading all reads for each category", flush=True)
cfdna_background_pool = load_all_reads(cfdna_background_files)
healthy_liver_pool = load_all_reads(healthy_liver_files)
cirrhosis_pool = load_all_reads(cirrhosis_files)
tumour_pool = load_all_reads(tumour_files)
print("Finished loading all reads into memory", flush=True)

# === GENERATE SAMPLES ===
metadata = []

# Healthy-only samples
for i in range(N_HEALTHY):
    reads = sample_reads_from_pool(cfdna_background_pool, READS_PER_SAMPLE)

    out_path = os.path.join(OUTPUT_DIR, f"synthetic_healthy_{i+1:03d}.bed.gz")
    write_sample(reads, out_path)
    metadata.append((os.path.basename(out_path), 'healthy', 1.0, 0.0, 0.0))

# Technical control against overfitting: healthy cfDNA + low fractions of healthy liver tissue
for i in range(N_LIVER_CONTROL):
    frac_liver_tissue = HEALTHY_LIVER_DIST()
    n_liver_tissue = int(READS_PER_SAMPLE * frac_liver_tissue)
    n_healthy = READS_PER_SAMPLE - n_liver_tissue

    reads = sample_reads_from_pool(healthy_liver_pool, n_liver_tissue) + sample_reads_from_pool(cfdna_background_pool, n_healthy)
    random.shuffle(reads)

    out_path = os.path.join(OUTPUT_DIR, f"healthy_liver_control_{i+1:03d}.bed.gz")
    write_sample(reads, out_path)
    metadata.append((os.path.basename(out_path), 'liver_control', 1 - frac_liver_tissue, 0.0, 0.0))


# Cirrhosis-mixed samples
for i in range(N_CIRRHOSIS):
    frac_cirrhosis = CIRRHOSIS_DIST() # sample from the uniform distribution
    n_cirr = int(READS_PER_SAMPLE * frac_cirrhosis)
    n_healthy = READS_PER_SAMPLE - n_cirr

    reads = sample_reads_from_pool(cirrhosis_pool, n_cirr) + sample_reads_from_pool(cfdna_background_pool, n_healthy)
    random.shuffle(reads)

    out_path = os.path.join(OUTPUT_DIR, f"synthetic_cirrhosis_{i+1:03d}.bed.gz")
    write_sample(reads, out_path)
    metadata.append((os.path.basename(out_path), 'cirrhosis', 1 - frac_cirrhosis, frac_cirrhosis, 0.0))

# Tumour-mixed samples
for i in range(N_TUMOUR):
    frac_tumour = TUMOUR_DIST() #sample from beta distribution (more probability to sample lower tumour fractions)
    n_tumour = int(READS_PER_SAMPLE * frac_tumour)
    n_bg = READS_PER_SAMPLE - n_tumour

    include_cirrhosis = random.choice([0, 1])  # 50% chance to include cirrhosis
    if include_cirrhosis:
        frac_cirrhosis = np.random.uniform(0.01, 0.05)  # small fraction (1%–5%)
        n_cirrhosis = int(n_bg * frac_cirrhosis)
        n_healthy = n_bg - n_cirrhosis
    else:
        frac_cirrhosis = 0.0
        n_cirrhosis = 0
        n_healthy = n_bg
    frac_healthy = 1.0 - frac_tumour - frac_cirrhosis
    
    healthy_reads = sample_reads_from_pool(cfdna_background_pool, n_healthy)
    cirrhosis_reads = sample_reads_from_pool(cirrhosis_pool, n_cirrhosis)
    tumour_reads = sample_reads_from_pool(tumour_pool, n_tumour)
    reads = healthy_reads + cirrhosis_reads + tumour_reads
    random.shuffle(reads)

    out_path = os.path.join(OUTPUT_DIR, f"synthetic_tumour_{i+1:03d}.bed.gz")
    write_sample(reads, out_path)
    metadata.append((os.path.basename(out_path), 'tumour', frac_healthy, frac_cirrhosis, frac_tumour)) 

# Write metadata
with open(META_FILE, 'w') as meta:
    meta.write("sample\ttype\thealthy_fraction\tcirrhosis_fraction\ttumour_fraction\n")
    for m in metadata:
        meta.write("\t".join(map(str, m)) + "\n")

print(f"Generated {N_HEALTHY + N_LIVER_CONTROL + N_CIRRHOSIS + N_TUMOUR} synthetic samples in {OUTPUT_DIR}", flush=True)
