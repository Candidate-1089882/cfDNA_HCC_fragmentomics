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
N_HEALTHY = 150
N_TUMOUR = 150


def TUMOUR_DIST():
    min_frac = 0.0001  # 0.01%
    max_frac = 0.15     # 15%
    sample = np.random.beta(a=2, b=8)  # skew toward low values
    scaled = min_frac + sample * (max_frac - min_frac)
    return scaled

# === COLLECT FILES ===
print("Start collecting files", flush=True)
cfdna_background_files = glob(os.path.join(cfDNA_dir, "*.bed.gz"))
tissue_files = glob(os.path.join(tissue_dir, "*.bed.gz"))
tumour_files = [f for f in tissue_files if "Tumour" in f]

assert cfdna_background_files, "No cfdna background files found"
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

# Load tumour reads individually with purity
tumour_purity_map = {
    "CD564934_Liver-Tumour_md.per-read.bed.gz": 0.3637,
    "CD564208_Liver-Tumour_md.per-read.bed.gz": 0.513,
    "CD564146_Liver-Tumour_md.per-read.bed.gz": 0.7285,
    "CD563176_Liver-Tumour_md.per-read.bed.gz": 0.4123,
}

tumour_sample_dict = {}  # {filename: (purity, reads)}
for f in tumour_files:
    fname = os.path.basename(f)
    purity = tumour_purity_map[fname]
    with gzip.open(f, 'rt') as infile:
        reads = infile.readlines()
    tumour_sample_dict[fname] = (purity, reads)
    print(f"Loaded tumour file {fname} with purity {purity} and {len(reads):,} reads", flush=True)


# === GENERATE SAMPLES ===
metadata = []

# Healthy-only samples
for i in range(N_HEALTHY):
    reads = sample_reads_from_pool(cfdna_background_pool, READS_PER_SAMPLE)
    out_path = os.path.join(OUTPUT_DIR, f"synthetic_healthy_{i+1:03d}.bed.gz")
    write_sample(reads, out_path)
    metadata.append((os.path.basename(out_path), 'healthy', 1.0, 0.0))

# Tumour-mixed samples
tumour_sample_items = list(tumour_sample_dict.items())

for i in range(N_TUMOUR):
    target_effective_tumour_fraction = TUMOUR_DIST()
    n_tumour_reads_needed = int(READS_PER_SAMPLE * target_effective_tumour_fraction)

    # Use all tumour samples to sample from
    selected = tumour_sample_items
    purities = [s[1][0] for s in selected]

    # Weighted allocation of reads per sample based on purity
    weights = np.random.dirichlet(np.ones(len(selected)))
    adjusted_weights = [(target_effective_tumour_fraction / purity) * w for purity, w in zip(purities, weights)]

    total_adjusted = sum(adjusted_weights)
    reads_per_sample = [int((w / total_adjusted) * n_tumour_reads_needed) for w in adjusted_weights]

    tumour_reads = []
    for (fname, (purity, reads)), n_reads in zip(selected, reads_per_sample):
        tumour_reads += random.sample(reads, min(n_reads, len(reads)))

    # Remaining reads from healthy background
    n_bg = READS_PER_SAMPLE - len(tumour_reads)
    healthy_reads = sample_reads_from_pool(cfdna_background_pool, n_bg)

    reads = healthy_reads + tumour_reads
    random.shuffle(reads)

    out_path = os.path.join(OUTPUT_DIR, f"synthetic_tumour_{i+1:03d}.bed.gz")
    write_sample(reads, out_path)

    actual_tumour_fraction = len(tumour_reads) / READS_PER_SAMPLE
    metadata.append((os.path.basename(out_path), 'tumour', 1 - actual_tumour_fraction, actual_tumour_fraction))

# Write metadata
with open(META_FILE, 'w') as meta:
    meta.write("sample\ttype\thealthy_fraction\tcirrhosis_fraction\ttumour_fraction\n")
    for m in metadata:
        meta.write("\t".join(map(str, m)) + "\n")

print(f"Generated {N_HEALTHY + N_TUMOUR} synthetic samples in {OUTPUT_DIR}", flush=True)
