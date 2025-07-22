import sys
import pandas as pd
import statsmodels.api as sm
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns

input_file = sys.argv[1]
output_file = sys.argv[2]

df = pd.read_csv(input_file, sep="\t", header=None,
                 names=["chr", "start", "end", "gc", "hypo", "total"],
                 na_values="NA")  

# Prepare empty columns for output
df["resid_total"] = np.nan
df["resid_hypo"] = np.nan
df["resid_fraction"] = np.nan

# Linear regression: GC vs total and hypo
valid = df[["gc", "hypo", "total"]].dropna().index #only correct valid fractions
X = sm.add_constant(df.loc[valid, "gc"])

model_total = sm.OLS(df.loc[valid, "total"], X).fit()
model_hypo = sm.OLS(df.loc[valid, "hypo"], X).fit()

# Residuals
# Fill residuals only for valid rows
df.loc[valid, "resid_total"] = model_total.resid
df.loc[valid, "resid_hypo"] = model_hypo.resid

# Calculate corrected fraction
df.loc[valid, "resid_fraction"] = df.loc[valid, "resid_hypo"] / (df.loc[valid, "resid_total"] + 1e-6)

# Save corrected table
df[["chr", "start", "end","resid_fraction"]].to_csv(
    output_file, sep="\t", index=False, header=False, na_rep="NA"
)

# ----------------- PLOTTING -------------------

# Original fraction (hypo / total)
df["original_fraction"] = df["hypo"] / df["total"]

# Filter valid data for plotting (avoid divide by zero or NaNs)
plot_valid = df.dropna(subset=["gc", "original_fraction", "resid_fraction"])

plt.figure(figsize=(12, 6))

# Plot original fraction vs GC content
plt.subplot(1, 2, 1)
sns.scatterplot(x="gc", y="original_fraction", data=plot_valid, alpha=0.4, s=20)
sns.regplot(x="gc", y="original_fraction", data=plot_valid, scatter=False, lowess=True, color='r')
plt.title("Original Fraction vs GC Content")
plt.xlabel("GC Content")
plt.ylabel("Hypo / Total Fraction")

# Plot corrected fraction vs GC content
plt.subplot(1, 2, 2)
sns.scatterplot(x="gc", y="resid_fraction", data=plot_valid, alpha=0.4, s=20)
sns.regplot(x="gc", y="resid_fraction", data=plot_valid, scatter=False, lowess=True, color='r')
plt.title("Corrected Fraction vs GC Content")
plt.xlabel("GC Content")
plt.ylabel("Corrected Fraction (Residual)")

# saving
# Get base filename (without extension)
base_name = os.path.splitext(os.path.basename(output_file))[0]

# Make plots subdirectory under the output directory
plot_dir = os.path.join(os.path.dirname(output_file), "plots")
os.makedirs(plot_dir, exist_ok=True)

# Save plot
plot_file = os.path.join(plot_dir, f"{base_name}_gc_plot.png")
plt.tight_layout()
plt.savefig(plot_file, dpi=300)