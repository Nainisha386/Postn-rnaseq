"""Aggregates kallisto transcript output to gene level.
Writes raw counts (input to DESeq2) and TPM (for reading expression levels)."""
import os, pandas as pd
home = os.path.expanduser("~")
t2g = pd.read_csv(home + "/Downloads/t2g.txt", sep="\t", header=None,
                  usecols=[0, 1, 2], names=["tx", "gene_id", "symbol"])
qdir = home + "/postn_rnaseq/quant"
counts, tpms = {}, {}
for s in sorted(os.listdir(qdir)):
    f = os.path.join(qdir, s, "abundance.tsv")
    if not os.path.exists(f): continue
    a = pd.read_csv(f, sep="\t")
    m = a.merge(t2g, left_on="target_id", right_on="tx", how="left")
    g = m.groupby(["gene_id", "symbol"], dropna=False)[["est_counts", "tpm"]].sum()
    counts[s], tpms[s] = g["est_counts"], g["tpm"]
pd.DataFrame(counts).to_csv(home + "/postn_rnaseq/gene_counts.csv")
pd.DataFrame(tpms).to_csv(home + "/postn_rnaseq/gene_tpm.csv")
print("genes:", len(counts[list(counts)[0]]), " libraries:", len(counts))
