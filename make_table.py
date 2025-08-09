import pandas as pd
import sys

df = pd.read_csv("2.tsv", sep="\t")
SQ_loci_meta = sys.argv[1]
meta = pd.read_csv(SQ_loci_meta)

df['locus_trimmed'] = df['locus'].str.replace('_nointergenome_regions', '', regex=False)

merged = df.merge(
    meta,
    left_on='locus_trimmed',
    right_on='cluster',
    how='left'
)

result = merged[[
    'assemble',
    'locus',
    'tpm',
    'log_tpm',
    'Genome',
    'Species',
    'Pathway',
    'Lineage'
]]

result.to_csv("3.tsv", sep="\t", index=False)

