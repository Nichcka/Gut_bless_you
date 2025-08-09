import pandas as pd
import sys

SraRunTable = sys.argv[1]

df = pd.read_csv("3.tsv", sep="\t")

meta = pd.read_csv(SraRunTable, sep="\t")

merged = df.merge(
    meta,
    left_on="assemble",
    right_on="Run",
    how="left"
)

merged.to_csv("final_table.tsv", sep="\t", index=False)

