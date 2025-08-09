from scipy.stats import wilcoxon
from lists import get_data
import sys

ibs_data = sys.argv[1]
hc_data = sys.argv[2]

common_names = set(ibs_data.keys()) & set(hc_data.keys())

results = []

for name in sorted(common_names):
    vals_ibs = ibs_data[name]
    vals_hc = hc_data[name]

    min_len = min(len(vals_ibs), len(vals_hc))
    if min_len == 0:
        continue
    pairs_ibs = vals_ibs[:min_len]
    pairs_hc = vals_hc[:min_len]

    try:
        _, p_value = wilcoxon(pairs_ibs, pairs_hc)
    except ValueError:

        p_value = 1.0

    results.append((name, p_value))

print("locus\tp-value")
for name, pval in results:
    print(f"{name}\t{pval}")

