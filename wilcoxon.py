from scipy.stats import wilcoxon
from lists import get_data
from lists import build_exclude_names
import sys

exclude_names = build_exclude_names("all_ibs.tsv", "all_hc.tsv")
ibs_data = get_data("all_ibs.tsv", exclude_names)
hc_data = get_data("all_hc.tsv", exclude_names)


common_names = set(ibs_data.keys()) & set(hc_data.keys())

results = []

threshold = 6.5
for name in sorted(common_names):
    vals_ibs = [x for x in ibs_data[name] if x > threshold]
    vals_hc = [x for x in hc_data[name] if x > threshold]

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


print("локус\tp-value")
for name, pval in results:
    print(f"{name}\t{pval}")

