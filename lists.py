import math

exclude_names = {
    "MGYG000001539_14:1200701:1214584_nointergenome_regions",
    "MGYG000003682_6:80423:96047_nointergenome_regions",
    "MGYG000006849.fa_35:19289:32391_nointergenome_regions",
    "MGYG000086822_124:24542:36192_nointergenome_regions",
    "MGYG000117568_24:25311:35480_nointergenome_regions",
    "MGYG000182605_28:523:10749_nointergenome_regions",
    "MGYG000209896_54:0:10859_nointergenome_regions",
    "MGYG000217771_68:1409:7804_nointergenome_regions"
}

def get_data(filename):
    data = {}
    with open(filename, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            name = parts[1]
            if name in exclude_names:
                continue
            try:
                value = float(parts[2])
            except ValueError:
                continue
            if value <= 0:
                continue 
            log_value = math.log(value)
            if name not in data:
                data[name] = []
            data[name].append(log_value)
    return data

