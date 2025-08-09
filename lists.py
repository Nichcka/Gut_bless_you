import math

def get_names_from_file(filename):
    names = set()
    with open(filename, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue
            names.add(parts[1])
    return names

def build_exclude_names(file1, file2):
    names1 = get_names_from_file(file1)
    names2 = get_names_from_file(file2)
    return names1 ^ names2  

def get_data(filename, exclude_names):
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
