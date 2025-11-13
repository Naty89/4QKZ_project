# parse_fpocket_runs.py
import os, glob, csv
fields = ["snapshot","Volume","Hydrophobicity","Polarity","DruggabilityScore","NumberAlphaSpheres"]
with open("fpocket_timeseries.csv","w",newline="") as out:
    writer = csv.DictWriter(out, fieldnames=fields)
    writer.writeheader()
    for d in sorted(glob.glob("fpocket_runs/*")):
        info = os.path.join(d,"out","*_info.txt")
        # fpocket layout differs by version; search for keys in all info files:
        info_files = glob.glob(os.path.join(d,"*info.txt")) + glob.glob(os.path.join(d,"out/*info.txt"))
        if not info_files:
            # try standard fpocket output folder name:
            info_files = glob.glob(os.path.join(d,"*_out/*_info.txt"))
        if not info_files:
            continue
        f = info_files[0]
        data = {"snapshot": os.path.basename(d)}
        with open(f) as fh:
            txt = fh.read()
        # crude extraction - robust enough for fpocket info format you showed earlier
        import re
        def find_float(key):
            m = re.search(rf"{key}\s*:\s*([-\d\.]+)", txt)
            return float(m.group(1)) if m else None
        def find_int(key):
            m = re.search(rf"{key}\s*:\s*([-\d]+)", txt)
            return int(m.group(1)) if m else None
        data["Volume"] = find_float("Volume")
        data["Hydrophobicity"] = find_float("Hydrophobicity score") or find_float("Hydrophobicity")
        data["Polarity"] = find_int("Polarity score") or find_int("Polarity")
        data["DruggabilityScore"] = find_float("Druggability Score") or find_float("Druggability")
        data["NumberAlphaSpheres"] = find_int("Number of Alpha Spheres")
        writer.writerow(data)
print("Wrote fpocket_timeseries.csv")
