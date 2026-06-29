from icecube import dataio, icetray, dataclasses
import numpy as np
from collections import defaultdict

# -----------------------------
# USER INPUT
# -----------------------------
filename = "/project/6008051/pone_simulation/geometry_subselections/70String_default_geometry/000002/results/linefit_2304.i3.zst"
print("working")
# -----------------------------
# STORAGE
# -----------------------------
om_counts = defaultdict(int)
om_times = defaultdict(list)

total_pulses = 0

# -----------------------------
# OPEN FILE
# -----------------------------
f = dataio.I3File(filename)

while f.more():
    frame = f.pop_frame()

    # skip non-physics frames
    if "Noise_K40" not in frame:
        continue

    k40_map = frame["Noise_K40"]

    # k40_map is an I3RecoPulseSeriesMap: OMKey -> vector<I3RecoPulse>
    for omkey, pulses in k40_map.items():

        n_pulses = len(pulses)
        om_counts[omkey] += n_pulses
        total_pulses += n_pulses

        for p in pulses:
            om_times[omkey].append(p.time)

# -----------------------------
# BASIC SUMMARY
# -----------------------------
print("Total K40 pulses:", total_pulses)
print("Number of PMTs:", len(om_counts))

# -----------------------------
# COMPUTE SIMPLE RATES
# -----------------------------
# You need the livetime; if not stored, estimate from pulse span
all_times = []

for tlist in om_times.values():
    all_times.extend(tlist)

t_min = min(all_times)
t_max = max(all_times)
livetime = t_max - t_min

print(f"Estimated livetime: {livetime:.3f} s")

om_rates = {}

for omkey, count in om_counts.items():
    om_rates[omkey] = count / livetime

# -----------------------------
# PRINT TOP OMs
# -----------------------------
print("\nTop 10 PMTs by rate:")
for omkey, count in sorted(om_counts.items(), key=lambda x: x[1], reverse=True)[:10]:
    rate = om_rates.get(omkey, 0.0)
    print(f"{omkey} | {count} hits | {rate:.3e} Hz")
