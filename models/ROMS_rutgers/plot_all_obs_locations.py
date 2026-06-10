from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import pydartdiags.obs_sequence.obs_sequence as obsq

obs_seq_file = Path("work/obs_seq_all.out")
fig_file = Path("all_obs_locations.png")

obs = obsq.ObsSequence(obs_seq_file)
df = obs.df.copy()

# Group colocated DART observation types by observing platform/product
groups = {
    "Satellite SST": ["SATELLITE_BLENDED_SST"],
    "Satellite SSH": ["SATELLITE_SSH"],
    "SVP Drifters": [
        "DRIFTER_TEMPERATURE",
        "DRIFTER_U_CURRENT_COMPONENT",
        "DRIFTER_V_CURRENT_COMPONENT",
    ],
    "Profiling Floats": [
        "FLOAT_TEMPERATURE",
        "FLOAT_SALINITY",
    ],
    "HF Radar Totals": [
        "HFRADAR_U_CURRENT_COMPONENT",
        "HFRADAR_V_CURRENT_COMPONENT",
    ],
    "HF Radar Radials": [
        "HFRADAR_RADIAL_VELOCITY",
    ],
}

colors = {
    "Satellite SST"   : "tab:green",
    "Satellite SSH"   : "tab:blue",
    "SVP Drifters"    : "tab:orange",
    "Profiling Floats": "tab:purple",
    "HF Radar Totals" : "tab:cyan",
    "HF Radar Radials": "tab:red",
}

sizes = {
    "Satellite SST"   : 5,
    "Satellite SSH"   : 14,
    "SVP Drifters"    : 65,
    "Profiling Floats": 60,
    "HF Radar Totals" : 25,
    "HF Radar Radials": 45,
}

alphas = {
    "Satellite SST"   : 0.25,
    "Satellite SSH"   : 0.85,
    "SVP Drifters"    : 0.95,
    "Profiling Floats": 0.95,
    "HF Radar Totals" : 0.85,
    "HF Radar Radials": 0.65,
}

zorders = {
    "Satellite SST"   : 1,
    "Satellite SSH"   : 4,
    "SVP Drifters"    : 6,
    "Profiling Floats": 7,
    "HF Radar Totals" : 9,
    "HF Radar Radials": 5,
}

proj = ccrs.PlateCarree()

fig = plt.figure(figsize=(12, 7))
ax = plt.axes(projection=proj)

ax.set_facecolor('aliceblue') 
ax.add_feature(cfeature.LAND, facecolor='whitesmoke', zorder=2) 
ax.add_feature(cfeature.NaturalEarthFeature('physical', 'coastline', '10m'), 
               edgecolor='black', facecolor='none', zorder=8)
ax.add_feature(cfeature.LAKES, facecolor='lightsteelblue', zorder=8)
ax.add_feature(cfeature.BORDERS, linewidth=0.5, zorder=8)

for group_name, obs_types in groups.items():

    this_ob = df[df["type"].isin(obs_types)]

    if len(this_ob) == 0:
        continue

    # Avoid over-counting colocated variables in the visual label
    unique_locations = this_ob[["longitude", "latitude"]].drop_duplicates()

    ax.scatter(
        this_ob["longitude"],
        this_ob["latitude"],
        s=sizes[group_name],
        alpha=alphas[group_name],
        color=colors[group_name],
        label=f"{group_name} (N={len(this_ob):,})",
        transform=proj,
        zorder=zorders[group_name],
        edgecolor="none",
    )

ax.set_extent([95, 145, -15, 15], crs=proj)

gl = ax.gridlines(draw_labels=True, linewidth=0.5, color="gray", alpha=0.2)
gl.top_labels = False
gl.right_labels = False

ax.legend(
    loc="center left",
    bbox_to_anchor=(1.02, 0.5),
    fontsize=9,
    frameon=False,
    markerscale=1.5,
)

ax.set_title(
    f"Observations Used for ROMS-DART\n"
    f"N = {len(df):,} observations, {df['type'].nunique()} observation types",
    fontsize=16,
)

fig.savefig(fig_file, dpi=300, bbox_inches="tight")
plt.show()
