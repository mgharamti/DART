import glob
import cmocean
import numpy             as np
import xarray            as xr
import matplotlib.pyplot as plt


# Configuration
ens_dir   = "/glade/u/home/gharamti/scratch/inacawo/roms_rst/roms_ens80"
var       = "salt"  # "temp", "salt"

xi0       = 1020
eta0      = 330
ref_depth = 280.0    # meters, positive down

files = sorted(glob.glob(f"{ens_dir}/roms_mem_*.nc"))


# Helper Functions
def nearest_k(depth, target):
    return int(np.argmin(np.abs(depth - target)))


def roms_z_rho(ds, xi, eta):
    h = ds["h"].isel(xi_rho=xi, eta_rho=eta).values
    zeta = ds["zeta"].isel(xi_rho=xi, eta_rho=eta).squeeze().values

    s_rho = ds["s_rho"].values
    Cs_r  = ds["Cs_r"].values
    hc    = float(ds["hc"].values)

    z0 = (hc * s_rho + h * Cs_r) / (hc + h)
    z  = zeta + (zeta + h) * z0

    return -z   # positive down


def gaspari_cohn(r):
    """
    Gaspari-Cohn localization, where r = distance / half-width.
    Zero at r >= 2.
    """
    r = np.asarray(r)
    gc = np.zeros_like(r, dtype=float)

    m1 = (r >= 0) & (r <= 1)
    m2 = (r >  1) & (r <  2)

    gc[m1] = (
        1
        - 5/3*r[m1]**2
        + 5/8*r[m1]**3
        + 0.5*r[m1]**4
        - 1/4*r[m1]**5
    )

    gc[m2] = (
        4
        - 5*r[m2]
        + 5/3*r[m2]**2
        + 5/8*r[m2]**3
        - 0.5*r[m2]**4
        + 1/12*r[m2]**5
        - 2/(3*r[m2])
    )

    return gc


######
# main body here:
# Use first member to find profile depth and reference level
with xr.open_dataset(files[0], decode_times=False) as ds0:
    depth0 = roms_z_rho(ds0, xi0, eta0)
    k0 = nearest_k(depth0, ref_depth)

    lon0 = ds0["lon_rho"].isel(xi_rho=xi0, eta_rho=eta0).values
    lat0 = ds0["lat_rho"].isel(xi_rho=xi0, eta_rho=eta0).values

    print(f"\nRequested depth: {ref_depth:.1f} m")
    print(f"Using ROMS level index: {k0}")
    print(f"Actual depth at reference point: {depth0[k0]:.1f} m")
    print(f"Reference lon/lat: {lon0:.2f}, {lat0:.2f}")


# Read ensemble vertical profiles
profiles = []

for f in files:
    with xr.open_dataset(f, decode_times=False) as ds:
        da = ds[var].squeeze()

        prof = da.isel(
            xi_rho=xi0,
            eta_rho=eta0
        ).values

        profiles.append(prof)

E = np.asarray(profiles)   # ens_size x nz
nens = E.shape[0]

print(f"\nRead {nens} members")
print(f"Ensemble array shape: {E.shape}")


# Ensemble anomalies
A = E - E.mean(axis=0, keepdims=True)

# Reference ensemble anomaly at chosen depth
x = A[:, k0]     # ens_size
Y = A            # ens_size x nz

num = np.sum(x[:, None] * Y, axis=0)
den = np.sqrt(np.sum(x**2) * np.sum(Y**2, axis=0))

corr = num / den


# Optional: plot DART vertical localization curves for comparison
# These are 0-weight distances in meters
vert_zero_weight = [100.0, 250.0, 500.0]

loc_curves = {}
dz = np.abs(depth0 - depth0[k0])

for zw in vert_zero_weight:
    half_width = zw / 2.0
    r = dz / half_width
    loc_curves[zw] = gaspari_cohn(r)

localized_corrs = {}

for zw in vert_zero_weight:
    localized_corrs[zw] = corr * loc_curves[zw]

# Plot
fig, ax = plt.subplots(figsize=(7, 9))

ax.plot(corr, depth0, "k-", linewidth=2.5, label="ensemble corr")

for zw in vert_zero_weight:
    ax.plot(
        localized_corrs[zw], depth0,
        linewidth=2,
        linestyle="--",
        label=f"corr × GC, {zw:.0f} m 0-weight"
    )

ax.axvline(0.0, color="0.5", linewidth=1)
ax.axvline(0.2, color="0.7", linewidth=1, linestyle=":")
ax.axvline(-0.2, color="0.7", linewidth=1, linestyle=":")

ax.plot(1.0, depth0[k0], "bo", markersize=8, label="reference depth")

ax.invert_yaxis()
ax.grid(True, linestyle="--", linewidth=0.5, alpha=0.5)

ax.set_xlabel("Correlation / localization weight", fontsize=14)
ax.set_ylabel("Depth (m)", fontsize=14)

ax.set_title(
    f"{var} vertical raw and localized correlations\n"
    f"Ref: lon={lon0:.2f}, lat={lat0:.2f}, depth≈{depth0[k0]:.1f} m",
    fontsize=18
)

ax.legend(loc="best", frameon=True)

plt.tight_layout()

outfile = f"{var}_vertical_localized_corr_xi{xi0}_eta{eta0}_z{k0}.png"
plt.savefig(outfile, dpi=200)

plt.show()
