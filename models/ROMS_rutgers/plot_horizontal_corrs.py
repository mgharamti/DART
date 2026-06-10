import glob
import cmocean
import numpy             as np
import xarray            as xr
import cartopy.crs       as ccrs
import cartopy.feature   as cfeature
import matplotlib.pyplot as plt


# Configuration
ens_dir   = "/glade/u/home/gharamti/scratch/inacawo/roms_rst/roms_ens80"
var       = "salt" # "temp", "salt"

xi0       = 1020
eta0      = 330
ref_depth = 10.   

files = sorted(glob.glob(f"{ens_dir}/roms_mem_*.nc"))

# Local window around reference point
# Don't compute correlations for the entire domain
dx = 250
dy = 250


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

    return -z

def circle_lonlat(lon0, lat0, radius_km, n=361):
    """
    Great-circle circle around lon0, lat0 with radius_km.
    Returns lon/lat arrays.
    """
    R = 6371.0
    ang = radius_km / R

    lon0r = np.deg2rad(lon0)
    lat0r = np.deg2rad(lat0)
    bearing = np.linspace(0, 2*np.pi, n)

    lat = np.arcsin(
        np.sin(lat0r) * np.cos(ang)
        + np.cos(lat0r) * np.sin(ang) * np.cos(bearing)
    )

    lon = lon0r + np.arctan2(
        np.sin(bearing) * np.sin(ang) * np.cos(lat0r),
        np.cos(ang) - np.sin(lat0r) * np.sin(lat)
    )

    return np.rad2deg(lon), np.rad2deg(lat)



######
# main body here: 
# Use first member to find vertical level nearest ref_depth
with xr.open_dataset(files[0], decode_times=False) as ds0:
    depth0 = roms_z_rho(ds0, xi0, eta0)
    k0 = nearest_k(depth0, ref_depth)

    print(f"\nRequested depth: {ref_depth:.1f} m")
    print(f"Using ROMS level index: {k0}")
    print(f"Actual depth at reference point: {depth0[k0]:.1f} m")

    nx = ds0.sizes["xi_rho"]
    ny = ds0.sizes["eta_rho"]

    x1 = max(0, xi0 - dx)
    x2 = min(nx, xi0 + dx + 1)
    y1 = max(0, eta0 - dy)
    y2 = min(ny, eta0 + dy + 1)

    lon = ds0["lon_rho"].isel(xi_rho=slice(x1, x2), eta_rho=slice(y1, y2)).values
    lat = ds0["lat_rho"].isel(xi_rho=slice(x1, x2), eta_rho=slice(y1, y2)).values


# Read ensemble horizontal slices
slices = []

for f in files:
    with xr.open_dataset(f, decode_times=False) as ds:
        da = ds[var].squeeze()

        field = da.isel(
            s_rho=k0,
            xi_rho=slice(x1, x2),
            eta_rho=slice(y1, y2)
        ).values

        slices.append(field)

E = np.asarray(slices)   # ens_size x eta_window x xi_window
nens = E.shape[0]

print(f"\nRead {nens} members")
print(f"Ensemble array shape: {E.shape}")


# Ensemble anomalies
A = E - E.mean(axis=0, keepdims=True)

# Reference ensemble anomaly at xi0, eta0 inside the window
iref = eta0 - y1
jref = xi0 - x1

x = A[:, iref, jref]     # ens_size
Y = A                    # ens_size x y x x

num = np.sum(x[:, None, None] * Y, axis=0)
den = np.sqrt(np.sum(x**2) * np.sum(Y**2, axis=0))

corr = num / den

# data for the loc circle: 
cutoff1, cutoff2 = 0.02355, 0.03923
earth_radius_km  = 6371.0

zero_weight_km1  = 2.0 * cutoff1 * earth_radius_km  
zero_weight_km2  = 2.0 * cutoff2 * earth_radius_km  

lon0 = lon[iref, jref]
lat0 = lat[iref, jref]

clon1, clat1 = circle_lonlat(lon0, lat0, zero_weight_km1)
clon2, clat2 = circle_lonlat(lon0, lat0, zero_weight_km2)

# Plot
proj = ccrs.PlateCarree()
fig  = plt.figure(figsize=(10,8))
ax   = plt.axes(projection=proj)

# ax.set_facecolor('aliceblue') 
ax.add_feature(cfeature.LAND, facecolor='whitesmoke') 
ax.add_feature(cfeature.NaturalEarthFeature('physical', 'coastline', '10m'), edgecolor='black', facecolor='none')
ax.add_feature(cfeature.LAKES, facecolor='lightsteelblue')
ax.add_feature(cfeature.BORDERS, linewidth=0.5)

pcm  = ax.pcolormesh(lon, lat, corr, 
                     cmap= cmocean.cm.curl, #'RdBu',
                     transform=proj, 
                     vmin=-1, vmax=1)

# plot reference point
ax.plot(lon[iref, jref], lat[iref, jref], "bo", markersize=10, label="reference point")

# draw loc circle
ax.plot(clon1, clat1, "k-", linewidth=2, label=f"{zero_weight_km1:.0f} km 0-weight")
ax.plot(clon2, clat2, "k--", linewidth=2, label=f"{zero_weight_km2:.0f} km 0-weight")

cbar = fig.colorbar(pcm, ax=ax, shrink=0.4, label="ensemble correlations")

ax.set_title(f"{var} horizontal ensemble correlations\n"
             f"Ref: lon={lon0:.1f}, lat={lat0:.1f}, depth≈{depth0[k0]:.1f} m",
             fontsize= 18)

ax.legend(
    loc='upper left', 
    bbox_to_anchor=(1.05, 1), 
    frameon=False)

gl = ax.gridlines(draw_labels=True, linewidth=0.5, color="gray", alpha=0.2)
gl.top_labels   = False
gl.right_labels = False

ax.set_extent([108, 130, -15, 2.5], crs=proj)

plt.tight_layout()

outfile = f"{var}_horizontal_corr_xi{xi0}_eta{eta0}_z{k0}.png"
plt.savefig(outfile, dpi=200)

plt.show()
