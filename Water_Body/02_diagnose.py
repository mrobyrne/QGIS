import rasterio
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from skimage.morphology import reconstruction

path = r"C:\QGIS\Water_Body\topo.tif"
with rasterio.open(path) as ds:
    dem = ds.read(1).astype("float64")
    res = ds.res[0]

nodata = ~np.isfinite(dem)
valid = ~nodata

# ---- slope (degrees), handling nan by filling with neighbour-ish mean first ----
filled_for_grad = dem.copy()
filled_for_grad[nodata] = np.nanmean(dem)
gy, gx = np.gradient(filled_for_grad, res, res)
slope = np.degrees(np.arctan(np.hypot(gx, gy)))
slope[nodata] = np.nan

# ---- depression fill (sink depth) via morphological reconstruction ----
# replace nodata with a HIGH value so voids act as barriers (handled separately)
work = dem.copy()
hi = np.nanmax(dem) + 100
work[nodata] = hi
seed = work.copy()
seed[1:-1, 1:-1] = work.max()
filled = reconstruction(seed, work, method="erosion")
sink_depth = filled - work
sink_depth[nodata] = 0.0

# ---- flatness mask ----
flat = (slope < 1.0) & valid
flat_in_sink = flat & (sink_depth > 0.05)

print("=== diagnostics ===")
print("slope deg  : min %.3f  median %.3f  p95 %.3f  max %.3f" % (
    np.nanmin(slope), np.nanmedian(slope), np.nanpercentile(slope, 95), np.nanmax(slope)))
print("cells slope<1deg      :", int(flat.sum()), "(%.2f%%)" % (100*flat.sum()/valid.sum()))
print("cells flat & in-sink  :", int(flat_in_sink.sum()))
print("sink_depth max        : %.3f m" % sink_depth.max())
print("nodata cells          :", int(nodata.sum()))

# ---- figure ----
fig, ax = plt.subplots(2, 2, figsize=(16, 17))
# hillshade
az, alt = np.radians(315), np.radians(45)
gy2, gx2 = np.gradient(filled_for_grad, res, res)
slp = np.arctan(np.hypot(gx2, gy2)); asp = np.arctan2(-gx2, gy2)
hs = np.sin(alt)*np.cos(slp) + np.cos(alt)*np.sin(slp)*np.cos(az - asp)
ax[0,0].imshow(hs, cmap="gray"); ax[0,0].set_title("Hillshade")
# elevation
el = np.where(valid, dem, np.nan)
im = ax[0,1].imshow(el, cmap="terrain"); ax[0,1].set_title("Elevation")
plt.colorbar(im, ax=ax[0,1], shrink=0.6)
# nodata voids
ax[1,0].imshow(hs, cmap="gray")
nd = np.ma.masked_where(~nodata, nodata)
ax[1,0].imshow(nd, cmap="autumn", alpha=0.9); ax[1,0].set_title("NoData voids (red)")
# candidate water: flat + in-sink
ax[1,1].imshow(hs, cmap="gray")
cw = np.ma.masked_where(~flat_in_sink, flat_in_sink)
ax[1,1].imshow(cw, cmap="cool", alpha=0.9)
ax[1,1].set_title("Candidate water: flat(<1deg) & in depression")
for a in ax.ravel(): a.axis("off")
plt.tight_layout()
out = r"C:\QGIS\Water_Body\diag.png"
plt.savefig(out, dpi=90, bbox_inches="tight")
print("wrote", out)
