import rasterio
from rasterio import features
import numpy as np
from skimage.morphology import reconstruction
from scipy import ndimage as ndi
import geopandas as gpd
from shapely.geometry import shape
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# stricter: must be a real enclosed, flat, single-level basin
SLOPE_MAX, SINK_GROW = 0.75, 0.20
RANGE_MAX, DEPTH_MIN, AREA_MIN = 0.50, 1.00, 15

path = r"C:\QGIS\Water_Body\topo.tif"
with rasterio.open(path) as ds:
    dem = ds.read(1).astype("float64")
    res, transform, crs = ds.res[0], ds.transform, ds.crs
nodata = ~np.isfinite(dem); valid = ~nodata
fg = dem.copy(); fg[nodata] = np.nanmean(dem)
gy, gx = np.gradient(fg, res, res)
slope = np.degrees(np.arctan(np.hypot(gx, gy)))
work = dem.copy(); work[nodata] = np.nanmax(dem) + 100
seed = work.copy(); seed[1:-1, 1:-1] = work.max()
sink = reconstruction(seed, work, method="erosion") - work; sink[nodata] = 0.0

cand = ndi.binary_fill_holes(valid & (slope < SLOPE_MAX) & (sink > SINK_GROW))
lbl, n = ndi.label(cand, structure=np.ones((3, 3)))
demf = np.where(valid, dem, np.nan)
keep = np.zeros_like(cand); rows = []
for i in range(1, n + 1):
    m = lbl == i; area = int(m.sum())
    if area < AREA_MIN: continue
    vals = demf[m]; erange = float(np.nanmax(vals) - np.nanmin(vals))
    depth = float(sink[m].max())
    if erange <= RANGE_MAX and depth >= DEPTH_MIN:
        keep |= m
        ys, xs = np.where(m)
        cx, cy = transform * (xs.mean() + 0.5, ys.mean() + 0.5)
        rows.append((i, area, area*res*res/1e4, erange, depth, cx, cy, ys, xs))

rows.sort(key=lambda r: -r[1])
print("high-confidence basins:", len(rows))
print("\n rank  ha    range  depth   easting     northing")
for k, r in enumerate(rows[:20], 1):
    print("%4d %6.2f %6.2f %6.2f  %10.1f %11.1f" % (k, r[2], r[3], r[4], r[5], r[6]))

geoms = [shape(g) for g, v in features.shapes(keep.astype("uint8"), mask=keep, transform=transform) if v == 1]
gdf = gpd.GeoDataFrame({"id": range(len(geoms))}, geometry=geoms, crs=crs)
gdf["area_ha"] = gdf.area / 1e4
gdf.to_file(r"C:\QGIS\Water_Body\water_highconf.gpkg", driver="GPKG")
print("\nwrote water_highconf.gpkg ->", len(gdf), "polygons, total",
      round(gdf.area_ha.sum(), 1), "ha")

# zoom panels on top 4
az, alt = np.radians(315), np.radians(45)
slp = np.arctan(np.hypot(gx, gy)); asp = np.arctan2(-gx, gy)
hs = np.sin(alt)*np.cos(slp) + np.cos(alt)*np.sin(slp)*np.cos(az - asp)
fig, axs = plt.subplots(1, 4, figsize=(20, 6))
for ax, r in zip(axs, rows[:4]):
    ys, xs = r[7], r[8]
    y0, y1 = max(ys.min()-25, 0), min(ys.max()+25, hs.shape[0])
    x0, x1 = max(xs.min()-25, 0), min(xs.max()+25, hs.shape[1])
    ax.imshow(hs[y0:y1, x0:x1], cmap="gray")
    mk = np.zeros(hs.shape, bool); mk[ys, xs] = True
    ax.imshow(np.ma.masked_where(~mk[y0:y1, x0:x1], mk[y0:y1, x0:x1]),
              cmap="cool", alpha=0.55)
    ax.set_title("%.1f ha  depth %.1fm" % (r[2], r[4])); ax.axis("off")
plt.tight_layout()
plt.savefig(r"C:\QGIS\Water_Body\top_basins.png", dpi=95, bbox_inches="tight")
print("wrote top_basins.png")
