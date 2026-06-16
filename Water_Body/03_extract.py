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

# ---------------- tunables ----------------
SLOPE_MAX   = 0.75   # deg  -- water is flat
SINK_MIN    = 0.20   # m    -- must sit in a depression
RANGE_MAX   = 0.75   # m    -- elevation spread across the whole blob (one water level)
AREA_MIN    = 12     # cells (*100 m2) -> ~1200 m2 minimum
# ------------------------------------------

path = r"C:\QGIS\Water_Body\topo.tif"
with rasterio.open(path) as ds:
    dem = ds.read(1).astype("float64")
    res = ds.res[0]; transform = ds.transform; crs = ds.crs
nodata = ~np.isfinite(dem); valid = ~nodata

# slope
fg = dem.copy(); fg[nodata] = np.nanmean(dem)
gy, gx = np.gradient(fg, res, res)
slope = np.degrees(np.arctan(np.hypot(gx, gy)))

# depression fill -> sink depth (voids act as high barriers)
work = dem.copy(); work[nodata] = np.nanmax(dem) + 100
seed = work.copy(); seed[1:-1, 1:-1] = work.max()
filled = reconstruction(seed, work, method="erosion")
sink_depth = filled - work; sink_depth[nodata] = 0.0

# candidate flat-low cells
cand = valid & (slope < SLOPE_MAX) & (sink_depth > SINK_MIN)
cand = ndi.binary_fill_holes(cand)

# connected components + region-level tests
lbl, n = ndi.label(cand, structure=np.ones((3, 3)))
print("raw components:", n)
keep = np.zeros_like(cand)
rows = []
demf = np.where(valid, dem, np.nan)
for i in range(1, n + 1):
    m = lbl == i
    area = int(m.sum())
    if area < AREA_MIN:
        continue
    vals = demf[m]
    erange = float(np.nanmax(vals) - np.nanmin(vals))
    depth = float(sink_depth[m].max())
    passed = erange <= RANGE_MAX
    rows.append((i, area, area * res * res / 1e4, erange, depth, passed))
    if passed:
        keep |= m

rows.sort(key=lambda r: -r[1])
print("\n  id   cells   ha     elev_range  maxdepth  keep")
for r in rows[:30]:
    print("%5d %6d %6.2f %9.2f %9.2f   %s" % (r[0], r[1], r[2], r[3], r[4], r[5]))
print("\nKEPT regions:", sum(1 for r in rows if r[5]),
      " total ha:", round(sum(r[2] for r in rows if r[5]), 1))

# vectorize kept mask
geoms = []
for geom, val in features.shapes(keep.astype("uint8"), mask=keep, transform=transform):
    if val == 1:
        geoms.append(shape(geom))
gdf = gpd.GeoDataFrame({"id": range(len(geoms))}, geometry=geoms, crs=crs)
gdf["area_ha"] = gdf.area / 1e4
out_gpkg = r"C:\QGIS\Water_Body\water_candidates.gpkg"
gdf.to_file(out_gpkg, driver="GPKG")
print("wrote", out_gpkg, "->", len(gdf), "polygons")

# overlay preview
az, alt = np.radians(315), np.radians(45)
slp = np.arctan(np.hypot(gx, gy)); asp = np.arctan2(-gx, gy)
hs = np.sin(alt)*np.cos(slp) + np.cos(alt)*np.sin(slp)*np.cos(az - asp)
fig, a = plt.subplots(figsize=(11, 12))
a.imshow(hs, cmap="gray")
ov = np.ma.masked_where(~keep, keep)
a.imshow(ov, cmap="cool", alpha=0.95)
a.set_title("Detected water (flat + depression + single-level + area)")
a.axis("off")
plt.savefig(r"C:\QGIS\Water_Body\detected.png", dpi=95, bbox_inches="tight")
print("wrote detected.png")
