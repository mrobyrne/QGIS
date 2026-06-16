"""Sentinel-2 water + marsh detector.
 open water : NDWI > 0
 marsh      : transitional NDWI, high moisture (NDMI), vegetated (not forest),
              on flat terrain  -> catches vegetated wetland RGB/DEM missed."""
import rasterio
from rasterio import features
import numpy as np
from scipy import ndimage as ndi
import geopandas as gpd
from shapely.geometry import shape
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- thresholds ----
OPEN_NDWI   = 0.00
MARSH_NDWI  = -0.30    # wetter than dry upland (NDWI above this)
MARSH_NDMI  = 0.35     # high vegetation/soil moisture
MARSH_NDVI_MAX = 0.80  # exclude dense forest canopy
SLOPE_MAX   = 3.0      # marsh sits on flat ground (deg)
MIN_AREA_M2 = 800

s2 = rasterio.open(r"C:\QGIS\Water_Body\sentinel2_aoi.tif")
G, Rr, NIR, SWIR, ndwi, ndmi = [s2.read(i).astype("float32") for i in range(1, 7)]
ndvi = (NIR - Rr) / (NIR + Rr + 1e-6)
transform, crs = s2.transform, s2.crs
px_m2 = abs(transform.a * transform.e)

with rasterio.open(r"C:\QGIS\Water_Body\topo.tif") as dem_ds:
    dem = dem_ds.read(1).astype("float64")
nod = ~np.isfinite(dem)
fg = dem.copy(); fg[nod] = np.nanmean(dem)
gy, gx = np.gradient(fg, 10.0, 10.0)
slope = np.degrees(np.arctan(np.hypot(gx, gy)))

open_water = ndwi > OPEN_NDWI
marsh = (ndwi > MARSH_NDWI) & (ndwi <= OPEN_NDWI) & (ndmi > MARSH_NDMI) \
        & (ndvi < MARSH_NDVI_MAX) & (slope < SLOPE_MAX)

cls = np.zeros(ndwi.shape, "uint8")
cls[marsh] = 2
cls[open_water] = 1            # open water wins over marsh
cls[nod] = 0

# clean + size filter per class
out = np.zeros_like(cls)
for c in (1, 2):
    m = ndi.binary_opening(cls == c, iterations=1)
    lbl, n = ndi.label(m, structure=np.ones((3, 3)))
    if n:
        sizes = ndi.sum(np.ones_like(lbl), lbl, range(1, n + 1))
        keep = {i + 1 for i, s in enumerate(sizes) if s * px_m2 >= MIN_AREA_M2}
        out[np.isin(lbl, list(keep))] = c

print("open-water ha: %.1f   marsh ha: %.1f"
      % ((out == 1).sum()*px_m2/1e4, (out == 2).sum()*px_m2/1e4))

# vectorize both classes
recs = []
for c, label in [(1, "open_water"), (2, "marsh")]:
    mask = out == c
    for g, v in features.shapes(mask.astype("uint8"), mask=mask, transform=transform):
        if v == 1:
            recs.append({"class": label, "geometry": shape(g)})
gdf = gpd.GeoDataFrame(recs, crs=crs)
gdf["area_ha"] = gdf.area / 1e4
gdf = gdf.sort_values("area_ha", ascending=False).reset_index(drop=True)
gdf.to_file(r"C:\QGIS\Water_Body\water_s2.gpkg", driver="GPKG")
print("wrote water_s2.gpkg ->", len(gdf), "polys (",
      (gdf['class'] == 'open_water').sum(), "open,",
      (gdf['class'] == 'marsh').sum(), "marsh )")

# overview figure
ds = 5
fig, ax = plt.subplots(1, 3, figsize=(22, 8))
ax[0].imshow(ndwi[::ds, ::ds], cmap="BrBG", vmin=-0.6, vmax=0.6)
ax[0].set_title("NDWI");
ax[1].imshow(ndmi[::ds, ::ds], cmap="YlGnBu", vmin=-0.2, vmax=0.6)
ax[1].set_title("NDMI (moisture)")
base = np.dstack([np.clip(Rr, 0, 3000)/3000, np.clip(G, 0, 3000)/3000,
                  np.clip(NIR, 0, 5000)/5000])[::ds, ::ds]   # false-colour
ax[2].imshow(base)
ov = out[::ds, ::ds]
ax[2].imshow(np.ma.masked_where(ov != 1, ov), cmap="cool", alpha=0.9)
ax[2].imshow(np.ma.masked_where(ov != 2, ov), cmap="autumn", alpha=0.9)
ax[2].set_title("open water (cyan) + marsh (red)")
for a in ax: a.axis("off")
plt.tight_layout(); plt.savefig(r"C:\QGIS\Water_Body\s2_detect.png", dpi=95,
                                bbox_inches="tight")
print("wrote s2_detect.png")
