"""Final detector:
 open_water = Sentinel-2 NDWI>0 pixels (river + ponds, high recall)
 marsh      = DEM flat-basin, not open, vegetated, high moisture (NDMI)
Marsh threshold calibrated from per-basin stats (NDMI splits wet vs dry veg)."""
import rasterio
from rasterio import features
from rasterio.features import geometry_mask
from rasterio.windows import from_bounds
import numpy as np
from scipy import ndimage as ndi
import geopandas as gpd
from shapely.geometry import shape, mapping
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OPEN_NDWI   = 0.0
MARSH_NDMI  = 0.18
MARSH_NDVI  = 0.40
OPEN_FRAC   = 0.30
MIN_AREA_M2 = 800

s2 = rasterio.open(r"C:\QGIS\Water_Body\sentinel2_aoi.tif")
G, Rr, NIR, SWIR, ndwi, ndmi = [s2.read(i).astype("float32") for i in range(1, 7)]
ndvi = (NIR - Rr) / (NIR + Rr + 1e-6)
transform, crs = s2.transform, s2.crs
px_m2 = abs(transform.a * transform.e)

# ---- open water: NDWI pixels ----
ow = ndwi > OPEN_NDWI
ow = ndi.binary_opening(ow, iterations=1)
lbl, n = ndi.label(ow, structure=np.ones((3, 3)))
sizes = ndi.sum(np.ones_like(lbl), lbl, range(1, n + 1))
ow = np.isin(lbl, [i+1 for i, s in enumerate(sizes) if s*px_m2 >= MIN_AREA_M2])

# ---- marsh: classify DEM basins ----
basins = gpd.read_file(r"C:\QGIS\Water_Body\water_candidates.gpkg")
basins = basins[basins.area > 2000].reset_index(drop=True)
marsh_geoms = []
for geom in basins.geometry:
    minx, miny, maxx, maxy = geom.bounds
    win = from_bounds(minx, miny, maxx, maxy, s2.transform)
    arr = s2.read(window=win)
    if arr.shape[1] == 0 or arr.shape[2] == 0:
        continue
    wt = s2.window_transform(win)
    m = geometry_mask([geom], (arr.shape[1], arr.shape[2]), wt, invert=True)
    if m.sum() < 4:
        continue
    nw, nm = arr[4][m], arr[5][m]
    nv = (arr[2][m]-arr[1][m])/(arr[2][m]+arr[1][m]+1e-6)
    open_frac = (nw > OPEN_NDWI).mean()
    if open_frac > OPEN_FRAC:
        continue                                   # already open water
    if np.median(nm) > MARSH_NDMI and np.median(nv) > MARSH_NDVI:
        marsh_geoms.append(geom)

print("open-water ha: %.1f" % (ow.sum()*px_m2/1e4))
print("marsh basins :", len(marsh_geoms),
      " ha %.1f" % (sum(g.area for g in marsh_geoms)/1e4))

recs = []
for g, v in features.shapes(ow.astype("uint8"), mask=ow, transform=transform):
    if v == 1:
        recs.append({"class": "open_water", "geometry": shape(g)})
for g in marsh_geoms:
    recs.append({"class": "marsh", "geometry": g})
gdf = gpd.GeoDataFrame(recs, crs=crs)
gdf["area_ha"] = gdf.area / 1e4
gdf.to_file(r"C:\QGIS\Water_Body\water_s2_v2.gpkg", driver="GPKG")
print("wrote water_s2_v2.gpkg ->", len(gdf), "polys")

# verify marsh candidates on Esri imagery
if marsh_geoms:
    img = rasterio.open(r"C:\QGIS\Water_Body\esri_rgb_reproj.tif")
    mg = sorted(marsh_geoms, key=lambda g: -g.area)
    nshow = min(8, len(mg)); cols = 4; rows = (nshow+cols-1)//cols
    fig, axs = plt.subplots(rows, cols, figsize=(20, 5*rows))
    for k, ax in enumerate(np.array(axs).ravel()):
        if k >= nshow:
            ax.axis("off"); continue
        geom = mg[k]; minx, miny, maxx, maxy = geom.bounds; pad = 60
        win = from_bounds(minx-pad, miny-pad, maxx+pad, maxy+pad, img.transform)
        rgb = img.read([1, 2, 3], window=win); wt = img.window_transform(win)
        ax.imshow(np.transpose(rgb, (1, 2, 0)),
                  extent=[wt.c, wt.c+rgb.shape[2]*wt.a, wt.f+rgb.shape[1]*wt.e, wt.f])
        xs, ys = geom.exterior.xy
        ax.plot(xs, ys, color="yellow", lw=2)
        ax.set_title("marsh %.2f ha" % (geom.area/1e4)); ax.set_xticks([]); ax.set_yticks([])
    plt.tight_layout(); plt.savefig(r"C:\QGIS\Water_Body\verify_marsh.png", dpi=90,
                                    bbox_inches="tight")
    print("wrote verify_marsh.png")
