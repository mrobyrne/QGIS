"""Characterise each DEM flat-basin by its Sentinel-2 signature so we can
separate open-water / marsh / dry from real numbers instead of guesses."""
import rasterio
from rasterio.features import geometry_mask
from rasterio.windows import from_bounds
import numpy as np
import geopandas as gpd

s2 = rasterio.open(r"C:\QGIS\Water_Body\sentinel2_aoi.tif")
basins = gpd.read_file(r"C:\QGIS\Water_Body\water_candidates.gpkg")
basins = basins[basins.area > 4000].reset_index(drop=True)   # >0.4 ha

rows = []
for k, geom in enumerate(basins.geometry):
    minx, miny, maxx, maxy = geom.bounds
    win = from_bounds(minx, miny, maxx, maxy, s2.transform)
    arr = s2.read(window=win)            # 6 bands
    if arr.shape[1] == 0 or arr.shape[2] == 0:
        continue
    wt = s2.window_transform(win)
    m = geometry_mask([geom], (arr.shape[1], arr.shape[2]), wt, invert=True)
    if m.sum() < 4:
        continue
    G, Rr, NIR, SWIR, ndwi, ndmi = [arr[i][m] for i in range(6)]
    ndvi = (NIR - Rr) / (NIR + Rr + 1e-6)
    rows.append((basins.area[k]/1e4, np.median(ndwi), (ndwi > 0).mean(),
                 np.median(ndmi), np.median(ndvi)))

import numpy as np
a = np.array(rows)
print("basins analysed:", len(a))
print("\n   ha   ndwi_med  open_frac  ndmi_med  ndvi_med")
order = np.argsort(-a[:, 0])
for r in a[order][:35]:
    print("%6.2f %8.2f %9.2f %9.2f %8.2f" % (r[0], r[1], r[2], r[3], r[4]))

# crude buckets to see structure
openw = a[a[:, 2] > 0.5]
marshy = a[(a[:, 2] <= 0.5) & (a[:, 3] > 0.2) & (a[:, 1] > -0.35)]
dry = a[(a[:, 2] <= 0.5) & ((a[:, 3] <= 0.2) | (a[:, 1] <= -0.35))]
print("\nopen-water basins (open_frac>0.5):", len(openw))
print("marshy basins (wet, veg)          :", len(marshy),
      " ndmi range %.2f-%.2f" % (marshy[:,3].min(), marshy[:,3].max()) if len(marshy) else "")
print("dry/false basins                  :", len(dry))
