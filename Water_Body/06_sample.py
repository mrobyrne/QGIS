import rasterio
from rasterio.windows import from_bounds
from rasterio.features import geometry_mask
import numpy as np
import geopandas as gpd

src = rasterio.open(r"C:\QGIS\Water_Body\esri_rgb_reproj.tif")
gdf = gpd.read_file(r"C:\QGIS\Water_Body\water_highconf.gpkg").sort_values(
    "area_ha", ascending=False).reset_index(drop=True)

def stats(rank):
    geom = gdf.geometry[rank]
    minx, miny, maxx, maxy = geom.bounds
    win = from_bounds(minx, miny, maxx, maxy, src.transform)
    rgb = src.read([1, 2, 3], window=win).astype(float)
    wt = src.window_transform(win)
    m = geometry_mask([geom], (rgb.shape[1], rgb.shape[2]), wt, invert=True)
    R, G, B = rgb[0][m], rgb[1][m], rgb[2][m]
    bright = (R + G + B) / 3
    exg = 2*G - R - B                      # excess green (vegetation)
    bmr = B - R                            # blue-minus-red (water tends >= 0)
    print("rank %d  %.1f ha  n=%d" % (rank+1, gdf.area_ha[rank], m.sum()))
    for name, v in [("bright", bright), ("exg", exg), ("B-R", bmr)]:
        print("   %-7s p10 %6.1f  median %6.1f  p90 %6.1f"
              % (name, np.percentile(v, 10), np.median(v), np.percentile(v, 90)))

for rk in [6, 1, 0, 7, 4]:   # lake, developed, grassy, pond, mixed
    stats(rk)
