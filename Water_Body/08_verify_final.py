import rasterio
from rasterio.windows import from_bounds
import numpy as np
import geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

src = rasterio.open(r"C:\QGIS\Water_Body\esri_rgb_reproj.tif")
gdf = gpd.read_file(r"C:\QGIS\Water_Body\water_final.gpkg").sort_values(
    "area_ha", ascending=False).reset_index(drop=True)
n = len(gdf)
cols = 4; rows = (n + cols - 1) // cols
fig, axs = plt.subplots(rows, cols, figsize=(20, 5*rows))
for k, ax in enumerate(np.array(axs).ravel()):
    if k >= n:
        ax.axis("off"); continue
    geom = gdf.geometry[k]
    minx, miny, maxx, maxy = geom.bounds
    pad = 50
    win = from_bounds(minx-pad, miny-pad, maxx+pad, maxy+pad, src.transform)
    rgb = src.read([1, 2, 3], window=win)
    wt = src.window_transform(win)
    ax.imshow(np.transpose(rgb, (1, 2, 0)),
              extent=[wt.c, wt.c + rgb.shape[2]*wt.a,
                      wt.f + rgb.shape[1]*wt.e, wt.f])
    for poly in (geom.geoms if geom.geom_type == "MultiPolygon" else [geom]):
        xs, ys = poly.exterior.xy
        ax.plot(xs, ys, color="red", lw=1.8)
    ax.set_title("%.2f ha" % gdf.area_ha[k]); ax.set_xticks([]); ax.set_yticks([])
plt.tight_layout()
plt.savefig(r"C:\QGIS\Water_Body\verify_final.png", dpi=90, bbox_inches="tight")
print("wrote verify_final.png,", n, "polygons")
