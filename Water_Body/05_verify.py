import rasterio
from rasterio.windows import from_bounds
import numpy as np
import geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

img_path = r"C:\QGIS\Water_Body\esri_rgb_reproj.tif"
gdf = gpd.read_file(r"C:\QGIS\Water_Body\water_highconf.gpkg").sort_values(
    "area_ha", ascending=False).reset_index(drop=True)

src = rasterio.open(img_path)
n = min(8, len(gdf))
fig, axs = plt.subplots(2, 4, figsize=(20, 10))
for k, ax in enumerate(axs.ravel()):
    if k >= n:
        ax.axis("off"); continue
    geom = gdf.geometry[k]
    minx, miny, maxx, maxy = geom.bounds
    pad = 60
    win = from_bounds(minx-pad, miny-pad, maxx+pad, maxy+pad, src.transform)
    rgb = src.read([1, 2, 3], window=win)
    wt = src.window_transform(win)
    ax.imshow(np.transpose(rgb, (1, 2, 0)),
              extent=[wt.c, wt.c + rgb.shape[2]*wt.a,
                      wt.f + rgb.shape[1]*wt.e, wt.f])
    xs, ys = geom.exterior.xy
    ax.plot(xs, ys, color="cyan", lw=2)
    ax.set_title("rank %d  %.1f ha" % (k+1, gdf.area_ha[k]))
    ax.set_xticks([]); ax.set_yticks([])
plt.tight_layout()
plt.savefig(r"C:\QGIS\Water_Body\verify.png", dpi=95, bbox_inches="tight")
print("wrote verify.png")
