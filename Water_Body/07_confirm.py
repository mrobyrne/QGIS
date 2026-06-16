"""Combined water detector: imagery colour test  AND  DEM flat-basin mask.
DEM basins remove water-coloured shadows/roofs outside depressions;
imagery removes dry/graded/vegetated basins. Intersection = high precision."""
import rasterio
from rasterio import features
import numpy as np
from scipy import ndimage as ndi
import geopandas as gpd
from shapely.geometry import shape
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# water colour rule (0-255), calibrated on lake vs developed/grass
BR_MIN, EXG_MAX, BRIGHT_MAX = 12, 45, 110
MIN_AREA_M2 = 800            # drop specks

img_path = r"C:\QGIS\Water_Body\esri_rgb_reproj.tif"
basins = gpd.read_file(r"C:\QGIS\Water_Body\water_candidates.gpkg")  # broad DEM set

with rasterio.open(img_path) as src:
    rgb = src.read([1, 2, 3]).astype(float)
    transform, crs, shape_hw = src.transform, src.crs, (src.height, src.width)
    px_m2 = abs(src.transform.a * src.transform.e)

R, G, B = rgb
bright = (R + G + B) / 3
exg = 2*G - R - B
water_col = (B - R >= BR_MIN) & (exg <= EXG_MAX) & (bright <= BRIGHT_MAX)

# rasterize DEM basins onto the imagery grid
basin_mask = features.rasterize(
    ((g, 1) for g in basins.geometry), out_shape=shape_hw,
    transform=transform, fill=0, dtype="uint8").astype(bool)

combined = water_col & basin_mask
combined = ndi.binary_opening(combined, iterations=1)
combined = ndi.binary_fill_holes(combined)

# size filter
lbl, n = ndi.label(combined, structure=np.ones((3, 3)))
sizes = ndi.sum(np.ones_like(lbl), lbl, range(1, n + 1))
keep_ids = {i + 1 for i, s in enumerate(sizes) if s * px_m2 >= MIN_AREA_M2}
final = np.isin(lbl, list(keep_ids))

print("water-colour px:", int(water_col.sum()))
print("inside basins  :", int(combined.sum()))
print("after size flt  :", int(final.sum()),
      " -> %.1f ha" % (final.sum() * px_m2 / 1e4))

geoms = [shape(g) for g, v in features.shapes(final.astype("uint8"),
         mask=final, transform=transform) if v == 1]
gdf = gpd.GeoDataFrame({"id": range(len(geoms))}, geometry=geoms, crs=crs)
gdf["area_ha"] = gdf.area / 1e4
gdf = gdf.sort_values("area_ha", ascending=False).reset_index(drop=True)
gdf.to_file(r"C:\QGIS\Water_Body\water_final.gpkg", driver="GPKG")
print("wrote water_final.gpkg ->", len(gdf), "polygons")
print(gdf.head(12)[["id", "area_ha"]].to_string(index=False))

# preview: imagery downsampled + final water overlaid
ds = 6
small = np.transpose(rgb[:, ::ds, ::ds], (1, 2, 0)).astype("uint8")
fig, ax = plt.subplots(figsize=(12, 13))
ax.imshow(small)
ov = final[::ds, ::ds]
ax.imshow(np.ma.masked_where(~ov, ov), cmap="autumn", alpha=0.8)
ax.set_title("Final water = imagery colour AND DEM basin"); ax.axis("off")
plt.savefig(r"C:\QGIS\Water_Body\final.png", dpi=95, bbox_inches="tight")
print("wrote final.png")
