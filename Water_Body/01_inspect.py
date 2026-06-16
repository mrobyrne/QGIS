import rasterio
import numpy as np

path = r"C:\QGIS\Water_Body\topo.tif"
with rasterio.open(path) as ds:
    print("driver   :", ds.driver)
    print("size     :", ds.width, "x", ds.height)
    print("bands    :", ds.count, "dtypes", ds.dtypes)
    print("crs      :", ds.crs)
    print("res      :", ds.res, "(units of CRS)")
    print("bounds   :", ds.bounds)
    print("nodata   :", ds.nodata)
    arr = ds.read(1, masked=True)

data = arr.compressed()
print("--- elevation stats (valid cells) ---")
print("valid cells :", data.size, "of", arr.size)
print("nodata cells:", int(np.ma.count_masked(arr)))
print("min/max     : %.3f / %.3f" % (data.min(), data.max()))
print("mean/median : %.3f / %.3f" % (data.mean(), np.median(data)))
for p in (0.1, 1, 5, 50, 95, 99, 99.9):
    print("  p%-5s = %.3f" % (p, np.percentile(data, p)))
