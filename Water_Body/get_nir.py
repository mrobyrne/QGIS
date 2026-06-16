"""
get_nir.py -- fetch free Sentinel-2 NIR/SWIR (+G,R) from AWS Earth Search
(public COGs, no auth) for the DEM's AOI, warp onto the exact DEM grid,
and save locally as a 4-band GeoTIFF (G, R, NIR, SWIR) plus NDWI/NDMI.

    python get_nir.py --dem topo.tif --out sentinel2_aoi.tif
"""
import argparse
import requests
import numpy as np
import rasterio
from rasterio.vrt import WarpedVRT
from rasterio.enums import Resampling
from rasterio.warp import transform_bounds

STAC = "https://earth-search.aws.element84.com/v1/search"
# Sentinel-2 L2A common-name assets on Earth Search
ASSETS = {"green": "green", "red": "red", "nir": "nir", "swir": "swir16"}
GDAL_ENV = dict(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
                CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif",
                GDAL_HTTP_MULTIRANGE="YES", GDAL_HTTP_MERGE_CONSECUTIVE_RANGES="YES")


def search(bbox, datetime, max_cloud, limit=30):
    body = {
        "collections": ["sentinel-2-l2a"],
        "bbox": list(bbox),
        "datetime": datetime,
        "query": {"eo:cloud_cover": {"lt": max_cloud}},
        "limit": limit,
        "sortby": [{"field": "properties.eo:cloud_cover", "direction": "asc"}],
    }
    r = requests.post(STAC, json=body, timeout=60)
    r.raise_for_status()
    return r.json()["features"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dem", required=True)
    ap.add_argument("--out", default="sentinel2_aoi.tif")
    ap.add_argument("--datetime",
                    default="2024-06-01T00:00:00Z/2025-09-30T23:59:59Z")
    ap.add_argument("--max-cloud", type=float, default=15)
    a = ap.parse_args()

    with rasterio.open(a.dem) as ds:
        dem_crs, dem_t = ds.crs, ds.transform
        dem_w, dem_h, dem_bounds = ds.width, ds.height, ds.bounds
    bbox4326 = transform_bounds(dem_crs, "EPSG:4326", *dem_bounds)
    print("AOI bbox (4326):", [round(x, 4) for x in bbox4326])

    items = search(bbox4326, a.datetime, a.max_cloud)
    print(f"{len(items)} scenes < {a.max_cloud}% cloud")
    if not items:
        raise SystemExit("no scenes found; widen --datetime or --max-cloud")
    it = items[0]
    print("chosen:", it["id"], "cloud %.1f%%" % it["properties"]["eo:cloud_cover"],
          it["properties"]["datetime"][:10])

    bands = {}
    with rasterio.Env(**GDAL_ENV):
        for name, key in ASSETS.items():
            href = it["assets"][key]["href"]
            with rasterio.open(href) as src:
                with WarpedVRT(src, crs=dem_crs, transform=dem_t,
                               width=dem_w, height=dem_h,
                               resampling=Resampling.bilinear) as vrt:
                    bands[name] = vrt.read(1).astype("float32")
            print("  read", name, key)

    G, Rr, NIR, SWIR = bands["green"], bands["red"], bands["nir"], bands["swir"]
    eps = 1e-6
    ndwi = (G - NIR) / (G + NIR + eps)          # McFeeters: open water
    ndvi = (NIR - Rr) / (NIR + Rr + eps)        # vegetation
    ndmi = (NIR - SWIR) / (NIR + SWIR + eps)    # veg/soil moisture (marsh)

    prof = dict(driver="GTiff", height=dem_h, width=dem_w, count=6,
                dtype="float32", crs=dem_crs, transform=dem_t,
                compress="deflate", tiled=True)
    with rasterio.open(a.out, "w", **prof) as dst:
        for i, (nm, arr) in enumerate(
                [("green", G), ("red", Rr), ("nir", NIR), ("swir", SWIR),
                 ("ndwi", ndwi), ("ndmi", ndmi)], start=1):
            dst.write(arr, i); dst.set_band_description(i, nm)
    print("wrote", a.out)

    def stat(n, v):
        print("  %-5s min %7.3f  med %7.3f  max %7.3f"
              % (n, np.nanmin(v), np.nanmedian(v), np.nanmax(v)))
    print("index stats:")
    for n, v in [("ndwi", ndwi), ("ndvi", ndvi), ("ndmi", ndmi)]:
        stat(n, v)
    print("NDWI>0 (open water) px:", int((ndwi > 0).sum()),
          " | NDMI>0.1 (wet) px:", int((ndmi > 0.1).sum()))


if __name__ == "__main__":
    main()
