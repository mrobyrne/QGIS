"""
esri_tiles.py  --  Download Esri World Imagery XYZ tiles for an AOI and
mosaic them into a georeferenced GeoTIFF.

Imagery (c) Esri, Maxar, Earthstar Geographics, and the GIS User Community.
Intended for internal analysis / water-body confirmation against the DEM.

Usage (driven by the DEM extent):
    python esri_tiles.py --dem topo.tif --zoom 16 --out esri_rgb.tif
    python esri_tiles.py --dem topo.tif --zoom 16 --plan          # dry-run, just report
    python esri_tiles.py --dem topo.tif --zoom 16 --to-dem-crs    # also reproject to DEM CRS

Or import and call download_xyz_mosaic(...) directly.
"""
import os
import math
import time
import argparse
import io

import numpy as np
import requests
from PIL import Image
import rasterio
from rasterio.transform import Affine
from rasterio.warp import transform_bounds, calculate_default_transform, reproject, Resampling

# Esri World Imagery basemap tiles: {z}/{y}/{x}  (web mercator, top-left origin)
ESRI_URL = ("https://server.arcgisonline.com/ArcGIS/rest/services/"
            "World_Imagery/MapServer/tile/{z}/{y}/{x}")
R = 20037508.342789244          # web-mercator half-extent (m)
TILE = 256
MAX_TILES = 6000                # safety guard
HEADERS = {"User-Agent": "qgis-waterbody-tool/1.0"}


def lonlat_to_tile(lon, lat, z):
    n = 2 ** z
    x = int((lon + 180.0) / 360.0 * n)
    lat_r = math.radians(lat)
    y = int((1.0 - math.asinh(math.tan(lat_r)) / math.pi) / 2.0 * n)
    return max(0, min(n - 1, x)), max(0, min(n - 1, y))


def res_at(zoom, lat):
    """ground resolution (m/px) at given latitude."""
    return 156543.03392 * math.cos(math.radians(lat)) / (2 ** zoom)


def download_xyz_mosaic(bounds_4326, zoom, out_tif, url_template=ESRI_URL,
                        cache_dir=None, to_crs=None, dst_res=None, plan_only=False):
    """bounds_4326 = (west, south, east, north) in lon/lat."""
    w, s, e, n = bounds_4326
    xmin, ymax = lonlat_to_tile(w, n, zoom)     # NW corner -> min x, min y
    xmax, ymin_ = lonlat_to_tile(e, s, zoom)    # SE corner -> max x, max y
    if xmax < xmin:
        xmin, xmax = xmax, xmin
    if ymin_ < ymax:
        ymax, ymin_ = ymin_, ymax
    nx, ny = (xmax - xmin + 1), (ymin_ - ymax + 1)
    ntiles = nx * ny
    midlat = (s + n) / 2.0
    print(f"zoom {zoom}: {nx} x {ny} = {ntiles} tiles, "
          f"~{res_at(zoom, midlat):.2f} m/px, "
          f"mosaic {nx*TILE} x {ny*TILE}px")
    if plan_only:
        return None
    if ntiles > MAX_TILES:
        raise SystemExit(f"{ntiles} tiles exceeds MAX_TILES={MAX_TILES}; "
                         f"lower the zoom or shrink the AOI.")

    if cache_dir is None:
        cache_dir = os.path.join(os.path.dirname(os.path.abspath(out_tif)),
                                 f"_tilecache_z{zoom}")
    os.makedirs(cache_dir, exist_ok=True)

    mosaic = np.zeros((ny * TILE, nx * TILE, 3), dtype=np.uint8)
    sess = requests.Session(); sess.headers.update(HEADERS)
    got = 0
    for ix in range(xmin, xmax + 1):
        for iy in range(ymax, ymin_ + 1):
            cpath = os.path.join(cache_dir, f"{zoom}_{iy}_{ix}.jpg")
            if os.path.exists(cpath):
                data = open(cpath, "rb").read()
            else:
                url = url_template.format(z=zoom, x=ix, y=iy)
                for attempt in range(4):
                    try:
                        r = sess.get(url, timeout=30)
                        if r.status_code == 200 and r.content:
                            data = r.content
                            open(cpath, "wb").write(data)
                            break
                    except requests.RequestException:
                        pass
                    time.sleep(0.5 * (attempt + 1))
                else:
                    print(f"  !! failed {zoom}/{iy}/{ix}, leaving black")
                    continue
                time.sleep(0.05)        # be polite
            try:
                img = np.asarray(Image.open(io.BytesIO(data)).convert("RGB"))
            except Exception:
                continue
            py, px = (iy - ymax) * TILE, (ix - xmin) * TILE
            mosaic[py:py + TILE, px:px + TILE] = img
            got += 1
            if got % 100 == 0:
                print(f"  {got}/{ntiles} tiles")
    print(f"  fetched {got}/{ntiles} tiles")

    # georeference in EPSG:3857
    ts = 2 * R / (2 ** zoom)                    # tile size in meters
    px_size = ts / TILE
    x0 = -R + xmin * ts                         # left edge
    y0 = R - ymax * ts                          # top edge
    transform = Affine(px_size, 0, x0, 0, -px_size, y0)

    profile = dict(driver="GTiff", height=mosaic.shape[0], width=mosaic.shape[1],
                   count=3, dtype="uint8", crs="EPSG:3857", transform=transform,
                   compress="deflate", photometric="rgb", tiled=True)
    with rasterio.open(out_tif, "w", **profile) as dst:
        for b in range(3):
            dst.write(mosaic[:, :, b], b + 1)
    print("wrote", out_tif)

    if to_crs is not None:
        warp_to(out_tif, to_crs, dst_res, out_tif.replace(".tif", "_reproj.tif"))
    return out_tif


def warp_to(src_tif, dst_crs, dst_res, out_tif):
    with rasterio.open(src_tif) as src:
        kw = dict(resolution=dst_res) if dst_res else {}
        transform, width, height = calculate_default_transform(
            src.crs, dst_crs, src.width, src.height, *src.bounds, **kw)
        prof = src.profile.copy()
        prof.update(crs=dst_crs, transform=transform, width=width, height=height)
        with rasterio.open(out_tif, "w", **prof) as dst:
            for b in range(1, src.count + 1):
                reproject(rasterio.band(src, b), rasterio.band(dst, b),
                          src_transform=src.transform, src_crs=src.crs,
                          dst_transform=transform, dst_crs=dst_crs,
                          resampling=Resampling.bilinear)
    print("wrote", out_tif)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dem", required=True, help="reference DEM (sets AOI + CRS)")
    ap.add_argument("--zoom", type=int, default=16)
    ap.add_argument("--out", default="esri_rgb.tif")
    ap.add_argument("--plan", action="store_true", help="report tile count only")
    ap.add_argument("--to-dem-crs", action="store_true",
                    help="also write a copy reprojected to the DEM CRS")
    a = ap.parse_args()

    with rasterio.open(a.dem) as ds:
        b4326 = transform_bounds(ds.crs, "EPSG:4326", *ds.bounds)
        dem_crs = ds.crs
    out = a.out if os.path.isabs(a.out) else os.path.join(
        os.path.dirname(os.path.abspath(a.dem)), a.out)
    download_xyz_mosaic(b4326, a.zoom, out, plan_only=a.plan,
                        to_crs=dem_crs if a.to_dem_crs else None)


if __name__ == "__main__":
    main()
