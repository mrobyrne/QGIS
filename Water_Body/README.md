# Automatic Water-Body Detection from a DEM + Imagery

A pipeline for detecting water bodies (open water **and** vegetated marsh) from a
gridded elevation model, public Esri imagery, and free Sentinel-2 NIR — built to
replace a colour-picker approach that produced many false positives (dark roofs,
shadows, asphalt) and missed marshy wetlands.

Developed/tested on a 10 m bare-earth DEM near Calgary, Alberta (EPSG:32611).

---

## Why this approach

A single signal is not enough:

| Signal | Catches | Fails on |
|--------|---------|----------|
| Imagery colour alone | water-coloured pixels | dark roofs, shadows, asphalt (false positives) |
| DEM shape alone | flat enclosed basins | graded lots, quarries, benches look identical to ponds |
| **Colour ∩ DEM basin** | real open ponds, **0 false positives** | rivers (not closed basins), **marshy/vegetated water** |
| **+ Sentinel-2 NIR (NDWI/NDMI)** | river, open ponds, **marshes** | features < ~30 m (10 m resolution floor) |

The final detector combines them:

- **Open water** = Sentinel-2 **NDWI > 0** (reliable; also catches rivers).
- **Marsh** = a flat DEM **depression** that is **vegetated and wet** —
  `NDMI_median > 0.18 AND NDVI_median > 0.40`. NDMI is the key discriminator:
  wet vegetation (marsh) ≈ 0.20–0.32, dry meadow ≈ 0.0–0.15, while both have
  strongly negative NDWI — so NDWI alone misses vegetated marsh.

---

## Requirements

Python 3.11+ with:

```
pip install rasterio geopandas shapely scikit-image scipy matplotlib requests pillow
```

`rasterio`'s pip wheel bundles its own GDAL (3.12+), so no separate GDAL/QGIS
Python install is needed. No POSIX shell required (developed on Windows /
PowerShell).

**Input:** a single-band elevation GeoTIFF named `topo.tif` (any projected CRS).
It is **not** included in this repo — drop your own DEM in the folder.

---

## Reusable modules

### `esri_tiles.py` — download Esri World Imagery (RGB)

Downloads XYZ tiles for the DEM's extent and mosaics them into a georeferenced
GeoTIFF. Tiles are disk-cached, so re-runs are free.

```bash
python esri_tiles.py --dem topo.tif --zoom 16 --plan          # dry-run: report tile count
python esri_tiles.py --dem topo.tif --zoom 16 --out esri_rgb.tif --to-dem-crs
```

- `--zoom` controls resolution (z16 ≈ 1.5 m/px at 51°N; z15 ≈ 3 m).
- `--plan` reports tile count / mosaic size without downloading.
- `--to-dem-crs` also writes a copy reprojected to the DEM CRS (pixel-usable
  with the DEM).

> Imagery © Esri, Maxar, Earthstar Geographics. For internal analysis only —
> do not redistribute the downloaded tiles.

### `get_nir.py` — fetch free Sentinel-2 NIR/SWIR

Queries the **AWS Earth Search STAC** (public Cloud-Optimized GeoTIFFs, no auth),
picks the least-cloudy scene in a date window, and warps Green/Red/NIR/SWIR onto
the **exact DEM grid**, also writing NDWI and NDMI.

```bash
python get_nir.py --dem topo.tif --out sentinel2_aoi.tif
python get_nir.py --dem topo.tif --datetime "2024-06-01T00:00:00Z/2025-09-30T23:59:59Z" --max-cloud 15
```

Output `sentinel2_aoi.tif` has 6 bands: `green, red, nir, swir, ndwi, ndmi`.

---

## Pipeline (run in order)

| Script | Purpose | Output |
|--------|---------|--------|
| `01_inspect.py` | DEM metadata + elevation stats | console |
| `02_diagnose.py` | hillshade, NoData voids, slope, candidate water | `diag.png` |
| `03_extract.py` | flat + depression + single-level region extraction | `water_candidates.gpkg`, `detected.png` |
| `04_strict.py` | high-confidence basins + zoom panels | `water_highconf.gpkg`, `top_basins.png` |
| `esri_tiles.py` | download Esri RGB mosaic | `esri_rgb.tif`, `esri_rgb_reproj.tif` |
| `05_verify.py` | overlay basins on imagery | `verify.png` |
| `06_sample.py` | calibrate RGB water-colour thresholds | console |
| `07_confirm.py` | RGB colour ∩ DEM basin (precision detector) | `water_final.gpkg`, `final.png` |
| `08_verify_final.py` | verify final RGB polygons on imagery | `verify_final.png` |
| `get_nir.py` | fetch Sentinel-2 NIR/SWIR | `sentinel2_aoi.tif` |
| `09_s2_detect.py` | first Sentinel-2 pass (open water + marsh) | `water_s2.gpkg`, `s2_detect.png` |
| `10_basin_s2_stats.py` | per-basin S2 signature (threshold calibration) | console |
| `11_s2_basin_detect.py` | **final detector**: NDWI open water + NDMI marsh | `water_s2_v2.gpkg`, `verify_marsh.png` |

**Primary output: `water_s2_v2.gpkg`** — polygons classed `open_water` / `marsh`,
loadable directly in QGIS.

---

## Key tunables

In `11_s2_basin_detect.py`:

```python
OPEN_NDWI   = 0.0    # NDWI threshold for open water
MARSH_NDMI  = 0.18   # min NDMI (moisture) for marsh
MARSH_NDVI  = 0.40   # min NDVI (vegetation) for marsh
OPEN_FRAC   = 0.30   # basin counted as open water above this NDWI>0 fraction
MIN_AREA_M2 = 800    # drop specks
```

DEM basin shape (`03_extract.py`): `SLOPE_MAX`, `SINK_MIN`, `RANGE_MAX`, `AREA_MIN`.
RGB colour rule (`07_confirm.py`): `BR_MIN`, `EXG_MAX`, `BRIGHT_MAX`.

Thresholds are scene-/sensor-specific — re-check with `06_sample.py` (RGB) and
`10_basin_s2_stats.py` (Sentinel-2) when moving to a new region.

---

## Limitations

- **Resolution floor:** 10 m DEM / Sentinel-2 → water bodies under ~30 m are missed.
- **Marsh boundaries are blocky** — they follow the 10 m DEM basin outline, not
  the true marsh edge. Refine by trimming to NDMI pixels if needed.
- **Single-date Sentinel-2** — a late-summer scene under-represents seasonal
  marsh. Multi-temporal max-NDWI / flooding-frequency is the durable upgrade.
- **NoData voids** in this DEM were data gaps, not water (region-dependent).
- Marsh class can occasionally include a lush irrigated field sitting in a
  depression — spot-check per region.

## Possible next steps

- Cross-check against authoritative water/wetland vectors (NHN, Ducks Unlimited
  Canada, Canadian Wetland Inventory).
- Multi-temporal Sentinel-2 stack for seasonal flooding frequency.
- Topographic Wetness Index (TWI) / HAND from the DEM to catch wet flats that
  are not closed depressions.
- Package `02`–`11` into a single `detect_water(dem, zoom)` entry point.
