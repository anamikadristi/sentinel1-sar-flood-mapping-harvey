#!/usr/bin/env bash
# Run the full pipeline. Large downloads and intermediate rasters go to
# $HARVEY_CACHE (default: ./cache inside the repo, which is git-ignored).
set -euo pipefail
cd "$(dirname "$0")"
for s in 01_download_data.py 02_preprocess_sar.py 03_hand_dem.py 04_flood_detection.py \
         05_validation.py 06_exposure_stats.py 07_figures.py; do
  echo "=== $s"
  python -I "$s"
done
