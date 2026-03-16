# CS Coordinate Converter — Usage Guide

## Available targets
- Jagger
- Lancer
- ArinaLrFor
- Stanley
- Spelt
- Mace
- SY_Mattis
- Julius
- Landmark
- Norin61

## Python

```python
from cs_coordinate_converter import CSCoordinateConverter

conv = CSCoordinateConverter("outputs/conversion_tables")

# Single position
tgt_chr, tgt_pos, density = conv.convert("chr1A", 500_000_000, "Mace")
print(f"chr1A:500Mb -> {tgt_chr}:{tgt_pos} (density={density})")

# Region
tgt_chr, start, end, density = conv.convert_region("chr1A", 490e6, 510e6, "Mace")

# BED file
df = conv.convert_bed("my_regions.bed", "SY_Mattis", output_file="converted.tsv")

# List targets
print(conv.available_targets())
```

## CLI (Python)
```bash
python cs_coordinate_converter.py --tables-dir outputs/conversion_tables \
    --target Mace --chr chr1A --pos 500000000

python cs_coordinate_converter.py --tables-dir outputs/conversion_tables \
    --target SY_Mattis --bed regions.bed --output converted.tsv
```

## R
```r
source("cs_coordinate_converter.R")
conv <- load_cs_converter("outputs/conversion_tables")

# Single position
res <- cs_convert(conv, "chr1A", 500e6, "Mace")

# Region
res <- cs_convert_region(conv, "chr1A", 490e6, 510e6, "Mace")

# BED data frame
result_df <- cs_convert_bed(conv, my_bed_df, "SY_Mattis")
```

## Confidence/quality
Each conversion returns `anchor_density` = number of gene anchors within ±5 Mb.
- **≥10**: high confidence
- **5–10**: moderate confidence (sparse region, possibly pericentromeric)
- **<5**: low confidence (extrapolated or poorly mapped region)
