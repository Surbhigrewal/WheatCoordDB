# CS Coordinate Converter — Usage Guide

Generated from 24 target assemblies.

## Available targets
- Aikang58
- ArinaLrFor
- Attraktion
- CS_CAU
- CS_IAAS
- CS_v1
- Cadenza_v2
- Chuanmai104
- Fielder
- JIN50
- Jagger
- Julius
- Kariega
- Lancer
- Landmark
- MOV
- Mace
- Norin61
- Paragon_v3
- Renan_v2
- SY_Mattis
- Spelt
- Stanley
- Sumai3

## Confidence metric — anchor_fraction (0.0–1.0)
Proportion of CS v2.1 genes within ±5 Mb recovered as anchors in the target assembly.
Normalised by CS gene density so pericentromeric regions are not unfairly penalised.

| anchor_fraction | Confidence | Interpretation |
|---|---|---|
| ≥ 0.8 | High | Well-anchored, reliable interpolation |
| 0.5–0.8 | Moderate | Reduced gene recovery — introgression, assembly gap, or diverged region |
| < 0.5 | Low | Poor recovery — coordinates unreliable |

## Python

```python
from cs_coordinate_converter import CSCoordinateConverter

conv = CSCoordinateConverter("outputs/conversion_tables")

# Single position
tgt_chr, tgt_pos, conf = conv.convert("Chr1A", 500_000_000, "Mace")
print(f"Chr1A:500Mb -> {tgt_chr}:{tgt_pos} ({conf:.0%} gene recovery)")

# Region
tgt_chr, start, end, conf = conv.convert_region("Chr1A", 490e6, 510e6, "Mace")

# BED file
df = conv.convert_bed("my_regions.bed", "SY_Mattis", output_file="converted.tsv")

# List targets
print(conv.available_targets())
```

## CLI
```bash
# Single position
python cs_coordinate_converter.py --tables-dir outputs/conversion_tables \
    --target Mace --chr Chr1A --pos 500000000

# Region
python cs_coordinate_converter.py --tables-dir outputs/conversion_tables \
    --target Mace --chr Chr1A --start 490000000 --end 510000000

# BED file
python cs_coordinate_converter.py --tables-dir outputs/conversion_tables \
    --target SY_Mattis --bed regions.bed --output converted.tsv
```

## R
```r
source("cs_coordinate_converter.R")
conv <- load_cs_converter("outputs/conversion_tables")

# Single position
res <- cs_convert(conv, "Chr1A", 500e6, "Mace")
cat(res$tgt_chr, res$tgt_pos, res$anchor_fraction, "\n")

# Region
res <- cs_convert_region(conv, "Chr1A", 490e6, 510e6, "Mace")

# BED data frame
result_df <- cs_convert_bed(conv, my_bed_df, "SY_Mattis")
```
