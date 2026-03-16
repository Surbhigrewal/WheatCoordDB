# WheatCoordDB

**Continuous coordinate conversion from IWGSC CS RefSeq v2.1 to 24 hexaploid wheat assemblies.**

🌾 **Live tool:** https://surbhigrewal.github.io/WheatCoordDB/

---

## What it does

WheatCoordDB converts any genomic position (QTL intervals, KASP marker positions, introgression breakpoints, or arbitrary regions) from IWGSC Chinese Spring RefSeq v2.1 coordinates to equivalent positions across 24 chromosome-scale wheat assemblies. Unlike gene-homology tools, WheatCoordDB supports continuous coordinate interpolation: any position, genic or intergenic, can be converted.
Conversion uses approximately 94,000 orthologous gene anchors per assembly (projected by Liftoff) as reference points, with PCHIP spline interpolation at 1 kb resolution across all 21 wheat chromosomes. An anchor density score is returned with every query as a confidence proxy.

---

## Use the web tool

No installation needed. Visit **https://surbhigrewal.github.io/WheatCoordDB/**

- Single region, BED batch upload, and multi-assembly comparison
- All computation runs in your browser.
- Translocation-aware: ArinaLrFor and SY_Mattis 5B/7B queries return 
  split results with correct inverted-orientation coordinates
- Interactive synteny dotplot per query, downloadable as PNG

---

## Run the pipeline for your own species

The pipeline is species-agnostic. Any organism with a reference genome 
annotation and one or more target assemblies can be processed. It has been 
designed for large polyploid genomes but works equally well for barley, 
Arabidopsis, rice, or any other species.

**Requirements:**
- Reference genome FASTA + GFF3 annotation
- One or more target assembly FASTAs
- Conda (mamba recommended)
- SLURM HPC (or adapt the shell scripts to run locally)

**Setup:**
```bash
git clone https://github.com/Surbhigrewal/WheatCoordDB
cd WheatCoordDB/scripts

# 1. Edit config.sh to set paths to your reference and target assemblies
#    and update assembly names/FASTAs in 1_run_rename_all.sh

# 2. Create the conda environment
bash 0_run_setup_conda_env.sh

# 3. Standardise chromosome names
sbatch --array=0-N 1_run_rename_all.sh

# 4. Pre-build the annotation database (prevents array job conflicts)
DBID=$(sbatch --parsable 2a_prebuild_db.sh)

# 5. Run Liftoff + anchor extraction + conversion tables
sbatch --array=0-N --dependency=afterok:$DBID 2_run_liftoff_array.sh

# 6. Generate Python and R converter scripts
python3 03_generate_converter.py \
    --conversion-dir outputs/conversion_tables \
    --target-names Jagger Lancer Mace \  # list your assemblies
    --output-dir outputs/converter

# 7. Compile master files
python3 04_postprocess.py
```

**Programmatic access** 
First generate the converter scripts using step 5 above, then:

**(Python):**
```python
from cs_coordinate_converter import CSCoordinateConverter

conv = CSCoordinateConverter('outputs/conversion_tables')

# Single position
tgt_chr, tgt_pos, density = conv.convert('Chr5A', 557_234_900, 'Kariega')

# Region
tgt_chr, tgt_start, tgt_end, density = conv.convert_region('Chr1B', 500e6, 550e6, 'Mace')

# BED file
df = conv.convert_bed('my_qtl_intervals.bed', 'Jagger', output_file='converted.tsv')

# List available assemblies
print(conv.available_targets())
```

**(R):**
```r
source("cs_coordinate_converter.R")
conv <- load_cs_converter("/path/to/outputs/conversion_tables")
result <- cs_convert(conv, "Chr5A", 557234900, "Kariega")
result <- cs_convert_region(conv, "Chr1B", 500e6, 550e6, "Mace")
result_df <- cs_convert_bed(conv, my_bed_df, "Jagger")
```

---
**Confidence guide:**

| anchor_density | Confidence | Typical region |
|---|---|---|
| 10 or above | High | Gene-rich chromosome arms |
| 5 to 10 | Moderate | Interstitial regions |
| Below 5 | Low | Pericentromeric heterochromatin |

## Repository contents
```
scripts/                          # Pipeline scripts (species-agnostic)
  config.sh                       # Paths and parameters — edit this first
  0_run_setup_conda_env.sh        # Create conda environment
  rename_chromosomes.py           # Standardise chromosome naming
  1_run_rename_all.sh             # SLURM: rename all assemblies
  2a_prebuild_db.sh               # Pre-build gffutils database
  2_run_liftoff_array.sh          # SLURM array: Liftoff + anchors + tables
  01_extract_anchors.py           # Anchor pair extraction
  02_build_conversion_table.py    # PCHIP spline construction
  03_generate_converter.py        # Generate Python class + R functions
  04_postprocess.py               # Master workbook compilation
  translocation_map.tsv           # Translocation breakpoint definitions
data/anchors/                     # Per-assembly anchor TSV files (24 assemblies)
index.html                        # Web application (GitHub Pages)
outputs/
  conversion_tables/      504 gzip-compressed conversion tables (21 chr x 24 assemblies)
  plots/                  528 synteny dotplots (21-chr overviews + per-chromosome)
  master_conversion_summary.xlsx   Per-chromosome summary across all assemblies

# master_anchors.xlsx is not included due to file size.
# Regenerate with: python3 scripts/04_postprocess.py
```

---

## Assemblies included

| Group | Assemblies |
|---|---|
| 10-wheat panel | Jagger, Lancer, ArinaLrFor, Stanley, Spelt, Mace, SY_Mattis, Julius, Landmark, Norin61 |
| CS versions | CS_v1, CS_IAAS, CS_CAU |
| Additional cultivars | Kariega, Renan_v2, Fielder, Attraktion, Paragon_v3, Cadenza_v2, Aikang58, Chunmai104, Sumai3, JIN50, MOV |

---

## Citation

> Grewal, S. WheatCoordDB: a gene-anchored coordinate conversion resource 
> for hexaploid wheat assemblies. *Scientific Data* (in preparation).

---

## Dependencies

Liftoff v1.6.3 · minimap2 v2.24 · Python 3.10 (numpy, pandas, scipy, 
matplotlib, openpyxl) · R (data.table)

---

## AI disclosure

Pipeline development, script writing, and web application implementation were conducted with assistance from Claude AI (Anthropic). All scientific  decisions, assembly selection, biological interpretation, feature specification, site design, and validation were performed by the author.

---

*Built by the [Grewal lab](https://grewallab.com/), University of Nottingham.*
