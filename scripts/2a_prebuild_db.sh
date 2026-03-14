#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --job-name=liftoff_prebuild_db
#SBATCH --partition=your_partition
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=32g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your@email.com
#SBATCH --output=outputs/logs/prebuild_db_%j.out
#SBATCH --error=outputs/logs/prebuild_db_%j.err
#
# Pre-builds the gffutils SQLite database from the CS GFF3 annotation ONCE,
# before submitting the Liftoff array. Without this, simultaneous array jobs
# each attempt to build the same database file concurrently, causing
# corruption and random job failures.
#
# RECOMMENDED USAGE — submit array only after database is ready:
#
#   DBID=$(sbatch --parsable 2a_prebuild_db.sh)
#   sbatch --array=0-23 --dependency=afterok:$DBID 2_run_liftoff_array.sh
#
# ALTERNATIVE — resubmit any failed array jobs individually:
#
#   sbatch --array=3,7,12 2_run_liftoff_array.sh
#
# =============================================================================

set -eo pipefail
source config.sh
eval "$(conda shell.bash hook)"
conda activate "$CONDA_ENV"

CS_DB="$OUTDIR/renamed/CS/CS.gff3.db"

echo "Pre-building gffutils database..."
echo "  Input GFF : $CS_GFF"
echo "  Output DB : $CS_DB"
echo "  Started   : $(date)"

python3 - << PYEOF
import gffutils
db = gffutils.create_db(
    "$CS_GFF",
    dbfn="$CS_DB",
    force=True,
    keep_order=True,
    merge_strategy="merge",
    sort_attribute_values=True,
    disable_infer_genes=True,
    disable_infer_transcripts=True,
)
print(f"Database built: $CS_DB")
PYEOF

echo "Done: $(date)"
echo "Database: $CS_DB"
