#!/bin/bash
#SBATCH --job-name=nf-t2t-polish
#SBATCH --partition=norm
#SBATCH --cpus-per-task=2
#SBATCH --mem=8g
#SBATCH --time=10-00:00:00
#SBATCH --chdir=.
#SBATCH --output=logs/nextflow_%j.log
#SBATCH --error=logs/nextflow_%j.log

# ---------------------------------------------------------------------------
# Nextflow head-job wrapper for T2T-Polish on Biowulf.
#
# Usage:
#   sbatch run.sh <config_file> [-resume]
#
# The Nextflow process itself is lightweight (2 CPU / 8 GB is plenty).
# It submits and monitors all pipeline jobs via Slurm — it must stay alive
# for the entire duration of the run, hence the 3-day wall time.
# Expand if needed.
#
# Flags:
#   -resume          — resume from the last successful checkpoint [OPTIONAL]
#   -ansi-log false  — plain line-by-line output (readable in log files)
#   -with-trace      — TSV of per-task CPU/mem/time/I/O usage (timestamped filename)
#   -with-report     — HTML summary report (timestamped filename)
#   -with-dag        — workflow DAG (timestamped filename; requires graphviz for .svg)
# ---------------------------------------------------------------------------

set -euo pipefail

mkdir -p logs

# Graceful shutdown: forward SIGTERM (Slurm walltime/scancel) to Nextflow
# so it can flush the cache before exiting.  Without this, a walltime kill
# corrupts the cache the same way an interactive Ctrl+C does.
_term() {
    echo "Caught SIGTERM — sending SIGINT to Nextflow for graceful shutdown..."
    kill -INT "$NF_PID" 2>/dev/null
    wait "$NF_PID"
}
trap _term SIGTERM

module load nextflow || true
module load singularity || true   # needed if any process uses a container

echo "Usage: sbatch run.sh <config_file> [-resume]"

if [ $# -lt 1 ]; then
    echo "Error: Missing config file argument."
    exit 1
fi

ts=$(date +%Y%m%d_%H%M%S)

nextflow run $tools/T2T-Polish/nextflow/main.nf \
    -c $1 \
    ${2:-} \
    -ansi-log false \
    -with-trace  logs/trace_${ts}.txt \
    -with-report logs/report_${ts}.html \
    -with-dag    logs/dag_${ts}.svg &
NF_PID=$!
wait "$NF_PID"
