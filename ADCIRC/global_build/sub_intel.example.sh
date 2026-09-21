#!/usr/bin/env bash
# Copy to sub_intel.sh and adapt the partition/module names for your cluster.
# ADCIRC_BIN must contain executable adcprep and padcirc files.

#SBATCH --nodes=1
#SBATCH --ntasks=12

set -euo pipefail

: "${ADCIRC_BIN:?Set ADCIRC_BIN to the directory containing adcprep and padcirc}"
NPROC="${SLURM_NTASKS:-12}"
MPI_LAUNCHER="${MPI_LAUNCHER:-mpirun}"

rm -rf PE[0-9][0-9][0-9][0-9]
rm -f partmesh.txt metis_graph.txt fort.16 fort.18 fort.33 fort.61 fort.62 \
  fort.63 fort.64 fort.67 fort.68 fort.71 fort.72 fort.73 fort.74 fort.80 \
  adcprep.log padcirc_log.txt

"${ADCIRC_BIN}/adcprep" --np "$NPROC" --partmesh > adcprep.log
"${ADCIRC_BIN}/adcprep" --np "$NPROC" --prepall >> adcprep.log

if grep -Rqi "DOES NOT LIE WITHIN ANY ELEMENT" fort.16 PE*/fort.16 2>/dev/null; then
  echo "ERROR: at least one recording station is outside the mesh." >&2
  exit 2
fi

"$MPI_LAUNCHER" -np "$NPROC" "${ADCIRC_BIN}/padcirc" > padcirc_log.txt
