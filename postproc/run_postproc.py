import sys
import subprocess

pid = subprocess.Popen(['sbatch', 'DPS_VSCI_postproc.sh'])
# pid = subprocess.Popen(['sbatch', 'DPS_VSVI_postproc.sh'])
pid = subprocess.Popen(['sbatch', 'VSCI_CSCI_postproc.sh'])
pid = subprocess.Popen(['sbatch', 'VSVI_CSVI_postproc.sh'])
