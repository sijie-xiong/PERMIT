import sys
import subprocess

i_eps = int(sys.argv[1]) # num of eps
i_rho = int(sys.argv[2]) # num of rho
shaper = str(sys.argv[3])

pid = subprocess.Popen(['sbatch', 'shaper_opt.sh', str(i_eps), str(i_rho), shaper])
