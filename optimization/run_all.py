import sys
import subprocess

n_eps = int(sys.argv[1]) # num of eps
n_rho = int(sys.argv[2]) # num of rho
shaper = str(sys.argv[3])

for i_eps in range(n_eps):
    for i_rho in range(n_rho):
        pid = subprocess.Popen(['sbatch', 'shaper_opt.sh', str(i_eps), str(i_rho), shaper])
