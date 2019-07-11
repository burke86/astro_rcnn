srun --partition=solo  --pty --nodes=1 --ntasks-per-node=12 --cores-per-socket=12 --gres=gpu:v100:1 --mem-per-cpu=1500 --time=72:00:00 --wait=0 --export=ALL /bin/bash
