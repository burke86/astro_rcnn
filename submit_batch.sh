sbatch --partition=solo --nodes=4 --ntasks-per-node=12 --cores-per-socket=12 --gres=gpu:v100:16 --mem-per-cpu=1500 --time=12:00:00 --export=ALL ./astro_rcnn.sh
