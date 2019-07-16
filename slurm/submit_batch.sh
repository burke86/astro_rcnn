sbatch --partition=solo --nodes=1 --ntasks-per-node=12 --cores-per-socket=12 --gres=gpu:v100:4 --mem-per-cpu=1500 --time=72:00:00 --export=ALL ../astro_rcnn.sh
