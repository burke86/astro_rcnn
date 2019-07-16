sbatch --partition=solo --nodes=1 --ntasks-per-node=144 --sockets-per-node=1 --cores-per-socket=18 --threads-per-core=4 --mem-per-cpu=1500 --time=72:00:00 --export=ALL ../simulate.sh
