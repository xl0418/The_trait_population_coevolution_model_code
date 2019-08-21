#!/bin/bash

#gamma=(0 0.001 0.01 0.1 0.5 1)
#a=(0 0.001 0.01 0.1 0.5 1)

for((s = 0; s <= 1; s++))
do
for((m = 0; m <= 1; m++))
do
echo "#!/bin/bash" > zMLjob$s$m
echo "#SBATCH --time=05:59:00" >> zMLjob$s$m
echo "#SBATCH --ntasks=1" >> zMLjob$s$m
echo "#SBATCH --nodes=1" >> zMLjob$s$m
echo "#SBATCH --cpus-per-task=24" >> zMLjob$s$m
echo "module load Python" >> zMLjob$s$m
echo "srun python3 evo.py "$s" "$m"" >> zMLjob$s$m
echo "rm zMLjob$s$m" >> zMLjob$s$m

#sbatch --partition=regular --mem=12GB --job-name=ML$s --mail-type=FAIL,TIME_LIMIT --mail-user=glaudanno@gmail.com zMLjob$s
sbatch --partition=regular --mem-per-cpu=1GB --job-name=smc$s$m --output=smc$s$m.log --mail-type=FAIL,TIME_LIMIT --mail-user=xl0418@gmail.com zMLjob$s$m

done
done
