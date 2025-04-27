#!/bin/bash
# IMa3_run_3pop_ghost.sh - Run IMa3 for 3-population model (ghost) on a PBS cluster

#PBS -V
#PBS -l nodes=2:ppn=20:mpi
#PBS -N IMa3_3pop_ghost
#PBS -joe
#PBS -q batch
#PBS -l walltime=500:00:00

cd $PBS_O_WORKDIR
NCORES=$(wc -w < $PBS_NODEFILE)

echo "Running 3-population IMa3 model (with ghost) using $NCORES cores"

for FILE in *.u
do
    echo "Running IMa3 on $FILE (ghost model)"
    mpirun -np $NCORES IMa3 -i "$FILE" \
        -o "${FILE%.u}_ghost.out" \
        -q100 -m5.5 -t5.5 \
        -b100000 -L5000 -d200 -p 2 -r245 \
        -hn 120 -ha 0.99 -hb 0.9 -j1
done

