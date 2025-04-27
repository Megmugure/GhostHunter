#!/bin/bash
# IMa3_run_2pop.sh - Run IMa3 for 2-population model on a PBS cluster

#PBS -V
#PBS -l nodes=1:ppn=20:mpi
#PBS -N IMa3_2pop
#PBS -joe
#PBS -q batch
#PBS -l walltime=500:00:00

cd $PBS_O_WORKDIR
NCORES=$(wc -w < $PBS_NODEFILE)

echo "Running 2-population IMa3 model with $NCORES cores"

for FILE in *_2pop.u
do
    echo "Running IMa3 on $FILE"
    mpirun -np $NCORES IMa3 -i "$FILE" \
        -o "${FILE%.u}.out" \
        -q100 -m5.5 -t5.5 \
        -b10000 -L5000 -d200 -p 2 -r245 \
        -hn 120 -ha 0.99 -hb 0.9
done

