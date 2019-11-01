#!/bin/sh
#PBS -q workq
#PBS -o stdout.$PBS_JOBID
#PBS -e stderr.$PBS_JOBID
#PBS -V
#PBS -N <name>
#PBS -l nodes=1:ppn=32

cd $PBS_O_WORKDIR
cp $PBS_NODEFILE nodelist.$PBS_JOBID

DIR=$HOME/git-flac/flac/src

### Recording state of the code
cp $DIR/snapshot.diff .

ulimit -c unlimited
### Execute the model
export OMP_NUM_THREADS=$PBS_NP
$DIR/flac subduction.inp.s.1
cp _contents.save _contents.rs

#$DIR/flac subduction.inp.s.2
#cp _contents.save _contents.rs

~/git-flac/flac/util/flac2vtk.py ./

# ~~~~ submit command ~~~~
# qsub < [script]

