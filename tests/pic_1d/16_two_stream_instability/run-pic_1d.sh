#!/bin/sh

if [[ $# -eq 0 ]] ; then
    echo 'The first argument should be the number of MPI processes.'
    exit 1
fi

n=$1
target=pic_1d
MPIEXEC=mpiexec

ninja $target || exit 1
find ../data -type f -print -delete || exit 1
$MPIEXEC -n $n ./src/$target/$target --wd ../data -record_particle_at_init -save --outer_Nt=5000 --load=false || exit 1
$MPIEXEC -n $n ./src/$target/$target --wd ../data -record_particle_at_init -save --outer_Nt=5000 --load=true  || exit 1
