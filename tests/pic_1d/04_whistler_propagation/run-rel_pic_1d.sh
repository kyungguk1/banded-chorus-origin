#!/bin/sh

if [[ $# -eq 0 ]] ; then
    echo 'The first argument should be the number of MPI processes.'
    exit 1
fi

n=$1
target=rel_pic_1d

ninja $target || exit 1
find ../data -type f -print -delete || exit 1
/Library/com.kyungguk.usr.localized/mpich/3.2.1/bin/mpiexec -n $n ./src/$target/$target --wd ../data -record_particle_at_init -save --outer_Nt 320 --load=false || exit 1
for i in {2..5}; do
    /Library/com.kyungguk.usr.localized/mpich/3.2.1/bin/mpiexec -n $n ./src/$target/$target --wd ../data -record_particle_at_init -save --outer_Nt 320 --load=true || exit 1
done
