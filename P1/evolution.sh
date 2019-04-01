#!/bin/bash
input="input"
bin="bin"
out="output/evolution"

algs=(localSearch localSearchDet)
files=(MDG-a_21_n2000_m200.txt SOM-b_13_n400_m40.txt)
for alg in "${algs[@]}"
do
  for file in "${files[@]}"
  do
    echo "Executing $alg"
    rm $out/${alg}Evolution.dat
    $bin/$alg < $input/$file > $out/${alg}Evolution${file}.dat
  done
done

exit 0
