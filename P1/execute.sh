#!/bin/bash
input="input"
bin="bin"
out_local="output/localSearch.dat"
out_greedy="output/greedy.dat"

rm ${out_local}
rm ${out_greedy}
touch ${out_local}
touch ${out_greedy}

for i in {1..10}
do
  echo "Executing GKD-c_$i"
  echo "Executing GKD-c_$i" >> ${out_local}
  ${bin}/p1 < ${input}/GKD-c_${i}_* >> ${out_local}
done

for i in {11..20}
do
  echo "Executing SOM-b_$i"
  echo "Executing SOM-b_$i" >> ${out_local}
  ${bin}/p1 < ${input}/SOM-b_${i}_* >> ${out_local}
done

for i in {21..30}
do
  echo "Executing MDG-a_$i"
  echo "Executing MDG-a_$i" >> ${out_local}
  ${bin}/p1 < ${input}/MDG-a_${i}_* >> ${out_local}
done

exit 0