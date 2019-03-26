#!/bin/bash
input="input"
bin="bin"
out="output/measureLS.dat"
out_temp="output/temp.dat"

iterations=1

rm ${out}
rm ${out_temp}
touch ${out_temp}
touch ${out}

echo ${iterations} >> ${out_temp}

for i in {21..30}
do
  echo "Executing MDG-a_$i"
  for k in {1..${iterations}}
  do
    ${bin}/localSearch < ${input}/MDG-a_${i}_* >> ${out_temp}
  done
done

for i in {11..20}
do
  echo "Executing SOM-b_$i"
  for k in {1..${iterations}}
  do
    ${bin}/localSearch < ${input}/SOM-b_${i}_* >> ${out_temp}
  done
done

for i in {1..10}
do
  echo "Executing GKD-c_$i"
  for k in {1..${iterations}}
  do
    ${bin}/localSearch < ${input}/GKD-c_${i}_* >> ${out_temp}
  done
done

echo "Executing stats"
${bin}/stats < ${out_temp} >> ${out}

# rm ${out_temp}

exit 0
