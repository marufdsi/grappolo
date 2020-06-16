#!/bin/bash

#export OMP_NUM_THREADS=18
#export KMP_AFFINITY="explicit,proclist=[0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34],verbose"

dr="datasets"
f=5
for g in "Cage15.graph" "uk2002.graph" "333SP.graph" "AS365.graph" "eu-2005.graph" "in-2004.graph" "M6.graph" "NLR.graph" "road_usa.graph" "af_shell10.graph" "cnr-2000.graph" "G_n_pin_pout.graph" "NACA0015.graph" "road_central.graph" "smallworld.graph"
do
    ./driverForGraphClustering -f $f  $dr"/"$g -b 0
done

snapf=8

for g in "Amazon0505.txt" "as-skitter.txt" "ca-HepPh.txt" "ca-HepTh.txt" "cit-Patents.txt" "com-lj.ungraph.txt" "com-youtube.ungraph.txt" "roadNet-CA.txt" "web-Google.txt" "com-orkut.ungraph.txt" "soc-LiveJournal1" "wiki-Talk" "web-Stanford"
do
    ./driverForGraphClustering -f $snapf  $dr"/"$g -b 0
done
