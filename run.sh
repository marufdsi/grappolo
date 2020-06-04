#!/bin/bash

dr="datasets"
f=5
for g in "333SP.graph" "AS365.graph" "eu-2005.graph" "in-2004.graph" "M6.graph" "NLR.graph" "road_usa.graph" "af_shell10.graph" "cnr-2000.graph" "G_n_pin_pout.graph" "NACA0015.graph" "road_central.graph" "smallworld.graph"
do
    ./driverForGraphClustering -f $f  $dr"/"$g -b 0
done
