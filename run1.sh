cp datasets/karate.graph datasets/karateApprox100.graph
cp datasets/karate.graph datasets/karateApprox90.graph
cp datasets/karate.graph datasets/karateApprox80.graph
cp datasets/karate.graph datasets/karateApprox70.graph
cp datasets/karate.graph datasets/karateMap.graph
cp datasets/karate.graph datasets/karateNoMap.graph
for i in {1..10};
do
#echo $i
         
#numactl --interleave=all

  ./driverForGraphClusteringApprox -p 100  -f 5  datasets/karateApprox100.graph -b 0 -o >> karate.graph_clustInfoApprox100.out$i
mv karateApprox100.graph_clustInfo results/karateApprox100.nocolor.rand.out$i
  ./driverForGraphClusteringApprox -p 90  -f 5  datasets/karateApprox90.graph -b 0 -o >> karate.graph_clustInfoApprox90.out$i
mv karateApprox90.graph_clustInfo results/karateApprox90.nocolor.rand.out$i
 ./driverForGraphClusteringApprox -p 80  -f 5  datasets/karateApprox80.graph -b 0 -o >> karate.graph_clustInfoApprox80.out$i
mv karateApprox80.graph_clustInfo results/karateApprox80.nocolor.rand.out$i
  ./driverForGraphClusteringApprox -p 70  -f 5  datasets/karateApprox70.graph -b 0 -o >> karate.graph_clustInfoApprox70.out$i
mv karateApprox70.graph_clustInfo results/karateApprox70.nocolor.rand.out$i

done
./driverForGraphClustering -f 5  datasets/karateMap.graph -b 0 -o >> karate.graph_clustInfoMap.out
mv karateMap.graph_clustInfo results/karateMap.nocolor.rand.out1
  ./driverForGraphClustering -f 5  datasets/karateNoMap.graph -b 1 -o >> karate.graph_clustInfoNoMap.out
mv karateNoMap.graph_clustInfo results/karateNoMap.nocolor.rand.out1 
