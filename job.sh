#!/bin/bash                                                                                                                                                                                                                                                                                                                                                               

qsub -q copperhead -N "grappolo_clustering" -l walltime=24:00:00 -l nodes=1:ppn=$1:skylake -d $(pwd) run.sh
