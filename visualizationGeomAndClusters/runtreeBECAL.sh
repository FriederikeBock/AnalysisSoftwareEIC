for i in {0..100};
do 
  root -x -q -l -b 'plot2DBecalClusters.C+("/home/fbock/eic/Singularity/Fun4All_G4_FullDetectorModular/1000_BECAL-SingleGamma_12.5GeV/G4EICDetector_eventtree.root","pdf","SingleGamma",12.5,'$i',"")'
done
