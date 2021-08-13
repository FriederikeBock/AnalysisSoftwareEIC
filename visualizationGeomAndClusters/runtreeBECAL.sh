# for i in {0..100};
# do 
#   root -x -q -l -b 'plot2DBecalClusters.C+("/home/fbock/eic/Singularity/Fun4All_G4_FullDetectorModular/1000_BECAL-SingleGamma_12.5GeV/G4EICDetector_eventtree.root","pdf","SingleGamma",12.5,'$i',"")'
# done

for i in {0..100};
do 
  root -x -q -l -b 'plot2DBecalClusters.C+("/media/fbock/2fba62ae-79f7-469c-bca3-385d6f6261da/eic-outputs/BECALcorr/1000_BECAL-SingleElectron_100GeV/G4EICDetector_eventtree.root","pdf","SingleElectron",100,'$i',"")'
done
