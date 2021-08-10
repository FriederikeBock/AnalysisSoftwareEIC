for i in {0..100};
do 
  root -b -x -q -l 'plot2DHCalInClusters.C+("/media/fbock/2fba62ae-79f7-469c-bca3-385d6f6261da/eic-outputs/CHCAL/1000_CHCAL-SinglePion_25GeV/G4EICDetector_eventtree.root","pdf","SinglePion",25,'$i')'
  root -b -x -q -l 'plot2DHCalOutClusters.C+("/media/fbock/2fba62ae-79f7-469c-bca3-385d6f6261da/eic-outputs/CHCAL/1000_CHCAL-SinglePion_25GeV/G4EICDetector_eventtree.root","pdf","SinglePion",25,'$i')'
done
