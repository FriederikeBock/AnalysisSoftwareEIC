for i in {0..100};
do 
  root -x -q -l -b 'plotBarrelClusters3D.C+("/home/ewa/alice/EIC/G4EICDetector_eventtree.root","/home/ewa/alice/EIC/geometry_ideal_BECALproj_EEMC_LFHCAL_corrected.root","pdf","SingleElectron",12.5,'$i',6,0,"")'
done
