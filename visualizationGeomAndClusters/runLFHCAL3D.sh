for i in {0..50};
do 
  root -x -q -l -b 'plotLFHCal3D.C+("/home/ewa/alice/EIC/DST_HFandJets_pythia8_ep-18x275-q2-100_000_0124000_02000_g4event_eval.root","/home/ewa/alice/EIC/geometry_ideal_BECALproj_EEMC_LFHCAL_corrected.root","pdf","PythiaMB",-1.0,'$i',0,"")'
done
