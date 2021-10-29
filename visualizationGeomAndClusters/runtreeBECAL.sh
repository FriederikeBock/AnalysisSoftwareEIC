# for i in {0..100};
# do 
#   root -x -q -l -b 'plot2DBecalClusters.C+("/home/fbock/eic/Singularity/Fun4All_G4_FullDetectorModular/1000_BECAL-SingleGamma_12.5GeV/G4EICDetector_eventtree.root","pdf","SingleGamma",12.5,'$i',"")'
# done

for i in {0..100};
# for i in {0..10};
do 
#   root -b -x -q -l 'plotBarrelClusters3D.C("/media/fbock/Samsung_T5/simulationOutputEIC/ECCE2ndCampaign/particleGun/singlePion/eval_00001/DST_General_particleGun_singlePion_000_1666000_02000_g4event_eval.root","/media/fbock/Samsung_T5/simulationOutputEIC/ECCE2ndCampaign/geometry2ndCampaign.root","pdf","SinglePion",-1,'$i',7)'
  root -b -x -q -l 'plotBarrelClusters3D.C("/media/fbock/Samsung_T5/simulationOutputEIC/ECCE2ndCampaign/pythia8/ep-18x275-q2-100/eval_00002/DST_HFandJets_pythia8_ep-18x275-q2-100_000_2202000_02000_g4event_eval.root","/media/fbock/Samsung_T5/simulationOutputEIC/ECCE2ndCampaign/geometry2ndCampaign.root","pdf","PythiaQ2_100",-1,'$i',7)'
  root -b -x -q -l 'plotBarrelClusters3D.C("/media/fbock/Samsung_T5/simulationOutputEIC/ECCE2ndCampaign/pythia8/ep-18x275-q2-100/eval_00002/DST_HFandJets_pythia8_ep-18x275-q2-100_000_2202000_02000_g4event_eval.root","/media/fbock/Samsung_T5/simulationOutputEIC/ECCE2ndCampaign/geometry2ndCampaign.root","pdf","PythiaQ2_100",-1,'$i',6)'
  root -b -x -q -l 'plotBarrelClusters3D.C("/media/fbock/Samsung_T5/simulationOutputEIC/ECCE2ndCampaign/pythia8/ep-18x275-q2-100/eval_00002/DST_HFandJets_pythia8_ep-18x275-q2-100_000_2202000_02000_g4event_eval.root","/media/fbock/Samsung_T5/simulationOutputEIC/ECCE2ndCampaign/geometry2ndCampaign.root","pdf","PythiaQ2_100",-1,'$i',10)'
  
  
done
