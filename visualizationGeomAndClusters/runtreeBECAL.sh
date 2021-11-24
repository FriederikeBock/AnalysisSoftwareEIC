# for i in {0..100};
# do 
#   root -x -q -l -b 'plot2DBecalClusters.C+("/home/fbock/eic/Singularity/Fun4All_G4_FullDetectorModular/1000_BECAL-SingleGamma_12.5GeV/G4EICDetector_eventtree.root","pdf","SingleGamma",12.5,'$i',"")'
# done

if [ $1 == "BARREL" ]; then
  for i in {0..50};
  # for i in {0..10};
  do 
  #   root -b -x -q -l 'plotBarrelClusters3D.C("/media/fbock/Samsung_T5/simulationOutputEIC/ECCE2ndCampaign/particleGun/singlePion/eval_00001/DST_General_particleGun_singlePion_000_1666000_02000_g4event_eval.root","/media/fbock/Samsung_T5/simulationOutputEIC/ECCE2ndCampaign/geometry2ndCampaign.root","pdf","SinglePion",-1,'$i',7)'
#     root -b -x -q -l 'plotBarrelClusters3D.C("/media/fbock/Samsung_T5/simulationOutputEIC/ECCE2ndCampaign/pythia8/ep-18x275-q2-100/eval_00002/DST_HFandJets_pythia8_ep-18x275-q2-100_000_2202000_02000_g4event_eval.root","/media/fbock/Samsung_T5/simulationOutputEIC/ECCE2ndCampaign/geometry2ndCampaign.root","pdf","PythiaQ2_100",-1,'$i',7)'
#     root -b -x -q -l 'plotBarrelClusters3D.C("/media/fbock/Samsung_T5/simulationOutputEIC/ECCE2ndCampaign/pythia8/ep-18x275-q2-100/eval_00002/DST_HFandJets_pythia8_ep-18x275-q2-100_000_2202000_02000_g4event_eval.root","/media/fbock/Samsung_T5/simulationOutputEIC/ECCE2ndCampaign/geometry2ndCampaign.root","pdf","PythiaQ2_100",-1,'$i',6)'
#     root -b -x -q -l 'plotBarrelClusters3D.C("/media/fbock/Samsung_T5/simulationOutputEIC/ECCE2ndCampaign/pythia8/ep-18x275-q2-100/eval_00002/DST_HFandJets_pythia8_ep-18x275-q2-100_000_2202000_02000_g4event_eval.root","/media/fbock/Samsung_T5/simulationOutputEIC/ECCE2ndCampaign/geometry2ndCampaign.root","pdf","PythiaQ2_100",-1,'$i',10)'
#     root -b -x -q -l 'plotBarrelClusters3D.C("/media/fbock/Samsung_T5/simulationOutputEIC/ECCESimulations/SingleParticleProductions/output_TTLGEO_7_HITS_MultFix_SimpleMultiPion/output_TTLGEO_7_HITS_MultFix_SimpleMultiPion.root","/media/fbock/Samsung_T5/simulationOutputEIC/ECCE2ndCampaign/geometry2ndCampaign.root","pdf","SinglePion",-1,'$i',10)'
#     root -b -x -q -l 'plotBarrelClusters3D.C("/media/fbock/Samsung_T5/simulationOutputEIC/ECCESimulations/SingleParticleProductions/output_TTLGEO_7_HITS_MultFix_SimpleMultiPion/output_TTLGEO_7_HITS_MultFix_SimpleMultiPion.root","/media/fbock/Samsung_T5/simulationOutputEIC/ECCE2ndCampaign/geometry2ndCampaign.root","pdf","SinglePion",-1,'$i',7)'
#     root -b -x -q -l 'plotBarrelClusters3D.C("/media/fbock/Samsung_T5/simulationOutputEIC/ECCESimulations/SingleParticleProductions/output_TTLGEO_7_HITS_MultFix_SimpleMultiPion/output_TTLGEO_7_HITS_MultFix_SimpleMultiPion.root","/media/fbock/Samsung_T5/simulationOutputEIC/ECCE2ndCampaign/geometry2ndCampaign.root","pdf","SinglePion",-1,'$i',6)'
    
    root -b -x -q -l 'plotBarrelClusters3D.C("/media/fbock/Samsung_T5/simulationOutputEIC/ECCESimulations/SingleParticleProductions/merged_output_TTLGEO_7_HITS_SimplePiZero/merged_output_TTLGEO_7_HITS_SimplePiZero.root","/media/fbock/Samsung_T5/simulationOutputEIC/ECCE2ndCampaign/geometry2ndCampaign.root","pdf","SinglePi0",-1,'$i',10)'
    
    
  done
elif [ $1 == "BCKW" ]; then
 for i in {0..100};
#  for i in {21..100};
  do 
#     root -b -x -q -l 'plotClustersEndCapCalo.C("/media/fbock/Samsung_T5/simulationOutputEIC/ECCE2ndCampaign/pythia8/ep-18x275-q2-100/eval_00002/DST_HFandJets_pythia8_ep-18x275-q2-100_000_2202000_02000_g4event_eval.root","/media/fbock/Samsung_T5/simulationOutputEIC/ECCE2ndCampaign/geometry2ndCampaign.root","pdf","PythiaQ2_100",-1,'$i',3)'
#     root -b -x -q -l 'plotClustersEndCapCalo.C("/media/fbock/Samsung_T5/simulationOutputEIC/ECCE2ndCampaign/pythia8/ep-18x275-q2-100/eval_00002/DST_HFandJets_pythia8_ep-18x275-q2-100_000_2202000_02000_g4event_eval.root","/media/fbock/Samsung_T5/simulationOutputEIC/ECCE2ndCampaign/geometry2ndCampaign.root","pdf","PythiaQ2_100",-1,'$i',5)'
#     root -b -x -q -l 'plotClustersEndCapCalo.C("/media/fbock/Samsung_T5/simulationOutputEIC/ECCESimulations/Sartre-ePb-18x108-e-IP8_10k_eventeval.root","/media/fbock/Samsung_T5/simulationOutputEIC/ECCE2ndCampaign/geometry2ndCampaign.root","pdf","Satre18x108",-1,'$i',3)'
#     root -b -x -q -l 'plotClustersEndCapCalo.C("/media/fbock/Samsung_T5/simulationOutputEIC/ECCE2ndCampaign/particleGun/singleElectron/eval_00001/singleE.root","/media/fbock/Samsung_T5/simulationOutputEIC/ECCE2ndCampaign/geometry2ndCampaign.root","pdf","SingleE",-1,'$i',3)'
#     root -b -x -q -l 'plotClustersEndCapCalo.C("/media/fbock/Samsung_T5/simulationOutputEIC/ECCESimulations/fixed_energy_productions_Cadesrunning/output_TTLGEO_7_HITS_13GeV_SimpleMultiElectron/output_TTLGEO_7_HITS_13GeV_SimpleMultiElectron.root","/media/fbock/Samsung_T5/simulationOutputEIC/ECCE2ndCampaign/geometry2ndCampaign.root","pdf","SingleEWMat",-1,'$i',3)'

    root -b -x -q -l 'plotClustersEndCapCalo.C("/media/fbock/Samsung_T5/simulationOutputEIC/ECCESimulations/SingleParticleProductions/merged_output_TTLGEO_7_HITS_SimplePiZero/merged_output_TTLGEO_7_HITS_SimplePiZero.root","/media/fbock/Samsung_T5/simulationOutputEIC/ECCE2ndCampaign/geometry2ndCampaign.root","pdf","SingleEWMat",-1,'$i',3)'


  done
elif [ $1 == "FWD" ]; then
 for i in {0..400};
  do 
#     root -b -x -q -l 'plotClustersEndCapCalo.C("/media/fbock/Samsung_T5/simulationOutputEIC/ECCE2ndCampaign/pythia8/ep-18x275-q2-100/eval_00002/DST_HFandJets_pythia8_ep-18x275-q2-100_000_2202000_02000_g4event_eval.root","/media/fbock/Samsung_T5/simulationOutputEIC/ECCE2ndCampaign/geometry2ndCampaign.root","pdf","PythiaQ2_100",-1,'$i',1)'
#     root -b -x -q -l 'plotClustersEndCapCalo.C("/media/fbock/Samsung_T5/simulationOutputEIC/ECCE2ndCampaign/particleGun/singleElectron/eval_00001/singleE.root","/media/fbock/Samsung_T5/simulationOutputEIC/ECCE2ndCampaign/geometry2ndCampaign.root","pdf","SingleE",-1,'$i',1)'
    root -b -x -q -l 'plotClustersEndCapCalo.C("/media/fbock/Samsung_T5/simulationOutputEIC/ECCESimulations/SingleParticleProductions/merged_output_TTLGEO_7_HITS_SimplePiZero/merged_output_TTLGEO_7_HITS_SimplePiZero.root","/media/fbock/Samsung_T5/simulationOutputEIC/ECCE2ndCampaign/geometry2ndCampaign.root","pdf","SingleEWMat",-1,'$i',1)'
  done
fi

