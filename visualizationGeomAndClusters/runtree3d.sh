# for i in {0..99};
# do 
#   root -x -q -l -b 'plot3DTowerstree.C("/home/fbock/eic/Singularity/Fun4All_G4_FullDetectorModular/G4EICDetector_eventtree.root","pdf","SinglePion",10,'$i')'
# done

# for i in {0..99};
# do 
#   root -x -q -l -b 'plot3DTowerstree.C("/home/fbock/eic/Singularity/Fun4All_G4_FullDetectorModular/1000_SinglePion_30GeV/G4EICDetector_eventtree.root","pdf","SinglePion",30,'$i')'
# done

# for i in {0..99};
# do 
#   root -x -q -l -b 'plot3DTowerstree.C("/home/fbock/eic/Singularity/Fun4All_G4_FullDetectorModular/1000_SinglePion_2GeV/G4EICDetector_eventtree.root","pdf","SinglePion",2,'$i')'
# done
# 
# for i in {0..99};
# do 
#   root -x -q -l -b 'plot3DTowerstree.C("/home/fbock/eic/Singularity/Fun4All_G4_FullDetectorModular/1000_SingleProton_2GeV/G4EICDetector_eventtree.root","pdf","SingleProton",2,'$i')'
# done

# for i in {0..99};
# do 
#   root -x -q -l -b 'plot3DTowerstree.C("/home/fbock/eic/Singularity/Fun4All_G4_FullDetectorModular/1000_SingleProton_30GeV/G4EICDetector_eventtree.root","pdf","SingleProton",30,'$i')'
# done

# for i in {0..99};
# do 
#   root -x -q -l -b 'plot3DTowerstree.C("/home/fbock/eic/Singularity/Fun4All_G4_FullDetectorModular/1000_SingleProton_100GeV/G4EICDetector_eventtree.root","pdf","SingleProton",100,'$i')'
# done

# for i in {0..20};
# do 
#   root -x -q -l -b 'plot3DTowerstree.C+("/media/fbock/2fba62ae-79f7-469c-bca3-385d6f6261da/eic-outputs/LFHCAL/1000_LFHCAL-SingleElectron_30GeV/G4EICDetector_eventtree.root","pdf","SingleElectron",30,'$i',true)'
# done
# 
# for i in {0..20};
# do 
#   root -x -q -l -b 'plot3DTowerstree.C+("/media/fbock/2fba62ae-79f7-469c-bca3-385d6f6261da/eic-outputs/LFHCAL/1000_LFHCAL-SingleProton_30GeV/G4EICDetector_eventtree.root","pdf","SingleProton",30,'$i',true)'
# done
# 
# for i in {0..20};
# do 
#   root -x -q -l -b 'plot3DTowerstree.C+("/media/fbock/2fba62ae-79f7-469c-bca3-385d6f6261da/eic-outputs/LFHCAL/1000_LFHCAL-SinglePion_30GeV/G4EICDetector_eventtree.root","pdf","SinglePion",30,'$i',true)'
# done
# 
# for i in {0..20};
# do 
#   root -x -q -l -b 'plot3DTowerstree.C+("/media/fbock/2fba62ae-79f7-469c-bca3-385d6f6261da/eic-outputs/LFHCAL/additional/1000_LFHCAL-SingleNeutron_30GeV/G4EICDetector_eventtree.root","pdf","SingleNeutron",30,'$i',true)'
# done
# 
# for i in {0..20};
# do 
#   root -x -q -l -b 'plot3DTowerstree.C+("/media/fbock/2fba62ae-79f7-469c-bca3-385d6f6261da/eic-outputs/LFHCAL/additional/1000_LFHCAL-SingleKaon_30GeV/G4EICDetector_eventtree.root","pdf","SingleKaon",30,'$i',true)'
# done
# 
# for i in {0..10};
# do 
#   root -x -q -l -b 'plot3DTowerstree.C+("/media/fbock/2fba62ae-79f7-469c-bca3-385d6f6261da/eic-outputs/LFHCAL/1000_LFHCAL-SingleElectron_100GeV/G4EICDetector_eventtree.root","pdf","SingleElectron",100,'$i',true)'
#   root -x -q -l -b 'plot3DTowerstree.C+("/media/fbock/2fba62ae-79f7-469c-bca3-385d6f6261da/eic-outputs/LFHCAL/1000_LFHCAL-SingleProton_100GeV/G4EICDetector_eventtree.root","pdf","SingleProton",100,'$i',true)'
#   root -x -q -l -b 'plot3DTowerstree.C+("/media/fbock/2fba62ae-79f7-469c-bca3-385d6f6261da/eic-outputs/LFHCAL/1000_LFHCAL-SinglePion_100GeV/G4EICDetector_eventtree.root","pdf","SinglePion",100,'$i',true)'
#   root -x -q -l -b 'plot3DTowerstree.C+("/media/fbock/2fba62ae-79f7-469c-bca3-385d6f6261da/eic-outputs/LFHCAL/additional/1000_LFHCAL-SingleNeutron_100GeV/G4EICDetector_eventtree.root","pdf","SingleNeutron",100,'$i',true)'
#   root -x -q -l -b 'plot3DTowerstree.C+("/media/fbock/2fba62ae-79f7-469c-bca3-385d6f6261da/eic-outputs/LFHCAL/additional/1000_LFHCAL-SingleKaon_100GeV/G4EICDetector_eventtree.root","pdf","SingleKaon",100,'$i',true)'
# done
# 
for i in {0..400};
do 
#   root -b -x -l -q 'plot3DTowerstree.C("/media/fbock/Samsung_T5/simulationOutputEIC/ECCE2ndCampaign/pythia8/ep-18x275-q2-100/eval_00002/DST_HFandJets_pythia8_ep-18x275-q2-100_000_2202000_02000_g4event_eval.root","/media/fbock/Samsung_T5/simulationOutputEIC/ECCE2ndCampaign/geometry2ndCampaign.root","pdf","PythiaQ2_100",-1,'$i')'
  root -b -x -l -q 'plot3DTowerstree.C("/media/fbock/Samsung_T5/simulationOutputEIC/ECCE2ndCampaign/particleGun/singlePion/eval_00001/singlepi.root","/media/fbock/Samsung_T5/simulationOutputEIC/ECCE2ndCampaign/geometry2ndCampaign.root","pdf","SinglePion",-1,'$i')'
done
