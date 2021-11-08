inputFileTR="/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/SingleParticleNew/output_TRKEFF.root"
inputFileCS="/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/SingleParticleNew/output_CS.root"
# inputFileTR="/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/SinglePartileFullSummed/output_TRKEFF.root"
# inputFileCS="/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/SinglePartileFullSummed/output_CS.root"
  root -x -q -l -b 'clustereffi.C("'$inputFileTR'","'$inputFileCS'","LFHCAL","pdf","SingleParticle","",0,0)'
  root -x -q -l -b 'clustereffi.C("'$inputFileTR'","'$inputFileCS'","HCALOUT","pdf","SingleParticle","",0,0)'
  root -x -q -l -b 'clustereffi.C("'$inputFileTR'","'$inputFileCS'","HCALIN","pdf","SingleParticle","",0,0)'
  root -x -q -l -b 'clustereffi.C("'$inputFileTR'","'$inputFileCS'","EHCAL","pdf","SingleParticle","",0,0)'
  root -x -q -l -b 'clustereffi.C("'$inputFileTR'","'$inputFileCS'","FEMC","pdf","SingleParticle","",0,0)'
  root -x -q -l -b 'clustereffi.C("'$inputFileTR'","'$inputFileCS'","BECAL","pdf","SingleParticle","",0,0)'
  root -x -q -l -b 'clustereffi.C("'$inputFileTR'","'$inputFileCS'","EEMC","pdf","SingleParticle","",0,0)'
