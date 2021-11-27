# inputFileTR="/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/SingleParticleMultClusterizers/output_TRKEFF.root"
# inputFileCS="/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/SingleParticleMultClusterizers/output_CS.root"

inputFileTR="/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/SingleParticleNewGeom/output_TRKEFF.root"
inputFileCS="/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/SingleParticleNewGeom/output_CS.root"
# inputFileTR="/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/SinglePartileFullSummed/output_TRKEFF.root"
# inputFileCS="/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/SinglePartileFullSummed/output_CS.root"
  root -x -q -l -b 'clustereffi.C("'$inputFileTR'","'$inputFileCS'","LFHCAL","pdf","SinglePart","",0,0)'
  root -x -q -l -b 'clustereffi.C("'$inputFileTR'","'$inputFileCS'","HCALOUT","pdf","SinglePart","",0,0)'
  root -x -q -l -b 'clustereffi.C("'$inputFileTR'","'$inputFileCS'","HCALIN","pdf","SinglePart","",0,0)'
  root -x -q -l -b 'clustereffi.C("'$inputFileTR'","'$inputFileCS'","EHCAL","pdf","SinglePart","",0,0)'
  root -x -q -l -b 'clustereffi.C("'$inputFileTR'","'$inputFileCS'","FEMC","pdf","SinglePart","",0,0)'
  root -x -q -l -b 'clustereffi.C("'$inputFileTR'","'$inputFileCS'","BECAL","pdf","SinglePart","",0,0)'
  root -x -q -l -b 'clustereffi.C("'$inputFileTR'","'$inputFileCS'","EEMC","pdf","SinglePart","",0,0)'
