# inputfile=/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/SinglePartileFullSummed/output_CRH.root
inputfile=/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/SingleParticleNew/output_CRH.root

if [ $1 == "energy" ] ; then
  # root -b -x -l -q 'energyResolutionCalorimeters.C("/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/SingleElectron500K/output_CRH.root","FEMC-wMat","MA","pdf",kFALSE)'
  # root -b -x -l -q 'energyResolutionCalorimeters.C("/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/SingleElectron500K/output_CRH.root","BECAL-wMat","MA","pdf",kFALSE)'
  # root -b -x -l -q 'energyResolutionCalorimeters.C("/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/SingleElectron500K/output_CRH.root","EEMC-wMat","MA","pdf",kFALSE)'

  # root -b -x -l -q 'energyResolutionCalorimeters.C("/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/SingleElectron500K/output_CRH.root","FEMC-wMat","MA","pdf",kTRUE)'
  # root -b -x -l -q 'energyResolutionCalorimeters.C("/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/SingleElectron500K/output_CRH.root","BECAL-wMat","MA","pdf",kTRUE)'
  # root -b -x -l -q 'energyResolutionCalorimeters.C("/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/SingleElectron500K/output_CRH.root","EEMC-wMat","MA","pdf",kTRUE)'

  suffix="pdf";
#   suffix="png";
  
  root -b -x -l -q 'energyResolutionCalorimeters.C("'$inputfile'","FEMC-wMat","MA","'$suffix'",kTRUE)'
  root -b -x -l -q 'energyResolutionCalorimeters.C("'$inputfile'","BECAL-wMat","MA","'$suffix'",kTRUE)'
  root -b -x -l -q 'energyResolutionCalorimeters.C("'$inputfile'","EEMC-wMat","MA","'$suffix'",kTRUE)'

  root -b -x -l -q 'energyResolutionCalorimeters.C("'$inputfile'","LFHCAL-wMat","MA","'$suffix'",kTRUE)'
  root -b -x -l -q 'energyResolutionCalorimeters.C("'$inputfile'","CHCAL","MA","'$suffix'",kTRUE)'
  root -b -x -l -q 'energyResolutionCalorimeters.C("'$inputfile'","HCALIN-wMat","MA","'$suffix'",kTRUE)'
  root -b -x -l -q 'energyResolutionCalorimeters.C("'$inputfile'","EHCAL-wMat","MA","'$suffix'",kTRUE)'

  root -b -x -l -q 'energyResolutionCalorimeters.C("'$inputfile'","LFHCAL-FEMC","MA","'$suffix'",kTRUE)'
  root -b -x -l -q 'energyResolutionCalorimeters.C("'$inputfile'","CHCAL-comb","MA","'$suffix'",kTRUE)'
  root -b -x -l -q 'energyResolutionCalorimeters.C("'$inputfile'","EHCAL-EEMC","MA","'$suffix'",kTRUE)'
elif [ $1 == "position" ];  then
  root -b -x -q -l 'positionResolutionCalorimeters.C("'$inputfile'","LFHCAL-wMat","MA","pdf",kFALSE)'
  root -b -x -q -l 'positionResolutionCalorimeters.C("'$inputfile'","EHCAL-wMat","MA","pdf",kFALSE)'
  root -b -x -q -l 'positionResolutionCalorimeters.C("'$inputfile'","CHCAL-wMat","MA","pdf",kFALSE)'
  root -b -x -q -l 'positionResolutionCalorimeters.C("'$inputfile'","EEMC-wMat","MA","pdf",kFALSE)'
  root -b -x -q -l 'positionResolutionCalorimeters.C("'$inputfile'","FEMC-wMat","MA","pdf",kFALSE)'
  root -b -x -q -l 'positionResolutionCalorimeters.C("'$inputfile'","BECAL-wMat","MA","pdf",kFALSE)'
fi
