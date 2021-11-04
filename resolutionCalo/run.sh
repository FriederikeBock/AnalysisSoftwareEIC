if [ $1 == "energy" ] ; then
  # root -b -x -l -q 'energyResolutionCalorimeters.C("/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/SingleElectron500K/output_CRH.root","FEMC-wMat","MA","pdf",kFALSE)'
  # root -b -x -l -q 'energyResolutionCalorimeters.C("/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/SingleElectron500K/output_CRH.root","BECAL-wMat","MA","pdf",kFALSE)'
  # root -b -x -l -q 'energyResolutionCalorimeters.C("/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/SingleElectron500K/output_CRH.root","EEMC-wMat","MA","pdf",kFALSE)'

  # root -b -x -l -q 'energyResolutionCalorimeters.C("/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/SingleElectron500K/output_CRH.root","FEMC-wMat","MA","pdf",kTRUE)'
  # root -b -x -l -q 'energyResolutionCalorimeters.C("/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/SingleElectron500K/output_CRH.root","BECAL-wMat","MA","pdf",kTRUE)'
  # root -b -x -l -q 'energyResolutionCalorimeters.C("/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/SingleElectron500K/output_CRH.root","EEMC-wMat","MA","pdf",kTRUE)'

  suffix="pdf";
#   suffix="png";
  
  root -b -x -l -q 'energyResolutionCalorimeters.C("/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/SinglePartileFullSummed/output_CRH.root","FEMC-wMat","MA","'$suffix'",kTRUE)'
  root -b -x -l -q 'energyResolutionCalorimeters.C("/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/SinglePartileFullSummed/output_CRH.root","BECAL-wMat","MA","'$suffix'",kTRUE)'
  root -b -x -l -q 'energyResolutionCalorimeters.C("/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/SinglePartileFullSummed/output_CRH.root","EEMC-wMat","MA","'$suffix'",kTRUE)'

  root -b -x -l -q 'energyResolutionCalorimeters.C("/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/SinglePartileFullSummed/output_CRH.root","LFHCAL-wMat","MA","'$suffix'",kTRUE)'
  root -b -x -l -q 'energyResolutionCalorimeters.C("/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/SinglePartileFullSummed/output_CRH.root","CHCAL","MA","'$suffix'",kTRUE)'
  root -b -x -l -q 'energyResolutionCalorimeters.C("/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/SinglePartileFullSummed/output_CRH.root","HCALIN-wMat","MA","'$suffix'",kTRUE)'
  root -b -x -l -q 'energyResolutionCalorimeters.C("/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/SinglePartileFullSummed/output_CRH.root","EHCAL-wMat","MA","'$suffix'",kTRUE)'

  root -b -x -l -q 'energyResolutionCalorimeters.C("/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/SinglePartileFullSummed/output_CRH.root","LFHCAL-FEMC","MA","'$suffix'",kTRUE)'
  root -b -x -l -q 'energyResolutionCalorimeters.C("/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/SinglePartileFullSummed/output_CRH.root","CHCAL-comb","MA","'$suffix'",kTRUE)'
  root -b -x -l -q 'energyResolutionCalorimeters.C("/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/SinglePartileFullSummed/output_CRH.root","EHCAL-EEMC","MA","'$suffix'",kTRUE)'
elif [ $1 == "position" ];  then
  root -b -x -q -l 'positionResolutionCalorimeters.C("/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/SinglePartileFullSummed/output_CRH.root","LFHCAL-wMat","MA","pdf",kFALSE)'
  root -b -x -q -l 'positionResolutionCalorimeters.C("/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/SinglePartileFullSummed/output_CRH.root","EHCAL-wMat","MA","pdf",kFALSE)'
  root -b -x -q -l 'positionResolutionCalorimeters.C("/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/SinglePartileFullSummed/output_CRH.root","CHCAL-wMat","MA","pdf",kFALSE)'
  root -b -x -q -l 'positionResolutionCalorimeters.C("/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/SinglePartileFullSummed/output_CRH.root","EEMC-wMat","MA","pdf",kFALSE)'
  root -b -x -q -l 'positionResolutionCalorimeters.C("/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/SinglePartileFullSummed/output_CRH.root","FEMC-wMat","MA","pdf",kFALSE)'
  root -b -x -q -l 'positionResolutionCalorimeters.C("/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/SinglePartileFullSummed/output_CRH.root","BECAL-wMat","MA","pdf",kFALSE)'
fi
