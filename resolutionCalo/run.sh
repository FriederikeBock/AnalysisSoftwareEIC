inputfile=/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/SingleParticleNewGeom/output_CRH.root

if [ $1 == "energy" ] ; then
  suffix="pdf";
  root -b -x -l -q 'energyResolutionCalorimeters.C("'$inputfile'","FEMC-wMat","MA","'$suffix'",kTRUE)'
  root -b -x -l -q 'energyResolutionCalorimeters.C("'$inputfile'","BECAL-wMat","MA","'$suffix'",kTRUE)'
  root -b -x -l -q 'energyResolutionCalorimeters.C("'$inputfile'","EEMC-wMat","MA","'$suffix'",kTRUE)'

  root -b -x -l -q 'energyResolutionCalorimeters.C("'$inputfile'","LFHCAL-wMat","MA","'$suffix'",kTRUE)'
  root -b -x -l -q 'energyResolutionCalorimeters.C("'$inputfile'","CHCAL","MA","'$suffix'",kTRUE)'
  root -b -x -l -q 'energyResolutionCalorimeters.C("'$inputfile'","EHCAL-wMat","MA","'$suffix'",kTRUE)'
# 
  root -b -x -l -q 'energyResolutionCalorimeters.C("'$inputfile'","LFHCAL-FEMC","MA","'$suffix'",kTRUE)'
  root -b -x -l -q 'energyResolutionCalorimeters.C("'$inputfile'","CHCAL-comb","MA","'$suffix'",kTRUE)'
  root -b -x -l -q 'energyResolutionCalorimeters.C("'$inputfile'","EHCAL-EEMC","MA","'$suffix'",kTRUE)'

elif [ $1 == "special" ]; then
  suffix="pdf";
  inputfile=/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/SinglePionSTANDALONE/output_CRH.root
  root -b -x -l -q 'energyResolutionCalorimeters.C("'$inputfile'","LFHCAL","MA","'$suffix'",kTRUE)'
  root -b -x -l -q 'energyResolutionCalorimeters.C("'$inputfile'","CHCAL","MA","'$suffix'",kTRUE)'
  root -b -x -l -q 'energyResolutionCalorimeters.C("'$inputfile'","HCALIN","MA","'$suffix'",kTRUE)'
  root -b -x -l -q 'energyResolutionCalorimeters.C("'$inputfile'","EHCALt","MA","'$suffix'",kTRUE)'

  inputfile=/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/SingleElectronStandaloneWMaterial_new/output_CRH.root
  inputfileWO=/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/SingleElectronStandaloneWOMaterial_new/output_CRH.root
  root -b -x -l -q 'energyResolutionCalorimeters.C("'$inputfile'","FEMC-wMat","MA","'$suffix'",kTRUE)'
  root -b -x -l -q 'energyResolutionCalorimeters.C("'$inputfile'","BECAL-wMat","MA","'$suffix'",kTRUE)'
  root -b -x -l -q 'energyResolutionCalorimeters.C("'$inputfile'","EEMC-wMat","MA","'$suffix'",kTRUE)'
  root -b -x -l -q 'energyResolutionCalorimeters.C("'$inputfileWO'","FEMC","MA","'$suffix'",kTRUE)'
  root -b -x -l -q 'energyResolutionCalorimeters.C("'$inputfileWO'","BECAL","MA","'$suffix'",kTRUE)'
  root -b -x -l -q 'energyResolutionCalorimeters.C("'$inputfileWO'","EEMC","MA","'$suffix'",kTRUE)'
  
elif [ $1 == "position" ];  then
  root -b -x -q -l 'positionResolutionCalorimeters.C("'$inputfile'","LFHCAL-wMat","MA","pdf",kFALSE)'
  root -b -x -q -l 'positionResolutionCalorimeters.C("'$inputfile'","EHCAL-wMat","MA","pdf",kFALSE)'
  root -b -x -q -l 'positionResolutionCalorimeters.C("'$inputfile'","CHCAL-wMat","MA","pdf",kFALSE)'
  root -b -x -q -l 'positionResolutionCalorimeters.C("'$inputfile'","EEMC-wMat","MA","pdf",kFALSE)'
  root -b -x -q -l 'positionResolutionCalorimeters.C("'$inputfile'","FEMC-wMat","MA","pdf",kFALSE)'
  root -b -x -q -l 'positionResolutionCalorimeters.C("'$inputfile'","BECAL-wMat","MA","pdf",kFALSE)'
elif [ $1 == "material" ];  then
  root -b -x -l -q 'compare_EnergyResol.C("compareMaterialBECAL.txt","MaterialStudies","BECAL","pdf",1)'
  root -b -x -l -q 'compare_EnergyResol.C("compareMaterialBECAL.txt","MaterialStudies","BECAL","pdf",3)'
  root -b -x -l -q 'compare_EnergyResol.C("compareMaterialBECAL.txt","MaterialStudies","BECAL","pdf",15)'
  root -b -x -l -q 'compare_EnergyResol.C("compareMaterialBECAL.txt","MaterialStudies","BECAL","pdf",20)'
  root -b -x -l -q 'compare_EnergyResol.C("compareMaterialEEMC.txt","MaterialStudies","EEMC","pdf",1)'
  root -b -x -l -q 'compare_EnergyResol.C("compareMaterialEEMC.txt","MaterialStudies","EEMC","pdf",3)'
  root -b -x -l -q 'compare_EnergyResol.C("compareMaterialEEMC.txt","MaterialStudies","EEMC","pdf",15)'
  root -b -x -l -q 'compare_EnergyResol.C("compareMaterialEEMC.txt","MaterialStudies","EEMC","pdf",20)'
  root -b -x -l -q 'compare_EnergyResol.C("compareMaterialFEMC.txt","MaterialStudies","FEMC","pdf",1)'
  root -b -x -l -q 'compare_EnergyResol.C("compareMaterialFEMC.txt","MaterialStudies","FEMC","pdf",3)'
  root -b -x -l -q 'compare_EnergyResol.C("compareMaterialFEMC.txt","MaterialStudies","FEMC","pdf",15)'
  root -b -x -l -q 'compare_EnergyResol.C("compareMaterialFEMC.txt","MaterialStudies","FEMC","pdf",20)'

fi
