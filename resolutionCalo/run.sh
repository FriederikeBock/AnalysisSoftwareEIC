inputfile=/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/SingleParticleNewGeom/output_CRH.root

if [ $1 == "energy" ] ; then
  suffix="pdf";
#   inputfile=/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/SingleParticleNewGeom2/output_CRH.root
#   inputfile=/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/StandaloneSinglePiCALComb/output_CRH.root
#   inputfile=/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/SingleParticleEEMCPCarbonEOP/output_CRH.root
#   
#   root -b -x -l -q 'energyResolutionCalorimeters.C("'$inputfile'","FEMC-wMat","MA","'$suffix'",kTRUE)'
#   root -b -x -l -q 'energyResolutionCalorimeters.C("'$inputfile'","BECAL-wMat","MA","'$suffix'",kTRUE)'
#   root -b -x -l -q 'energyResolutionCalorimeters.C("'$inputfile'","EEMC-wMat","MA","'$suffix'",kTRUE)'

#   root -b -x -l -q 'energyResolutionCalorimeters.C("'$inputfile'","LFHCAL-wMat","MA","'$suffix'",kTRUE)'
#   root -b -x -l -q 'energyResolutionCalorimeters.C("'$inputfile'","CHCAL","MA","'$suffix'",kTRUE)'
#   root -b -x -l -q 'energyResolutionCalorimeters.C("'$inputfile'","EHCAL-wMat","MA","'$suffix'",kTRUE)'
# 
#   root -b -x -l -q 'energyResolutionCalorimeters.C("'$inputfile'","LFHCAL-FEMC","MA","'$suffix'",kTRUE)'
#   root -b -x -l -q 'energyResolutionCalorimeters.C("'$inputfile'","CHCAL-comb","MA","'$suffix'",kTRUE)'
#   root -b -x -l -q 'energyResolutionCalorimeters.C("'$inputfile'","EHCAL-EEMC","MA","'$suffix'",kTRUE)'

#   inputfile=/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/SingleParticle_EOP_PbGlasG4/output_CRH.root
#   inputfile=/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/SingleParticle_EOP_PbGlasTF1/output_CRH.root
#   root -b -x -l -q 'energyResolutionCalorimeters.C("'$inputfile'","BECAL-wMat","MA","'$suffix'",kTRUE)'
  
#   inputfile=/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/SingleElectron_NewBEMCGeom/output_CRH.root
#   root -b -x -l -q 'energyResolutionCalorimeters.C("'$inputfile'","BECAL-wMat","MA","'$suffix'",kTRUE)'
#   inputfile=/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/SingleElectron_NewBEMCGeom_NoGaps/output_CRH.root
  inputfile=/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/SingleElectron_NewBEMCGeom_NoGaps1cell/output_CRH.root
  root -b -x -l -q 'energyResolutionCalorimeters.C("'$inputfile'","BECAL-wMat","MA","'$suffix'",kTRUE)'
elif [ $1 == "special" ]; then
  suffix="pdf";
#   inputfile=/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/SinglePionSTANDALONE/output_CRH.root
#   inputfile=/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/StandaloneSinglePiWOCalib/output_CRH.root
#   inputfile=/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/StandaloneSinglePiEdepCalib/output_CRH.root
#   inputfile=/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/StandaloneSinglePiEatAll/output_CRH.root
#   inputfile=/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/HCALStandaloneSingleParticle/output_CRH.root
  
  
#   root -b -x -l -q 'energyResolutionCalorimeters.C("'$inputfile'","LFHCAL","MA","'$suffix'",kTRUE)'
#   root -b -x -l -q 'energyResolutionCalorimeters.C("'$inputfile'","CHCAL","MA","'$suffix'",kTRUE)'
#   root -b -x -l -q 'energyResolutionCalorimeters.C("'$inputfile'","HCALIN","MA","'$suffix'",kTRUE)'
#   root -b -x -l -q 'energyResolutionCalorimeters.C("'$inputfile'","EHCAL","MA","'$suffix'",kTRUE)'
  # proper 0.5mm carbon between towers instead of 2mm
#   inputfileWO=/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/StandaloneEEMCPCarbon/output_CRH.root
#   inputfile=/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/SingleElectronEEMCPCarbon/output_CRH.root

#   inputfile=/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/SinglePionEEMCPCarbonCalib/output_CRH.root
  inputfile=/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/SinglePionEEMCNCCalib/output_CRH.root
#   inputfileWO=/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/STANDALONEEEMC30cm/output_CRH.root
#   inputfile=/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/SingleElectronEEMC30cm/output_CRH.root
  
#   inputfile=/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/SinglePionEEMCImpCarbon_13GeV/output_CRH.root
  #no carbon between towers
#   inputfileWO=/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/STANDALONEEEMCNC/output_CRH.root
#   inputfile=/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/SingleElectronEEMCNC/output_CRH.root

#   inputfileWO=/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/SingleElectronStandaloneEdepCalib/output_CRH.root
# #   inputfileWO=/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/SingleElectronStandaloneWOCalib/output_CRH.root
# #   inputfile=/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/SingleElectronStandaloneWMaterial_new/output_CRH.root
# #   inputfileWO=/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/SingleElectronStandaloneWOMaterial_new/output_CRH.root
# #   root -b -x -l -q 'energyResolutionCalorimeters.C("'$inputfile'","FEMC-wMat","MA","'$suffix'",kTRUE)'
# #   root -b -x -l -q 'energyResolutionCalorimeters.C("'$inputfile'","BECAL-wMat","MA","'$suffix'",kTRUE)'
  root -b -x -l -q 'energyResolutionCalorimeters.C("'$inputfile'","EEMC-wMat","MA","'$suffix'",kTRUE)'
#   root -b -x -l -q 'energyResolutionCalorimeters.C("'$inputfileWO'","FEMC","MA","'$suffix'",kTRUE)'
#   root -b -x -l -q 'energyResolutionCalorimeters.C("'$inputfileWO'","BECAL","MA","'$suffix'",kTRUE)'
#   root -b -x -l -q 'energyResolutionCalorimeters.C("'$inputfileWO'","EEMC","MA","'$suffix'",kTRUE)'
  
elif [ $1 == "position" ];  then
#   inputfile=/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/SingleParticleEEMCPCarbonEOP/output_CRH.root
#   root -b -x -q -l 'positionResolutionCalorimeters.C("'$inputfile'","LFHCAL-wMat","MA","pdf",kFALSE)'
#   root -b -x -q -l 'positionResolutionCalorimeters.C("'$inputfile'","EHCAL-wMat","MA","pdf",kFALSE)'
#   root -b -x -q -l 'positionResolutionCalorimeters.C("'$inputfile'","CHCAL-wMat","MA","pdf",kFALSE)'
#   root -b -x -q -l 'positionResolutionCalorimeters.C("'$inputfile'","EEMC-wMat","MA","pdf",kFALSE)'
#   root -b -x -q -l 'positionResolutionCalorimeters.C("'$inputfile'","FEMC-wMat","MA","pdf",kFALSE)'
#   root -b -x -q -l 'positionResolutionCalorimeters.C("'$inputfile'","BECAL-wMat","MA","pdf",kFALSE)'

  inputfile=/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/SingleParticle_EOP_PbGlasG4/output_CRH.root
  root -b -x -q -l 'positionResolutionCalorimeters.C("'$inputfile'","BECAL-wMat","MA","pdf",kFALSE)'
elif [ $1 == "material" ];  then
#   root -b -x -l -q 'compare_EnergyResol.C("compareMaterialBECAL.txt","MaterialStudies","BECAL","pdf",1)'
#   root -b -x -l -q 'compare_EnergyResol.C("compareMaterialBECAL.txt","MaterialStudies","BECAL","pdf",3)'
#   root -b -x -l -q 'compare_EnergyResol.C("compareMaterialBECAL.txt","MaterialStudies","BECAL","pdf",15)'
#   root -b -x -l -q 'compare_EnergyResol.C("compareMaterialBECAL.txt","MaterialStudies","BECAL","pdf",20)'
#   root -b -x -l -q 'compare_EnergyResol.C("compareMaterialEEMC.txt","MaterialStudies","EEMC","pdf",1)'
#   root -b -x -l -q 'compare_EnergyResol.C("compareMaterialEEMC.txt","MaterialStudies","EEMC","pdf",3)'
#   root -b -x -l -q 'compare_EnergyResol.C("compareMaterialEEMC.txt","MaterialStudies","EEMC","pdf",15)'
#   root -b -x -l -q 'compare_EnergyResol.C("compareMaterialEEMC.txt","MaterialStudies","EEMC","pdf",20)'
root -b -x -l -q 'compare_EnergyResol.C("compareMaterialEEMC_matstart.txt","MaterialStudiesStart","EEMC","pdf",20)'
  root -b -x -l -q 'compare_EnergyResol.C("compareMaterialEEMC_mat.txt","MaterialStudiesSize","EEMC","pdf",20)'
    root -b -x -l -q 'compare_EnergyResol.C("compareMaterialEEMC_matProper.txt","MaterialStudiesProper","EEMC","pdf",20)'
#   root -b -x -l -q 'compare_EnergyResol.C("compareMaterialEEMC_0_5mmCarbon.txt","MaterialStudies05mmCarbon","EEMC","pdf",20)'
#   root -b -x -l -q 'compare_EnergyResol.C("compareMaterialFEMC.txt","MaterialStudies","FEMC","pdf",1)'
#   root -b -x -l -q 'compare_EnergyResol.C("compareMaterialFEMC.txt","MaterialStudies","FEMC","pdf",3)'
#   root -b -x -l -q 'compare_EnergyResol.C("compareMaterialFEMC.txt","MaterialStudies","FEMC","pdf",15)'
#   root -b -x -l -q 'compare_EnergyResol.C("compareMaterialFEMC.txt","MaterialStudies","FEMC","pdf",20)'
elif [ $1 == "compareReso" ]; then
#   root -b -x -q -l 'energyResolutionCalorimeters_comparePlot.C("compareEnergyResolECals.txt","ECal","pdf")'
#   root -b -x -q -l 'energyResolutionCalorimeters_comparePlot.C("compareEnergyResolBEMC.txt","BECal","pdf")'
#   root -b -x -q -l 'energyResolutionCalorimeters_comparePlot.C("compareEnergyResolHCals.txt","HCal","pdf")'
#   root -b -x -l -q 'compare_EnergyResol.C("compareMaterialGlass_BECAL.txt","GlassStudy","BECAL","pdf",1)'
#   root -b -x -l -q 'compare_EnergyResol.C("compareMaterialGlass_BECAL.txt","GlassStudy","BECAL","pdf",3)'
#   root -b -x -l -q 'compare_EnergyResol.C("compareMaterialGlass_BECAL.txt","GlassStudy","BECAL","pdf",15)'
#   root -b -x -l -q 'compare_EnergyResol.C("compareMaterialGlass_BECAL.txt","GlassStudy","BECAL","pdf",20)'

  root -b -x -l -q 'compare_EnergyResol.C("compareMaterialGlass_BECALNewGeom.txt","GeomComparsionWOGaps","BECAL","pdf",1)'
  root -b -x -l -q 'compare_EnergyResol.C("compareMaterialGlass_BECALNewGeom.txt","GeomComparsionWOGaps","BECAL","pdf",3)'
  root -b -x -l -q 'compare_EnergyResol.C("compareMaterialGlass_BECALNewGeom.txt","GeomComparsionWOGaps","BECAL","pdf",15)'
  root -b -x -l -q 'compare_EnergyResol.C("compareMaterialGlass_BECALNewGeom.txt","GeomComparsionWOGaps","BECAL","pdf",20)'

  
  
fi
