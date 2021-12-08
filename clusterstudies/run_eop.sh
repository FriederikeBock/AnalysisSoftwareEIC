#input=/media/nschmidt/local/AnalysisSoftwareEIC/treeAnalysis/treeProcessing/CENTRALSIM_PYTHIA_TTL8_EoP_MB/output_TMSTUD.root

input_elec=/media/nschmidt/local/AnalysisSoftwareEIC/treeAnalysis/treeProcessing/CADES_SINGLEMULTIELECTRON_TTLGEO7_EOP2/output_EOPSTUD.root
input_pion=/media/nschmidt/local/AnalysisSoftwareEIC/treeAnalysis/treeProcessing/CADES_SINGLEMULTIPION_TTLGEO7_EOP2/output_EOPSTUD.root
if [ $1 = "fbock" ]; then
  input_elec=/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/SingleElectronEEMCPCarbonEOP/output_EOPSTUD.root
#   input_pion=/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/SinglePionEEMCPCarbonEOP/output_EOPSTUD_merged.root
  input_pion=/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/SinglePionEEMCPCarbonEOP/output_EOPSTUD_merged_low_mid.root
#   input_pion=/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/SinglePionEEMCPCarbonEOP/output_EOPSTUD_merged_all.root
fi

root -x -l -b -q 'eoverpstudies_YR.C+("'$input_pion'","'$input_elec'","pdf",true,"SinglePart")'
