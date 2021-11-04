input=/media/nschmidt/local/AnalysisSoftwareEIC/treeAnalysis/treeProcessing/CENTRALSIM_PYTHIA_TTL8_EoP_MB/output_TMSTUD.root

root -x -l -b -q 'eoverpstudies.C("'$input'","pdf",true)'
