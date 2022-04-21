#input=/media/nschmidt/local/AnalysisSoftwareEIC/treeAnalysis/treeProcessing/Pi0SIM_TTL7/output_pi0.root
input=/media/nschmidt/local/AnalysisSoftwareEIC/treeAnalysis/treeProcessing/PI0SIM_ECCEGEOM/output_pi0.root
#input=/media/nschmidt/local/AnalysisSoftwareEIC/treeAnalysis/treeProcessing/CENTRALSIM_PYTHIA_TTL8/output_pi0.root
algo=MA
root -x -l -b -q 'pi0resolution.C("'$input'","'$algo'","pdf",true)'
