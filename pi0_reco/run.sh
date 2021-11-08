#input=/media/nschmidt/local/AnalysisSoftwareEIC/treeAnalysis/treeProcessing/Pi0SIM_TTL8/output_pi0.root
input=/media/nschmidt/local/AnalysisSoftwareEIC/treeAnalysis/treeProcessing/CENTRALSIM_PYTHIA_TTL8/output_pi0.root
calo=FEMC
algo=MA
root -x -l -b -q 'pi0resolution.C("'$input'","'$algo'","pdf",true)'
calo=EEMC
#root -x -l -b -q 'pi0resolution.C("'$input'","'$calo'","'$algo'","pdf",true)'
calo=BECAL
#root -x -l -b -q 'pi0resolution.C("'$input'","'$calo'","'$algo'","pdf",true)'

