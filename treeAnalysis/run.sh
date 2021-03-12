
fileLBLACLGADMB=/home/fbock/eic/simulationOutput/output_EVTTREE-ALLSILICON-FTTLS3LC-ETTL-CTTL-ACLGAD-TREXTOUT_epMB/G4EICDetector_eventtree.root
fileLBLACLGADpTHard=/home/fbock/eic/simulationOutput/output_EVTTREE-ALLSILICON-FTTLS3LC-ETTL-CTTL-ACLGAD-TREXTOUT_pTHard5/G4EICDetector_eventtree.root
fileLBLMB=/home/fbock/eic/simulationOutput/output_EVTTREE-ALLSILICON-TREXTOUT_epMB/G4EICDetector_eventtree.root
fileLBLpTHard=/home/fbock/eic/simulationOutput/output_EVTTREE-ALLSILICON-TREXTOUT_pTHard5/G4EICDetector_eventtree.root
fileSinglePionOnlyCalo=/home/fbock/eic/simulationOutput/output_CALOSTANDALONE_SimplePion/G4EICDetector_eventtree_all.root
     
# root -x -l -b -q 'treeProcessing.C("'$fileLBLACLGADpTHard'","LBLACLGAD_pTHard",true,true,true,true,true,true,true,true,true,true,true,true,false,-1,0)'
# root -x -l -b -q 'treeProcessing.C("'$fileLBLACLGADMB'","LBLACLGAD_MB",true,true,true,true,true,true,true,true,true,true,true,true,false,-1,0)'
# root -x -l -b -q 'treeProcessing.C("'$fileLBLpTHard'","LBL_pTHard",true,true,true,true,true,true,true,true,true,true,true,true,false,-1,0)'
# root -x -l -b -q 'treeProcessing.C("'$fileLBLMB'","LBL_MB",true,true,true,true,true,true,true,true,true,true,true,true,false,-1,0)'
root -x -l -b -q 'treeProcessing.C("'$fileSinglePionOnlyCalo'","CaloOnly_SinglePi",true,true,true,true,true,true,true,true,true,true,true,true,false,-1,0)'
