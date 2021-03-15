
fileLBLACLGADMB=/home/fbock/eic/simulationOutput/Modular/output_ALLSILICON-FTTLS3LC-ETTL-CTTL-ACLGAD-TREXTOUT_e10p250MB/G4EICDetector_eventtree.root
fileLBLACLGADpTHard=/home/fbock/eic/simulationOutput/Modular/output_ALLSILICON-FTTLS3LC-ETTL-CTTL-ACLGAD-TREXTOUT_e10p250pTHard5/G4EICDetector_eventtree.root
fileLBLMB=/home/fbock/eic/simulationOutput/Modular/output_ALLSILICON-TREXTOUT_e10p250MB/G4EICDetector_eventtree.root
fileLBLpTHard=/home/fbock/eic/simulationOutput/Modular/output_ALLSILICON-TREXTOUT_e10p250pTHard5/G4EICDetector_eventtree.root
fileSinglePionOnlyCalo=/home/fbock/eic/simulationOutput/output_CALOSTANDALONE_SimplePion/G4EICDetector_eventtree_all.root
fileLBLLGADMB=/home/fbock/eic/simulationOutput/Modular/output_ALLSILICON-FTTLS3LC-ETTL-CTTL-TREXTOUT_e10p250MB/G4EICDetector_eventtree.root
fileLBLLGADpTHard=/home/fbock/eic/simulationOutput/Modular/output_ALLSILICON-FTTLS3LC-ETTL-CTTL-TREXTOUT_e10p250pTHard5/G4EICDetector_eventtree.root
fileLBLLGADpTHard2x=/home/fbock/eic/simulationOutput/Modular/output_ALLSILICON-FTTLS3LC-ETTL-CTTL-TREXTOUT-HC2xEC2x_e10p250pTHard5/G4EICDetector_eventtree.root

fileLBLLGAD5100MB=/home/fbock/eic/simulationOutput/Modular/output_ALLSILICON-FTTLS3LC-ETTL-CTTL-TREXTOUT_e5p100MB/G4EICDetector_eventtree.root
fileLBLLGAD5100pTHard=/home/fbock/eic/simulationOutput/Modular/output_ALLSILICON-FTTLS3LC-ETTL-CTTL-TREXTOUT_e5p100pTHard5/G4EICDetector_eventtree.root
fileLBLLGAD18275MB=/home/fbock/eic/simulationOutput/Modular/output_ALLSILICON-FTTLS3LC-ETTL-CTTL-TREXTOUT_e18p275MB/G4EICDetector_eventtree.root
fileLBLLGAD18275pTHard=/home/fbock/eic/simulationOutput/Modular/output_ALLSILICON-FTTLS3LC-ETTL-CTTL-TREXTOUT_e18p275pTHard5/G4EICDetector_eventtree.root

file3TLBLACLGADMB=/home/fbock/eic/simulationOutput/ModularBeast/output_ALLSILICON-FTTLS3LC-ETTL-CTTL-ACLGAD-TREXTOUT_e10p250MB/G4EICDetector_eventtree.root
file3TLBLACLGADpTHard=/home/fbock/eic/simulationOutput/ModularBeast/output_ALLSILICON-FTTLS3LC-ETTL-CTTL-ACLGAD-TREXTOUT_e10p250pTHard5/G4EICDetector_eventtree.root
file3TLBLMB=/home/fbock/eic/simulationOutput/ModularBeast/output_ALLSILICON-TREXTOUT_e10p250MB/G4EICDetector_eventtree.root
file3TLBLpTHard=/home/fbock/eic/simulationOutput/ModularBeast/output_ALLSILICON-TREXTOUT_e10p250pTHard5/G4EICDetector_eventtree.root
file3TLBLLGADMB=/home/fbock/eic/simulationOutput/ModularBeast/output_ALLSILICON-FTTLS3LC-ETTL-CTTL-TREXTOUT_e10p250MB/G4EICDetector_eventtree.root
file3TLBLLGADpTHard=/home/fbock/eic/simulationOutput/ModularBeast/output_ALLSILICON-FTTLS3LC-ETTL-CTTL-TREXTOUT_e10p250pTHard5/G4EICDetector_eventtree.root
file3TLBLLGADpTHard2x=/home/fbock/eic/simulationOutput/ModularBeast/output_ALLSILICON-FTTLS3LC-ETTL-CTTL-TREXTOUT-HC2xEC2x_e10p250pTHard5/G4EICDetector_eventtree.root

file3TLBLLGAD5100MB=/home/fbock/eic/simulationOutput/ModularBeast/output_ALLSILICON-FTTLS3LC-ETTL-CTTL-TREXTOUT_e5p100MB/G4EICDetector_eventtree.root
file3TLBLLGAD5100pTHard=/home/fbock/eic/simulationOutput/ModularBeast/output_ALLSILICON-FTTLS3LC-ETTL-CTTL-TREXTOUT_e5p100pTHard5/G4EICDetector_eventtree.root
file3TLBLLGAD18275MB=/home/fbock/eic/simulationOutput/ModularBeast/output_ALLSILICON-FTTLS3LC-ETTL-CTTL-TREXTOUT_e18p275MB/G4EICDetector_eventtree.root
file3TLBLLGAD18275pTHard=/home/fbock/eic/simulationOutput/ModularBeast/output_ALLSILICON-FTTLS3LC-ETTL-CTTL-TREXTOUT_e18p275pTHard5/G4EICDetector_eventtree.root

# # # single pion
# root -x -l -b -q 'treeProcessing.C("'$fileSinglePionOnlyCalo'","CaloOnly_SinglePi",true,true,true,true,true,true,true,true,true,true,true,true,false,false,false-1,0)'
# # 
# # LBL AC-LGAD 
# root -x -l -b -q 'treeProcessing.C("'$fileLBLACLGADpTHard'","LBLACLGAD_pTHard",true,true,true,true,true,true,true,true,true,true,true,true,false,true,true,-1,0)'
# root -x -l -b -q 'treeProcessing.C("'$fileLBLACLGADMB'","LBLACLGAD_MB",true,true,true,true,true,true,true,true,true,true,true,true,false,true,true,-1,0)'
LBL 
# root -x -l -b -q 'treeProcessing.C("'$fileLBLpTHard'","LBL_pTHard",true,true,true,true,true,true,true,true,true,true,true,true,false,true,true,-1,0)'
# root -x -l -b -q 'treeProcessing.C("'$fileLBLMB'","LBL_MB",true,true,true,true,true,true,true,true,true,true,true,true,false,true,true,-1,0)'
# # LBL LGAD
# root -x -l -b -q 'treeProcessing.C("'$fileLBLLGADMB'","LBLLGAD_MB",true,true,true,true,true,true,true,true,true,true,true,true,false,true,true,-1,0)'
root -x -l -b -q 'treeProcessing.C("'$fileLBLLGADpTHard'","LBLLGAD_pTHard",true,true,true,true,true,true,true,true,true,true,true,true,false,true,true,-1,0)'
# # LBL LGAD double granul 
# root -x -l -b -q 'treeProcessing.C("'$fileLBLLGADpTHard2x'","LBLLGAD2x_pTHard",true,true,true,true,true,true,true,true,true,true,true,true,false,true,true,-1,0)'

# LBL LGAD other cms energies
# root -x -l -b -q 'treeProcessing.C("'$fileLBLLGAD5100MB'","LBLLGAD_e5p100_MB",true,true,true,true,true,true,true,true,true,true,true,true,false,true,true,-1,0)'
# root -x -l -b -q 'treeProcessing.C("'$fileLBLLGAD5100pTHard'","LBLLGAD_e5p100_pTHard",true,true,true,true,true,true,true,true,true,true,true,true,false,true,true,-1,0)'
# root -x -l -b -q 'treeProcessing.C("'$fileLBLLGAD18275MB'","LBLLGAD_e18p275_MB",true,true,true,true,true,true,true,true,true,true,true,true,false,true,true,-1,0)'
# root -x -l -b -q 'treeProcessing.C("'$fileLBLLGAD18275pTHard'","LBLLGAD_e18p275_pTHard",true,true,true,true,true,true,true,true,true,true,true,true,false,true,true,-1,0)'
# 
# LBL AC-LGAD 3T
# root -x -l -b -q 'treeProcessing.C("'$file3TLBLACLGADpTHard'","LBLACLGAD_3T_pTHard",true,true,true,true,true,true,true,true,true,true,true,true,false,true,true,-1,0)'
# root -x -l -b -q 'treeProcessing.C("'$file3TLBLACLGADMB'","LBLACLGAD_3T_MB",true,true,true,true,true,true,true,true,true,true,true,true,false,true,true,-1,0)'
# LBL 3T
# root -x -l -b -q 'treeProcessing.C("'$file3TLBLpTHard'","LBL_3T_pTHard",true,true,true,true,true,true,true,true,true,true,true,true,false,true,true,-1,0)'
# root -x -l -b -q 'treeProcessing.C("'$file3TLBLMB'","LBL_3T_MB",true,true,true,true,true,true,true,true,true,true,true,true,false,true,true,-1,0)'
# LBL LGAD 3T
# root -x -l -b -q 'treeProcessing.C("'$file3TLBLLGADMB'","LBLLGAD_3T_MB",true,true,true,true,true,true,true,true,true,true,true,true,false,true,true,-1,0)'
# root -x -l -b -q 'treeProcessing.C("'$file3TLBLLGADpTHard'","LBLLGAD_3T_pTHard",true,true,true,true,true,true,true,true,true,true,true,true,false,true,true,-1,0)'
# LBL LGAD double granul 3T
# root -x -l -b -q 'treeProcessing.C("'$file3TLBLLGADpTHard2x'","LBLLGAD2x_3T_pTHard",true,true,true,true,true,true,true,true,true,true,true,true,false,true,true,-1,0)'
# 
# # LBL LGAD 3T other cms energies
# root -x -l -b -q 'treeProcessing.C("'$file3TLBLLGAD5100MB'","LBLLGAD_3T_e5p100_MB",true,true,true,true,true,true,true,true,true,true,true,true,false,true,true,-1,0)'
# root -x -l -b -q 'treeProcessing.C("'$file3TLBLLGAD5100pTHard'","LBLLGAD_3T_e5p100_pTHard",true,true,true,true,true,true,true,true,true,true,true,true,false,true,true,-1,0)'
# root -x -l -b -q 'treeProcessing.C("'$file3TLBLLGAD18275MB'","LBLLGAD_3T_e18p275_MB",true,true,true,true,true,true,true,true,true,true,true,true,false,true,true,-1,0)'
# root -x -l -b -q 'treeProcessing.C("'$file3TLBLLGAD18275pTHard'","LBLLGAD_3T_e18p275_pTHard",true,true,true,true,true,true,true,true,true,true,true,true,false,true,true,-1,0)'
