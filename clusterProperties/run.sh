# steering script to get effi extraction and comparison running
if [ $1 == "effiZF" ] ; then
  inputFileTR="/home/fbock/EIC/Analysis/CaloNIM/treeAnalysis/treeProcessing/ZeroField-ElectronAndPion/output_TRKEFF.root"
  inputFileCS="/home/fbock/EIC/Analysis/CaloNIM/treeAnalysis/treeProcessing/ZeroField-ElectronAndPion/output_CS.root"

  root -x -q -l -b 'clustereffi.C("'$inputFileTR'","'$inputFileCS'","LFHCAL","pdf","SinglePart","",0,0)'
  root -x -q -l -b 'clustereffi.C("'$inputFileTR'","'$inputFileCS'","HCALOUT","pdf","SinglePart","",0,0)'
  root -x -q -l -b 'clustereffi.C("'$inputFileTR'","'$inputFileCS'","HCALIN","pdf","SinglePart","",0,0)'
  root -x -q -l -b 'clustereffi.C("'$inputFileTR'","'$inputFileCS'","EHCAL","pdf","SinglePart","",0,0)'
  root -x -q -l -b 'clustereffi.C("'$inputFileTR'","'$inputFileCS'","FEMC","pdf","SinglePart","",0,0)'
  root -x -q -l -b 'clustereffi.C("'$inputFileTR'","'$inputFileCS'","BECAL","pdf","SinglePart","",0,0)'
  root -x -q -l -b 'clustereffi.C("'$inputFileTR'","'$inputFileCS'","EEMC","pdf","SinglePart","",0,0)'
elif [ $1 == "effiProp" ] ; then
  # proposal files
  inputFileTR="/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/SingleParticleEEMCPCarbonEOP/output_TRKEFF.root"
  inputFileCS="/home/fbock/EIC/Analysis/2ndSimulationCampaign/treeAnalysis/treeProcessing/SingleParticleEEMCPCarbonEOP/output_CS.root"

  root -x -q -l -b 'clustereffi.C("'$inputFileTR'","'$inputFileCS'","LFHCAL","pdf","SinglePart","",0,0)'
  root -x -q -l -b 'clustereffi.C("'$inputFileTR'","'$inputFileCS'","HCALOUT","pdf","SinglePart","",0,0)'
  root -x -q -l -b 'clustereffi.C("'$inputFileTR'","'$inputFileCS'","HCALIN","pdf","SinglePart","",0,0)'
  root -x -q -l -b 'clustereffi.C("'$inputFileTR'","'$inputFileCS'","EHCAL","pdf","SinglePart","",0,0)'
  root -x -q -l -b 'clustereffi.C("'$inputFileTR'","'$inputFileCS'","FEMC","pdf","SinglePart","",0,0)'
  root -x -q -l -b 'clustereffi.C("'$inputFileTR'","'$inputFileCS'","BECAL","pdf","SinglePart","",0,0)'
  root -x -q -l -b 'clustereffi.C("'$inputFileTR'","'$inputFileCS'","EEMC","pdf","SinglePart","",0,0)'

elif [ $1 == "effiNF" ] ; then #nominal field 
  # proposal files
  inputFileTR="/home/fbock/EIC/Analysis/CaloNIM/treeAnalysis/treeProcessing/Default-ElectronAndPion/output_TRKEFF.root"
  inputFileCS="/home/fbock/EIC/Analysis/CaloNIM/treeAnalysis/treeProcessing/Default-ElectronAndPion/output_CS.root"

  root -x -q -l -b 'clustereffi.C("'$inputFileTR'","'$inputFileCS'","LFHCAL","pdf","SinglePart","",0,0)'
  root -x -q -l -b 'clustereffi.C("'$inputFileTR'","'$inputFileCS'","HCALOUT","pdf","SinglePart","",0,0)'
  root -x -q -l -b 'clustereffi.C("'$inputFileTR'","'$inputFileCS'","HCALIN","pdf","SinglePart","",0,0)'
  root -x -q -l -b 'clustereffi.C("'$inputFileTR'","'$inputFileCS'","EHCAL","pdf","SinglePart","",0,0)'
  root -x -q -l -b 'clustereffi.C("'$inputFileTR'","'$inputFileCS'","FEMC","pdf","SinglePart","",0,0)'
  root -x -q -l -b 'clustereffi.C("'$inputFileTR'","'$inputFileCS'","BECAL","pdf","SinglePart","",0,0)'
  root -x -q -l -b 'clustereffi.C("'$inputFileTR'","'$inputFileCS'","EEMC","pdf","SinglePart","",0,0)'

elif [ $1 == "effiCompZF" ]; then
  echo -e "/home/fbock/EIC/Analysis/CaloNIM/treeAnalysis/treeProcessing/ZeroField-ElectronAndPion/output_CS.root EEMC 0 EEMC  0" > compareECalsHCals.txt
  echo -e "/home/fbock/EIC/Analysis/CaloNIM/treeAnalysis/treeProcessing/ZeroField-ElectronAndPion/output_CS.root BECAL   1   BEMC 1" >> compareECalsHCals.txt
  echo -e "/home/fbock/EIC/Analysis/CaloNIM/treeAnalysis/treeProcessing/ZeroField-ElectronAndPion/output_CS.root FEMC    2   FEMC 2" >> compareECalsHCals.txt
  echo -e "/home/fbock/EIC/Analysis/CaloNIM/treeAnalysis/treeProcessing/ZeroField-ElectronAndPion/output_CS.root HCALOUT    1   OHCAL 5" >> compareECalsHCals.txt
  echo -e "/home/fbock/EIC/Analysis/CaloNIM/treeAnalysis/treeProcessing/ZeroField-ElectronAndPion/output_CS.root LFHCAL   2   LFHCAL  6" >> compareECalsHCals.txt  
  root -x -q -l -b 'compare_clustereffi2.C+("compareECalsHCals.txt","E+HCAL-ZeroField","pdf")'
elif [ $1 == "effiCompNF" ]; then
  echo -e "/home/fbock/EIC/Analysis/CaloNIM/treeAnalysis/treeProcessing/Default-ElectronAndPion/output_CS.root EEMC 0 EEMC  0" > compareECalsHCals.txt
  echo -e "/home/fbock/EIC/Analysis/CaloNIM/treeAnalysis/treeProcessing/Default-ElectronAndPion/output_CS.root BECAL   1   BEMC 1" >> compareECalsHCals.txt
  echo -e "/home/fbock/EIC/Analysis/CaloNIM/treeAnalysis/treeProcessing/Default-ElectronAndPion/output_CS.root FEMC    2   FEMC 2" >> compareECalsHCals.txt
  echo -e "/home/fbock/EIC/Analysis/CaloNIM/treeAnalysis/treeProcessing/Default-ElectronAndPion/output_CS.root HCALOUT    1   OHCAL 5" >> compareECalsHCals.txt
  echo -e "/home/fbock/EIC/Analysis/CaloNIM/treeAnalysis/treeProcessing/Default-ElectronAndPion/output_CS.root LFHCAL   2   LFHCAL  6" >> compareECalsHCals.txt  
  root -x -q -l -b 'compare_clustereffi2.C+("compareECalsHCals.txt","E+HCAL-Default","pdf")'

fi
