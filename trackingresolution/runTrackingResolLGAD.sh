fileLBL=rootfiles/histosTrackingPerformance_LBL.root
fileLBLLGAD=rootfiles/histosTrackingPerformance_LBLandFullTTL.root
fileLBLACLGAD=rootfiles/histosTrackingPerformance_LBLandFullTTLACLGAD.root
fileLBL2LTTL=rootfiles/histosTrackingPerformance_LBLandFTTLS2LC-ETTL-CTTL.root
fileLBLTTLC=rootfiles/histosTrackingPerformance_LBLandTTLCoarse.root
fileLBL1LTTL=rootfiles/histosTrackingPerformance_LBLandFTTLSE1LC-ETTLSE1-CTTLSE1.root
fileLBL2LBETTL=rootfiles/histosTrackingPerformance_LBLandFTTLSE2LC-ETTL-CTTLSE1.root

fileLANL=rootfiles/histosTrackingPerformance_LANL.root
fileLANLLGAD=rootfiles/histosTrackingPerformance_LANLandFullTTL.root
fileLANLACLGAD=rootfiles/histosTrackingPerformance_LANLandFullTTLACLGAD.root
fileLANLTTLC=rootfiles/histosTrackingPerformance_LANLandTTLCoarse.root

fileLBLMB=rootfiles/histosTrackingPerformance_LBL-MB.root
fileLBLLGADMB=rootfiles/histosTrackingPerformance_LBLandFullTTL-MB.root
  

if [ $1 = "plots" ]; then
  fileLBLLGAD=rootfiles/histosTrackingResol_LBLwithLGAD.root
  fileLBL=rootfiles/histosTrackingResol_defaultLBL.root
  fileFST=rootfiles/histosTrackingResol_defaultFST.root
  fileFSTLGAD=rootfiles/histosTrackingResol_FST-LGAD.root
        
  root -b -x -q -l trackingreso_flatPt.C\(\"$fileLBL\"\,\"pdf\"\,\"defaultLBL-Hist\"\,kFALSE\)
  root -b -x -q -l trackingreso_flatPt.C\(\"$fileLBL\"\,\"pdf\"\,\"defaultLBL-Fit\"\)

  root -b -x -q -l trackingreso_flatPt.C\(\"$fileLBLLGAD\"\,\"pdf\"\,\"LBLwithLGAD-Hist\"\,kFALSE\)
  root -b -x -q -l trackingreso_flatPt.C\(\"$fileLBLLGAD\"\,\"pdf\"\,\"LBLwithLGAD-Fit\"\)

  root -b -x -q -l trackingreso_flatPt.C\(\"$fileFST\"\,\"pdf\"\,\"defaultFST-Hist\"\,kFALSE\)
  root -b -x -q -l trackingreso_flatPt.C\(\"$fileFST\"\,\"pdf\"\,\"defaultFST-Fit\"\)

  root -b -x -q -l trackingreso_flatPt.C\(\"$fileFSTLGAD\"\,\"pdf\"\,\"FSTwithLGAD-Hist\"\,kFALSE\)
  root -b -x -q -l trackingreso_flatPt.C\(\"$fileFSTLGAD\"\,\"pdf\"\,\"FSTwithLGAD-Fit\"\)
elif [ $1 = "plotsPyt" ]; then
#   root -b -x -q -l trackingreso_Pythia.C\(\"$fileLBL\"\,\"pdf\"\,\"defaultLBL-Fit\"\)
#   root -b -x -q -l trackingreso_Pythia.C\(\"$fileLBLLGAD\"\,\"pdf\"\,\"LBLwithLGAD-Fit\"\)
#   root -b -x -q -l trackingreso_Pythia.C\(\"$fileLBLACLGAD\"\,\"pdf\"\,\"LBLwithACLGAD-Fit\"\)
#   root -b -x -q -l trackingreso_Pythia.C\(\"$fileLBL2LTTL\"\,\"pdf\"\,\"LBLwithFTTLS2LC-ETTL-CTTL-Fit\"\)
  #   root -b -x -q -l trackingreso_Pythia.C\(\"$fileLBLTTLC\"\,\"pdf\"\,\"LBLwithFTTLS3LVC-ETTLLC-CTTLLC-Fit\"\)
  #   root -b -x -q -l trackingreso_Pythia.C\(\"$fileLBL1LTTL\"\,\"pdf\"\,\"LBLwithFTTLSE1LC-ETTLSE1-CTTLSE1-Fit\"\)
  #   root -b -x -q -l trackingreso_Pythia.C\(\"$fileLBL2LBETTL\"\,\"pdf\"\,\"LBLwithFTTLSE2LC-ETTL-CTTLSE1-Fit\"\)

  root -b -x -q -l trackingreso_Pythia.C\(\"$fileLBL\"\,\"pdf\"\,\"defaultLBL-Hist-pTHard5GeV\"\,kFALSE\)
  root -b -x -q -l trackingreso_Pythia.C\(\"$fileLBLLGAD\"\,\"pdf\"\,\"LBLwithLGAD-Hist-pTHard5GeV\"\,kFALSE\)
  root -b -x -q -l trackingreso_Pythia.C\(\"$fileLBLACLGAD\"\,\"pdf\"\,\"LBLwithACLGAD-Hist-pTHard5GeV\"\,kFALSE\)
  root -b -x -q -l trackingreso_Pythia.C\(\"$fileLBL2LTTL\"\,\"pdf\"\,\"LBLwithFTTLS2LC-ETTL-CTTL-Hist-pTHard5GeV\"\,kFALSE\)
  root -b -x -q -l trackingreso_Pythia.C\(\"$fileLBLTTLC\"\,\"pdf\"\,\"LBLwithFTTLS3LVC-ETTLLC-CTTLLC-Hist-pTHard5GeV\"\,kFALSE\)
  root -b -x -q -l trackingreso_Pythia.C\(\"$fileLBL1LTTL\"\,\"pdf\"\,\"LBLwithFTTLSE1LC-ETTLSE1-CTTLSE1-Hist-pTHard5GeV\"\,kFALSE\)
  root -b -x -q -l trackingreso_Pythia.C\(\"$fileLBL2LBETTL\"\,\"pdf\"\,\"LBLwithFTTLSE2LC-ETTL-CTTLSE1-Hist-pTHard5GeV\"\,kFALSE\)
  root -b -x -q -l trackingreso_Pythia.C\(\"$fileLANL\"\,\"pdf\"\,\"defaultLANL-Hist-pTHard5GeV\"\,kFALSE\)
  root -b -x -q -l trackingreso_Pythia.C\(\"$fileLANLLGAD\"\,\"pdf\"\,\"LANLwithLGAD-Hist-pTHard5GeV\"\,kFALSE\)
  root -b -x -q -l trackingreso_Pythia.C\(\"$fileLANLACLGAD\"\,\"pdf\"\,\"LANLwithACLGAD-Hist-pTHard5GeV\"\,kFALSE\)
  root -b -x -q -l trackingreso_Pythia.C\(\"$fileLANLTTLC\"\,\"pdf\"\,\"LANLwithFTTLS3LVC-ETTLLC-CTTLLC-Hist-pTHard5GeV\"\,kFALSE\)
  root -b -x -q -l trackingreso_Pythia.C\(\"$fileLBLMB\"\,\"pdf\"\,\"defaultLBL-Hist-MB\"\,kFALSE\)
  root -b -x -q -l trackingreso_Pythia.C\(\"$fileLBLLGADMB\"\,\"pdf\"\,\"LBLwithLGAD-Hist-MB\"\,kFALSE\)
  
  root -b -x -q -l trackingreso_Pythia.C\(\"$fileLBL\"\,\"pdf\"\,\"defaultLBL-Hist-pTHard5GeV\"\,kFALSE\,1\)
  root -b -x -q -l trackingreso_Pythia.C\(\"$fileLBLLGAD\"\,\"pdf\"\,\"LBLwithLGAD-Hist-pTHard5GeV\"\,kFALSE\,1\)
  root -b -x -q -l trackingreso_Pythia.C\(\"$fileLBLACLGAD\"\,\"pdf\"\,\"LBLwithACLGAD-Hist-pTHard5GeV\"\,kFALSE\,1\)
  root -b -x -q -l trackingreso_Pythia.C\(\"$fileLBL2LTTL\"\,\"pdf\"\,\"LBLwithFTTLS2LC-ETTL-CTTL-Hist-pTHard5GeV\"\,kFALSE\,1\)
  root -b -x -q -l trackingreso_Pythia.C\(\"$fileLBLTTLC\"\,\"pdf\"\,\"LBLwithFTTLS3LVC-ETTLLC-CTTLLC-Hist-pTHard5GeV\"\,kFALSE\,1\)
  root -b -x -q -l trackingreso_Pythia.C\(\"$fileLBL1LTTL\"\,\"pdf\"\,\"LBLwithFTTLSE1LC-ETTLSE1-CTTLSE1-Hist-pTHard5GeV\"\,kFALSE\,1\)
  root -b -x -q -l trackingreso_Pythia.C\(\"$fileLBL2LBETTL\"\,\"pdf\"\,\"LBLwithFTTLSE2LC-ETTL-CTTLSE1-Hist-pTHard5GeV\"\,kFALSE\,1\)
  root -b -x -q -l trackingreso_Pythia.C\(\"$fileLANL\"\,\"pdf\"\,\"defaultLANL-Hist-pTHard5GeV\"\,kFALSE\,1\)
  root -b -x -q -l trackingreso_Pythia.C\(\"$fileLANLLGAD\"\,\"pdf\"\,\"LANLwithLGAD-Hist-pTHard5GeV\"\,kFALSE\,1\)
  root -b -x -q -l trackingreso_Pythia.C\(\"$fileLANLACLGAD\"\,\"pdf\"\,\"LANLwithACLGAD-Hist-pTHard5GeV\"\,kFALSE\,1\)
  root -b -x -q -l trackingreso_Pythia.C\(\"$fileLANLTTLC\"\,\"pdf\"\,\"LANLwithFTTLS3LVC-ETTLLC-CTTLLC-Hist-pTHard5GeV\"\,kFALSE\,1\)
  root -b -x -q -l trackingreso_Pythia.C\(\"$fileLBLMB\"\,\"pdf\"\,\"defaultLBL-Hist-MB\"\,kFALSE\,1\)
  root -b -x -q -l trackingreso_Pythia.C\(\"$fileLBLLGADMB\"\,\"pdf\"\,\"LBLwithLGAD-Hist-MB\"\,kFALSE\,1\)

  root -b -x -q -l trackingreso_Pythia.C\(\"$fileLBL\"\,\"pdf\"\,\"defaultLBL-Hist-pTHard5GeV\"\,kFALSE\,2\)
  root -b -x -q -l trackingreso_Pythia.C\(\"$fileLBLLGAD\"\,\"pdf\"\,\"LBLwithLGAD-Hist-pTHard5GeV\"\,kFALSE\,2\)
  root -b -x -q -l trackingreso_Pythia.C\(\"$fileLBLACLGAD\"\,\"pdf\"\,\"LBLwithACLGAD-Hist-pTHard5GeV\"\,kFALSE\,2\)
  root -b -x -q -l trackingreso_Pythia.C\(\"$fileLBL2LTTL\"\,\"pdf\"\,\"LBLwithFTTLS2LC-ETTL-CTTL-Hist-pTHard5GeV\"\,kFALSE\,2\)
  root -b -x -q -l trackingreso_Pythia.C\(\"$fileLBLTTLC\"\,\"pdf\"\,\"LBLwithFTTLS3LVC-ETTLLC-CTTLLC-Hist-pTHard5GeV\"\,kFALSE\,2\)
  root -b -x -q -l trackingreso_Pythia.C\(\"$fileLBL1LTTL\"\,\"pdf\"\,\"LBLwithFTTLSE1LC-ETTLSE1-CTTLSE1-Hist-pTHard5GeV\"\,kFALSE\,2\)
  root -b -x -q -l trackingreso_Pythia.C\(\"$fileLBL2LBETTL\"\,\"pdf\"\,\"LBLwithFTTLSE2LC-ETTL-CTTLSE1-Hist-pTHard5GeV\"\,kFALSE\,2\)
  root -b -x -q -l trackingreso_Pythia.C\(\"$fileLANL\"\,\"pdf\"\,\"defaultLANL-Hist-pTHard5GeV\"\,kFALSE\,2\)
  root -b -x -q -l trackingreso_Pythia.C\(\"$fileLANLLGAD\"\,\"pdf\"\,\"LANLwithLGAD-Hist-pTHard5GeV\"\,kFALSE\,2\)
  root -b -x -q -l trackingreso_Pythia.C\(\"$fileLANLACLGAD\"\,\"pdf\"\,\"LANLwithACLGAD-Hist-pTHard5GeV\"\,kFALSE\,2\)
  root -b -x -q -l trackingreso_Pythia.C\(\"$fileLANLTTLC\"\,\"pdf\"\,\"LANLwithFTTLS3LVC-ETTLLC-CTTLLC-Hist-pTHard5GeV\"\,kFALSE\,2\)
  root -b -x -q -l trackingreso_Pythia.C\(\"$fileLBLMB\"\,\"pdf\"\,\"defaultLBL-Hist-MB\"\,kFALSE\,2\)
  root -b -x -q -l trackingreso_Pythia.C\(\"$fileLBLLGADMB\"\,\"pdf\"\,\"LBLwithLGAD-Hist-MB\"\,kFALSE\,2\)
  
# elif [ $1 = "plotsPID" ]; then
  
  root -b -x -q -l pidreso_Pythia.C\(\"$fileLBLLGAD\"\,\"pdf\"\,\"LBLwithLGAD-Hist-pTHard5GeV\"\,kFALSE\)
  root -b -x -q -l pidreso_Pythia.C\(\"$fileLBLACLGAD\"\,\"pdf\"\,\"LBLwithACLGAD-Hist-pTHard5GeV\"\,kFALSE\)
  root -b -x -q -l pidreso_Pythia.C\(\"$fileLBL2LTTL\"\,\"pdf\"\,\"LBLwithFTTLS2LC-ETTL-CTTL-Hist-pTHard5GeV\"\,kFALSE\)
  root -b -x -q -l pidreso_Pythia.C\(\"$fileLBLTTLC\"\,\"pdf\"\,\"LBLwithFTTLS3LVC-ETTLLC-CTTLLC-Hist-pTHard5GeV\"\,kFALSE\)
  root -b -x -q -l pidreso_Pythia.C\(\"$fileLBL1LTTL\"\,\"pdf\"\,\"LBLwithFTTLSE1LC-ETTLSE1-CTTLSE1-Hist-pTHard5GeV\"\,kFALSE\)
  root -b -x -q -l pidreso_Pythia.C\(\"$fileLBL2LBETTL\"\,\"pdf\"\,\"LBLwithFTTLSE2LC-ETTL-CTTLSE1-Hist-pTHard5GeV\"\,kFALSE\)
  root -b -x -q -l pidreso_Pythia.C\(\"$fileLANLLGAD\"\,\"pdf\"\,\"LANLwithLGAD-Hist-pTHard5GeV\"\,kFALSE\)
  root -b -x -q -l pidreso_Pythia.C\(\"$fileLANLACLGAD\"\,\"pdf\"\,\"LANLwithACLGAD-Hist-pTHard5GeV\"\,kFALSE\)
  root -b -x -q -l pidreso_Pythia.C\(\"$fileLANLTTLC\"\,\"pdf\"\,\"LANLwithFTTLS3LVC-ETTLLC-CTTLLC-Hist-pTHard5GeV\"\,kFALSE\)
  root -b -x -q -l pidreso_Pythia.C\(\"$fileLBLLGADMB\"\,\"pdf\"\,\"LBLwithLGAD-Hist-MB\"\,kFALSE\)

  root -b -x -q -l pidreso_Pythia.C\(\"$fileLBLLGAD\"\,\"pdf\"\,\"LBLwithLGAD-Hist-pTHard5GeV\"\,kFALSE\,2\)
  root -b -x -q -l pidreso_Pythia.C\(\"$fileLBLACLGAD\"\,\"pdf\"\,\"LBLwithACLGAD-Hist-pTHard5GeV\"\,kFALSE\,2\)
  root -b -x -q -l pidreso_Pythia.C\(\"$fileLBL2LTTL\"\,\"pdf\"\,\"LBLwithFTTLS2LC-ETTL-CTTL-Hist-pTHard5GeV\"\,kFALSE\,2\)
  root -b -x -q -l pidreso_Pythia.C\(\"$fileLBLTTLC\"\,\"pdf\"\,\"LBLwithFTTLS3LVC-ETTLLC-CTTLLC-Hist-pTHard5GeV\"\,kFALSE\,2\)
  root -b -x -q -l pidreso_Pythia.C\(\"$fileLBL1LTTL\"\,\"pdf\"\,\"LBLwithFTTLSE1LC-ETTLSE1-CTTLSE1-Hist-pTHard5GeV\"\,kFALSE\,2\)
  root -b -x -q -l pidreso_Pythia.C\(\"$fileLBL2LBETTL\"\,\"pdf\"\,\"LBLwithFTTLSE2LC-ETTL-CTTLSE1-Hist-pTHard5GeV\"\,kFALSE\,2\)
  root -b -x -q -l pidreso_Pythia.C\(\"$fileLANLLGAD\"\,\"pdf\"\,\"LANLwithLGAD-Hist-pTHard5GeV\"\,kFALSE\,2\)
  root -b -x -q -l pidreso_Pythia.C\(\"$fileLANLACLGAD\"\,\"pdf\"\,\"LANLwithACLGAD-Hist-pTHard5GeV\"\,kFALSE\,2\)
  root -b -x -q -l pidreso_Pythia.C\(\"$fileLANLTTLC\"\,\"pdf\"\,\"LANLwithFTTLS3LVC-ETTLLC-CTTLLC-Hist-pTHard5GeV\"\,kFALSE\,2\)
  root -b -x -q -l pidreso_Pythia.C\(\"$fileLBLLGADMB\"\,\"pdf\"\,\"LBLwithLGAD-Hist-MB\"\,kFALSE\,2\)

elif [ $1 = "compare" ]; then
  root -q -x -l -b 'compare_trackingreso_Pythia.C("filesCompareReso_granularity.txt","Granularity-pTHard5GeV","pdf",kFALSE)'
  root -q -x -l -b 'compare_trackingreso_Pythia.C("filesCompareReso_Setups.txt","Setups-pTHard5GeV","pdf",kFALSE)'
  root -q -x -l -b 'compare_trackingreso_Pythia.C("filesCompareReso_SetupsRed.txt","SetupsRed-pTHard5GeV","pdf",kFALSE)'
  root -q -x -l -b 'compare_trackingreso_Pythia.C("filesCompareReso_granularity_LANL.txt","GranularityLANL-pTHard5GeV","pdf",kFALSE)'
  root -q -x -l -b 'compare_trackingreso_Pythia.C("filesCompareReso_Tracker.txt","Tracker-pTHard5GeV","pdf",kFALSE)'
  
elif [ $1 = "tree" ]; then
  
  root -x -b -q -l 'analyseTreeForTrackingResolAndPID.C("LBL","","files_ALLSILICON.txt",kFALSE,kFALSE)'
  root -x -b -q -l 'analyseTreeForTrackingResolAndPID.C("LBLandFullTTL","","files_ALLSILICON-FTTLS3LC-ETTL-CTTL.txt")'
  root -x -b -q -l 'analyseTreeForTrackingResolAndPID.C("LBLandFullTTLACLGAD","","files_ALLSILICON-FTTLS3LC-ETTL-CTTL-ACLGAD.txt")'
  root -x -b -q -l 'analyseTreeForTrackingResolAndPID.C("LBLandFTTLS2LC-ETTL-CTTL","","files_ALLSILICON-FTTLS2LC-ETTL-CTTL.txt")'
  root -x -b -q -l 'analyseTreeForTrackingResolAndPID.C("LBLandTTLCoarse","","files_ALLSILICON-FTTLS3LVC-ETTLLC-CTTLLC.txt")'
  root -x -b -q -l 'analyseTreeForTrackingResolAndPID.C("LBLandFTTLSE1LC-ETTLSE1-CTTLSE1","","files_ALLSILICON-FTTLSE1LC-ETTLSE1-CTTLSE1.txt")'
  root -x -b -q -l 'analyseTreeForTrackingResolAndPID.C("LBLandFTTLSE2LC-ETTL-CTTLSE1","","files_ALLSILICON-FTTLSE2LC-ETTL-CTTLSE1.txt")'
# 
  root -x -b -q -l 'analyseTreeForTrackingResolAndPID.C("LANL","","files_defaultLANL.txt",kTRUE,kFALSE)'
  root -x -b -q -l 'analyseTreeForTrackingResolAndPID.C("LANLandFullTTL","","files_LANL-LGAD.txt",kTRUE)'
  root -x -b -q -l 'analyseTreeForTrackingResolAndPID.C("LANLandFullTTLACLGAD","","files_LANL-ACLGAD.txt",kTRUE)'
  root -x -b -q -l 'analyseTreeForTrackingResolAndPID.C("LANLandTTLCoarse","","files_LANL-LGADCoarse.txt",kTRUE)'
  root -x -b -q -l 'analyseTreeForTrackingResolAndPID.C("LBL-MB","","files_ALLSILICON-epMB.txt",kFALSE,kFALSE)'
  root -x -b -q -l 'analyseTreeForTrackingResolAndPID.C("LBLandFullTTL-MB","","files_ALLSILICON_TTL-epMB.txt")'
  
elif [ $1 = "effi" ]; then
  fileACLGADTREFFIMB=/home/fbock/eic/Analysis/EventTree_20210309/treeAnalysis/treeProcessing/LBLACLGAD_MB/output_TRKEFF.root
  fileACLGADTREFFIPTHARD=/home/fbock/eic/Analysis/EventTree_20210309/treeAnalysis/treeProcessing/LBLACLGAD_pTHard/output_TRKEFF.root
  fileLBLTREFFIMB=/home/fbock/eic/Analysis/EventTree_20210309/treeAnalysis/treeProcessing/LBL_MB/output_TRKEFF.root
  fileLBLTREFFIPTHARD=/home/fbock/eic/Analysis/EventTree_20210309/treeAnalysis/treeProcessing/LBL_pTHard/output_TRKEFF.root
  root -b -x -q -b 'trackingeffi.C("'$fileLBLTREFFIMB'","pdf","LBL-MB")'
  root -b -x -q -b 'trackingeffi.C("'$fileLBLTREFFIPTHARD'","pdf","LBL-pTHard5GeV")'
  root -b -x -q -b 'trackingeffi.C("'$fileACLGADTREFFIMB'","pdf","LBLwithACLGAD-MB")'
  root -b -x -q -b 'trackingeffi.C("'$fileACLGADTREFFIPTHARD'","pdf","LBLwithACLGAD-pTHard5GeV")'

fi
