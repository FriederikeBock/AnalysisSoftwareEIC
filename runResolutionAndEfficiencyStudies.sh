
if [ $1 = "plots" ]; then
  fileLBLLGAD=rootfiles/histosTrackingResol_LBLwithLGAD.root
  fileLBL=rootfiles/histosTrackingResol_defaultLBL.root
  fileFST=rootfiles/histosTrackingResol_defaultFST.root
  fileFSTLGAD=rootfiles/histosTrackingResol_FST-LGAD.root
        
  root -b -x -q -l trackingresolution/trackingreso_flatPt.C\(\"$fileLBL\"\,\"pdf\"\,\"defaultLBL-Hist\"\,kFALSE\)
  root -b -x -q -l trackingresolution/trackingreso_flatPt.C\(\"$fileLBL\"\,\"pdf\"\,\"defaultLBL-Fit\"\)

  root -b -x -q -l trackingresolution/trackingreso_flatPt.C\(\"$fileLBLLGAD\"\,\"pdf\"\,\"LBLwithLGAD-Hist\"\,kFALSE\)
  root -b -x -q -l trackingresolution/trackingreso_flatPt.C\(\"$fileLBLLGAD\"\,\"pdf\"\,\"LBLwithLGAD-Fit\"\)

  root -b -x -q -l trackingresolution/trackingreso_flatPt.C\(\"$fileFST\"\,\"pdf\"\,\"defaultFST-Hist\"\,kFALSE\)
  root -b -x -q -l trackingresolution/trackingreso_flatPt.C\(\"$fileFST\"\,\"pdf\"\,\"defaultFST-Fit\"\)

  root -b -x -q -l trackingresolution/trackingreso_flatPt.C\(\"$fileFSTLGAD\"\,\"pdf\"\,\"FSTwithLGAD-Hist\"\,kFALSE\)
  root -b -x -q -l trackingresolution/trackingreso_flatPt.C\(\"$fileFSTLGAD\"\,\"pdf\"\,\"FSTwithLGAD-Fit\"\)
elif [ $1 = "plotsPyt" ]; then
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
    
#   root -b -x -q -l trackingresolution/trackingreso_Pythia.C\(\"$fileLBL\"\,\"pdf\"\,\"defaultLBL-Fit\"\)
#   root -b -x -q -l trackingresolution/trackingreso_Pythia.C\(\"$fileLBLLGAD\"\,\"pdf\"\,\"LBLwithLGAD-Fit\"\)
#   root -b -x -q -l trackingresolution/trackingreso_Pythia.C\(\"$fileLBLACLGAD\"\,\"pdf\"\,\"LBLwithACLGAD-Fit\"\)
#   root -b -x -q -l trackingresolution/trackingreso_Pythia.C\(\"$fileLBL2LTTL\"\,\"pdf\"\,\"LBLwithFTTLS2LC-ETTL-CTTL-Fit\"\)
  #   root -b -x -q -l trackingresolution/trackingreso_Pythia.C\(\"$fileLBLTTLC\"\,\"pdf\"\,\"LBLwithFTTLS3LVC-ETTLLC-CTTLLC-Fit\"\)
  #   root -b -x -q -l trackingresolution/trackingreso_Pythia.C\(\"$fileLBL1LTTL\"\,\"pdf\"\,\"LBLwithFTTLSE1LC-ETTLSE1-CTTLSE1-Fit\"\)
  #   root -b -x -q -l trackingresolution/trackingreso_Pythia.C\(\"$fileLBL2LBETTL\"\,\"pdf\"\,\"LBLwithFTTLSE2LC-ETTL-CTTLSE1-Fit\"\)

  root -b -x -q -l trackingresolution/trackingreso_Pythia.C\(\"$fileLBL\"\,\"pdf\"\,\"defaultLBL-Hist-pTHard5GeV\"\,kFALSE\)
  root -b -x -q -l trackingresolution/trackingreso_Pythia.C\(\"$fileLBLLGAD\"\,\"pdf\"\,\"LBLwithLGAD-Hist-pTHard5GeV\"\,kFALSE\)
  root -b -x -q -l trackingresolution/trackingreso_Pythia.C\(\"$fileLBLACLGAD\"\,\"pdf\"\,\"LBLwithACLGAD-Hist-pTHard5GeV\"\,kFALSE\)
  root -b -x -q -l trackingresolution/trackingreso_Pythia.C\(\"$fileLBL2LTTL\"\,\"pdf\"\,\"LBLwithFTTLS2LC-ETTL-CTTL-Hist-pTHard5GeV\"\,kFALSE\)
  root -b -x -q -l trackingresolution/trackingreso_Pythia.C\(\"$fileLBLTTLC\"\,\"pdf\"\,\"LBLwithFTTLS3LVC-ETTLLC-CTTLLC-Hist-pTHard5GeV\"\,kFALSE\)
  root -b -x -q -l trackingresolution/trackingreso_Pythia.C\(\"$fileLBL1LTTL\"\,\"pdf\"\,\"LBLwithFTTLSE1LC-ETTLSE1-CTTLSE1-Hist-pTHard5GeV\"\,kFALSE\)
  root -b -x -q -l trackingresolution/trackingreso_Pythia.C\(\"$fileLBL2LBETTL\"\,\"pdf\"\,\"LBLwithFTTLSE2LC-ETTL-CTTLSE1-Hist-pTHard5GeV\"\,kFALSE\)
  root -b -x -q -l trackingresolution/trackingreso_Pythia.C\(\"$fileLANL\"\,\"pdf\"\,\"defaultLANL-Hist-pTHard5GeV\"\,kFALSE\)
  root -b -x -q -l trackingresolution/trackingreso_Pythia.C\(\"$fileLANLLGAD\"\,\"pdf\"\,\"LANLwithLGAD-Hist-pTHard5GeV\"\,kFALSE\)
  root -b -x -q -l trackingresolution/trackingreso_Pythia.C\(\"$fileLANLACLGAD\"\,\"pdf\"\,\"LANLwithACLGAD-Hist-pTHard5GeV\"\,kFALSE\)
  root -b -x -q -l trackingresolution/trackingreso_Pythia.C\(\"$fileLANLTTLC\"\,\"pdf\"\,\"LANLwithFTTLS3LVC-ETTLLC-CTTLLC-Hist-pTHard5GeV\"\,kFALSE\)
  root -b -x -q -l trackingresolution/trackingreso_Pythia.C\(\"$fileLBLMB\"\,\"pdf\"\,\"defaultLBL-Hist-MB\"\,kFALSE\)
  root -b -x -q -l trackingresolution/trackingreso_Pythia.C\(\"$fileLBLLGADMB\"\,\"pdf\"\,\"LBLwithLGAD-Hist-MB\"\,kFALSE\)
  
  root -b -x -q -l trackingresolution/trackingreso_Pythia.C\(\"$fileLBL\"\,\"pdf\"\,\"defaultLBL-Hist-pTHard5GeV\"\,kFALSE\,1\)
  root -b -x -q -l trackingresolution/trackingreso_Pythia.C\(\"$fileLBLLGAD\"\,\"pdf\"\,\"LBLwithLGAD-Hist-pTHard5GeV\"\,kFALSE\,1\)
  root -b -x -q -l trackingresolution/trackingreso_Pythia.C\(\"$fileLBLACLGAD\"\,\"pdf\"\,\"LBLwithACLGAD-Hist-pTHard5GeV\"\,kFALSE\,1\)
  root -b -x -q -l trackingresolution/trackingreso_Pythia.C\(\"$fileLBL2LTTL\"\,\"pdf\"\,\"LBLwithFTTLS2LC-ETTL-CTTL-Hist-pTHard5GeV\"\,kFALSE\,1\)
  root -b -x -q -l trackingresolution/trackingreso_Pythia.C\(\"$fileLBLTTLC\"\,\"pdf\"\,\"LBLwithFTTLS3LVC-ETTLLC-CTTLLC-Hist-pTHard5GeV\"\,kFALSE\,1\)
  root -b -x -q -l trackingresolution/trackingreso_Pythia.C\(\"$fileLBL1LTTL\"\,\"pdf\"\,\"LBLwithFTTLSE1LC-ETTLSE1-CTTLSE1-Hist-pTHard5GeV\"\,kFALSE\,1\)
  root -b -x -q -l trackingresolution/trackingreso_Pythia.C\(\"$fileLBL2LBETTL\"\,\"pdf\"\,\"LBLwithFTTLSE2LC-ETTL-CTTLSE1-Hist-pTHard5GeV\"\,kFALSE\,1\)
  root -b -x -q -l trackingresolution/trackingreso_Pythia.C\(\"$fileLANL\"\,\"pdf\"\,\"defaultLANL-Hist-pTHard5GeV\"\,kFALSE\,1\)
  root -b -x -q -l trackingresolution/trackingreso_Pythia.C\(\"$fileLANLLGAD\"\,\"pdf\"\,\"LANLwithLGAD-Hist-pTHard5GeV\"\,kFALSE\,1\)
  root -b -x -q -l trackingresolution/trackingreso_Pythia.C\(\"$fileLANLACLGAD\"\,\"pdf\"\,\"LANLwithACLGAD-Hist-pTHard5GeV\"\,kFALSE\,1\)
  root -b -x -q -l trackingresolution/trackingreso_Pythia.C\(\"$fileLANLTTLC\"\,\"pdf\"\,\"LANLwithFTTLS3LVC-ETTLLC-CTTLLC-Hist-pTHard5GeV\"\,kFALSE\,1\)
  root -b -x -q -l trackingresolution/trackingreso_Pythia.C\(\"$fileLBLMB\"\,\"pdf\"\,\"defaultLBL-Hist-MB\"\,kFALSE\,1\)
  root -b -x -q -l trackingresolution/trackingreso_Pythia.C\(\"$fileLBLLGADMB\"\,\"pdf\"\,\"LBLwithLGAD-Hist-MB\"\,kFALSE\,1\)

  root -b -x -q -l trackingresolution/trackingreso_Pythia.C\(\"$fileLBL\"\,\"pdf\"\,\"defaultLBL-Hist-pTHard5GeV\"\,kFALSE\,2\)
  root -b -x -q -l trackingresolution/trackingreso_Pythia.C\(\"$fileLBLLGAD\"\,\"pdf\"\,\"LBLwithLGAD-Hist-pTHard5GeV\"\,kFALSE\,2\)
  root -b -x -q -l trackingresolution/trackingreso_Pythia.C\(\"$fileLBLACLGAD\"\,\"pdf\"\,\"LBLwithACLGAD-Hist-pTHard5GeV\"\,kFALSE\,2\)
  root -b -x -q -l trackingresolution/trackingreso_Pythia.C\(\"$fileLBL2LTTL\"\,\"pdf\"\,\"LBLwithFTTLS2LC-ETTL-CTTL-Hist-pTHard5GeV\"\,kFALSE\,2\)
  root -b -x -q -l trackingresolution/trackingreso_Pythia.C\(\"$fileLBLTTLC\"\,\"pdf\"\,\"LBLwithFTTLS3LVC-ETTLLC-CTTLLC-Hist-pTHard5GeV\"\,kFALSE\,2\)
  root -b -x -q -l trackingresolution/trackingreso_Pythia.C\(\"$fileLBL1LTTL\"\,\"pdf\"\,\"LBLwithFTTLSE1LC-ETTLSE1-CTTLSE1-Hist-pTHard5GeV\"\,kFALSE\,2\)
  root -b -x -q -l trackingresolution/trackingreso_Pythia.C\(\"$fileLBL2LBETTL\"\,\"pdf\"\,\"LBLwithFTTLSE2LC-ETTL-CTTLSE1-Hist-pTHard5GeV\"\,kFALSE\,2\)
  root -b -x -q -l trackingresolution/trackingreso_Pythia.C\(\"$fileLANL\"\,\"pdf\"\,\"defaultLANL-Hist-pTHard5GeV\"\,kFALSE\,2\)
  root -b -x -q -l trackingresolution/trackingreso_Pythia.C\(\"$fileLANLLGAD\"\,\"pdf\"\,\"LANLwithLGAD-Hist-pTHard5GeV\"\,kFALSE\,2\)
  root -b -x -q -l trackingresolution/trackingreso_Pythia.C\(\"$fileLANLACLGAD\"\,\"pdf\"\,\"LANLwithACLGAD-Hist-pTHard5GeV\"\,kFALSE\,2\)
  root -b -x -q -l trackingresolution/trackingreso_Pythia.C\(\"$fileLANLTTLC\"\,\"pdf\"\,\"LANLwithFTTLS3LVC-ETTLLC-CTTLLC-Hist-pTHard5GeV\"\,kFALSE\,2\)
  root -b -x -q -l trackingresolution/trackingreso_Pythia.C\(\"$fileLBLMB\"\,\"pdf\"\,\"defaultLBL-Hist-MB\"\,kFALSE\,2\)
  root -b -x -q -l trackingresolution/trackingreso_Pythia.C\(\"$fileLBLLGADMB\"\,\"pdf\"\,\"LBLwithLGAD-Hist-MB\"\,kFALSE\,2\)
  
elif [ $1 = "plotsPID" ]; then
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

elif [ $1 = "CompTrackingResol" ]; then
  root -q -x -l -b 'trackingresolution/compare_trackingreso_Pythia.C("filesCompareReso_granularity.txt","Granularity-pTHard5GeV","pdf",kFALSE)'
  root -q -x -l -b 'trackingresolution/compare_trackingreso_Pythia.C("filesCompareReso_Setups.txt","Setups-pTHard5GeV","pdf",kFALSE)'
  root -q -x -l -b 'trackingresolution/compare_trackingreso_Pythia.C("filesCompareReso_SetupsRed.txt","SetupsRed-pTHard5GeV","pdf",kFALSE)'
  root -q -x -l -b 'trackingresolution/compare_trackingreso_Pythia.C("filesCompareReso_granularity_LANL.txt","GranularityLANL-pTHard5GeV","pdf",kFALSE)'
  root -q -x -l -b 'trackingresolution/compare_trackingreso_Pythia.C("filesCompareReso_Tracker.txt","Tracker-pTHard5GeV","pdf",kFALSE)'
  
elif [ $1 = "treeOld" ]; then
  
  root -x -b -q -l 'trackingresolution/analyseTreeForTrackingResolAndPID.C("LBL","","files_ALLSILICON.txt",kFALSE,kFALSE)'
  root -x -b -q -l 'trackingresolution/analyseTreeForTrackingResolAndPID.C("LBLandFullTTL","","files_ALLSILICON-FTTLS3LC-ETTL-CTTL.txt")'
  root -x -b -q -l 'trackingresolution/analyseTreeForTrackingResolAndPID.C("LBLandFullTTLACLGAD","","files_ALLSILICON-FTTLS3LC-ETTL-CTTL-ACLGAD.txt")'
  root -x -b -q -l 'trackingresolution/analyseTreeForTrackingResolAndPID.C("LBLandFTTLS2LC-ETTL-CTTL","","files_ALLSILICON-FTTLS2LC-ETTL-CTTL.txt")'
  root -x -b -q -l 'trackingresolution/analyseTreeForTrackingResolAndPID.C("LBLandTTLCoarse","","files_ALLSILICON-FTTLS3LVC-ETTLLC-CTTLLC.txt")'
  root -x -b -q -l 'trackingresolution/analyseTreeForTrackingResolAndPID.C("LBLandFTTLSE1LC-ETTLSE1-CTTLSE1","","files_ALLSILICON-FTTLSE1LC-ETTLSE1-CTTLSE1.txt")'
  root -x -b -q -l 'trackingresolution/analyseTreeForTrackingResolAndPID.C("LBLandFTTLSE2LC-ETTL-CTTLSE1","","files_ALLSILICON-FTTLSE2LC-ETTL-CTTLSE1.txt")'
# 
  root -x -b -q -l 'trackingresolution/analyseTreeForTrackingResolAndPID.C("LANL","","files_defaultLANL.txt",kTRUE,kFALSE)'
  root -x -b -q -l 'trackingresolution/analyseTreeForTrackingResolAndPID.C("LANLandFullTTL","","files_LANL-LGAD.txt",kTRUE)'
  root -x -b -q -l 'trackingresolution/analyseTreeForTrackingResolAndPID.C("LANLandFullTTLACLGAD","","files_LANL-ACLGAD.txt",kTRUE)'
  root -x -b -q -l 'trackingresolution/analyseTreeForTrackingResolAndPID.C("LANLandTTLCoarse","","files_LANL-LGADCoarse.txt",kTRUE)'
  root -x -b -q -l 'trackingresolution/analyseTreeForTrackingResolAndPID.C("LBL-MB","","files_ALLSILICON-epMB.txt",kFALSE,kFALSE)'
  root -x -b -q -l 'trackingresolution/analyseTreeForTrackingResolAndPID.C("LBLandFullTTL-MB","","files_ALLSILICON_TTL-epMB.txt")'
  
elif [ $1 = "trackingEffi" ]; then
  fileACLGADMB=/home/fbock/eic/Analysis/EventTree_20210309/treeAnalysis/treeProcessing/LBLACLGAD_MB
  fileACLGADPTHARD=/home/fbock/eic/Analysis/EventTree_20210309/treeAnalysis/treeProcessing/LBLACLGAD_pTHard
  fileLBLMB=/home/fbock/eic/Analysis/EventTree_20210309/treeAnalysis/treeProcessing/LBL_MB
  fileLBLPTHARD=/home/fbock/eic/Analysis/EventTree_20210309/treeAnalysis/treeProcessing/LBL_pTHard
  fileLGADMB=/home/fbock/eic/Analysis/EventTree_20210309/treeAnalysis/treeProcessing/LBLLGAD_MB
  fileLGADPTHARD=/home/fbock/eic/Analysis/EventTree_20210309/treeAnalysis/treeProcessing/LBLLGAD_pTHard
  
  fileLBLPTHARD2x=/home/fbock/eic/Analysis/EventTree_20210309/treeAnalysis/treeProcessing/LBLLGAD2x_pTHard
  
  fileLBLLGAD18275MB=/home/fbock/eic/Analysis/EventTree_20210309/treeAnalysis/treeProcessing/LBLLGAD_e18p275_MB
  fileLBLLGAD18275PTHARD=/home/fbock/eic/Analysis/EventTree_20210309/treeAnalysis/treeProcessing/LBLLGAD_e18p275_pTHard
  fileLBLLGAD5100MB=/home/fbock/eic/Analysis/EventTree_20210309/treeAnalysis/treeProcessing/LBLLGAD_e5p100_MB
  fileLBLLGAD5100PTHARD=/home/fbock/eic/Analysis/EventTree_20210309/treeAnalysis/treeProcessing/LBLLGAD_e5p100_pTHard
  
  file3TLBLMB=/home/fbock/eic/Analysis/EventTree_20210309/treeAnalysis/treeProcessing/LBL_3T_MB
  file3TLBLPTHARD=/home/fbock/eic/Analysis/EventTree_20210309/treeAnalysis/treeProcessing/LBL_3T_pTHard
  file3TACLGADMB=/home/fbock/eic/Analysis/EventTree_20210309/treeAnalysis/treeProcessing/LBLACLGAD_3T_MB 
  file3TACLGADPTHARD=/home/fbock/eic/Analysis/EventTree_20210309/treeAnalysis/treeProcessing/LBLACLGAD_3T_pTHard
  file3TLGADMB=/home/fbock/eic/Analysis/EventTree_20210309/treeAnalysis/treeProcessing/LBLLGAD_3T_MB
  file3TLGADPTHARD=/home/fbock/eic/Analysis/EventTree_20210309/treeAnalysis/treeProcessing/LBLLGAD_3T_pTHard
  
  file3TLBLPTHARD2x=/home/fbock/eic/Analysis/EventTree_20210309/treeAnalysis/treeProcessing/LBLLGAD2x_3T_pTHard
  
  file3TLBLLGAD18275MB=/home/fbock/eic/Analysis/EventTree_20210309/treeAnalysis/treeProcessing/LBLLGAD_3T_e18p275_MB
  file3TLBLLGAD18275PTHARD=/home/fbock/eic/Analysis/EventTree_20210309/treeAnalysis/treeProcessing/LBLLGAD_3T_e18p275_pTHard
  file3TLBLLGAD5100MB=/home/fbock/eic/Analysis/EventTree_20210309/treeAnalysis/treeProcessing/LBLLGAD_3T_e5p100_MB
  file3TLBLLGAD5100PTHARD=/home/fbock/eic/Analysis/EventTree_20210309/treeAnalysis/treeProcessing/LBLLGAD_3T_e5p100_pTHard
                    
  root -b -x -q -b 'trackingresolution/trackingeffi.C("'$fileLBLMB'/output_TRKEFF.root","pdf","defaultLBL-MB")'
  root -b -x -q -b 'trackingresolution/trackingeffi.C("'$fileLBLPTHARD'/output_TRKEFF.root","pdf","defaultLBL-pTHard5GeV")'
  root -b -x -q -b 'trackingresolution/trackingeffi.C("'$fileACLGADMB'/output_TRKEFF.root","pdf","LBLwithACLGAD-MB")'
  root -b -x -q -b 'trackingresolution/trackingeffi.C("'$fileACLGADPTHARD'/output_TRKEFF.root","pdf","LBLwithACLGAD-pTHard5GeV")'
  root -b -x -q -b 'trackingresolution/trackingeffi.C("'$fileLGADMB'/output_TRKEFF.root","pdf","LBLwithLGAD-MB")'
  root -b -x -q -b 'trackingresolution/trackingeffi.C("'$fileLGADPTHARD'/output_TRKEFF.root","pdf","LBLwithLGAD-pTHard5GeV")'

  root -b -x -q -b 'trackingresolution/trackingeffi.C("'$fileLBLLGAD18275MB'/output_TRKEFF.root","pdf","LBLwithLGAD-MB-e18p275")'
  root -b -x -q -b 'trackingresolution/trackingeffi.C("'$fileLBLLGAD18275PTHARD'/output_TRKEFF.root","pdf","LBLwithLGAD-pTHard5GeV-e18p275")'
  root -b -x -q -b 'trackingresolution/trackingeffi.C("'$fileLBLLGAD5100MB'/output_TRKEFF.root","pdf","LBLwithLGAD-MB-e5p100")'
  root -b -x -q -b 'trackingresolution/trackingeffi.C("'$fileLBLLGAD5100PTHARD'/output_TRKEFF.root","pdf","LBLwithLGAD-pTHard5GeV-e5p100")'

  root -b -x -q -b 'trackingresolution/trackingeffi.C("'$file3TLBLMB'/output_TRKEFF.root","pdf","defaultLBL-MB-3T")'
  root -b -x -q -b 'trackingresolution/trackingeffi.C("'$file3TLBLPTHARD'/output_TRKEFF.root","pdf","defaultLBL-pTHard5GeV-3T")'
  root -b -x -q -b 'trackingresolution/trackingeffi.C("'$file3TACLGADMB'/output_TRKEFF.root","pdf","LBLwithACLGAD-MB-3T")'
  root -b -x -q -b 'trackingresolution/trackingeffi.C("'$file3TACLGADPTHARD'/output_TRKEFF.root","pdf","LBLwithACLGAD-pTHard5GeV-3T")'
  root -b -x -q -b 'trackingresolution/trackingeffi.C("'$file3TLGADMB'/output_TRKEFF.root","pdf","LBLwithLGAD-MB-3T")'
  root -b -x -q -b 'trackingresolution/trackingeffi.C("'$file3TLGADPTHARD'/output_TRKEFF.root","pdf","LBLwithLGAD-pTHard5GeV-3T")'

  root -b -x -q -b 'trackingresolution/trackingeffi.C("'$file3TLBLLGAD18275MB'/output_TRKEFF.root","pdf","LBLwithLGAD-MB-e18p275-3T")'
  root -b -x -q -b 'trackingresolution/trackingeffi.C("'$file3TLBLLGAD18275PTHARD'/output_TRKEFF.root","pdf","LBLwithLGAD-pTHard5GeV-e18p275-3T")'
  root -b -x -q -b 'trackingresolution/trackingeffi.C("'$file3TLBLLGAD5100MB'/output_TRKEFF.root","pdf","LBLwithLGAD-MB-e5p100-3T")'
  root -b -x -q -b 'trackingresolution/trackingeffi.C("'$file3TLBLLGAD5100PTHARD'/output_TRKEFF.root","pdf","LBLwithLGAD-pTHard5GeV-e5p100-3T")'
  
  
elif [ $1 = "CompTrackingEffi" ]; then
  root -b -x -q -l 'trackingresolution/compare_trackingeffi.C("filesCompareEffiGranularity.txt","LBLGranularity-pTHard5GeV","pdf")'
  root -b -x -q -l 'trackingresolution/compare_trackingeffi.C("filesCompareEffiBfields.txt","LBLBfield-pTHard5GeV","pdf")'
  root -b -x -q -l 'trackingresolution/compare_trackingeffi.C("filesCompareEffiEnergies.txt","LBLLGAD-EnergyDep-pTHard5GeV","pdf")'
elif [ $1 = "clusterEffi" ]; then  
  fileACLGADMB=/home/fbock/eic/Analysis/EventTree_20210309/treeAnalysis/treeProcessing/LBLACLGAD_MB
  fileACLGADPTHARD=/home/fbock/eic/Analysis/EventTree_20210309/treeAnalysis/treeProcessing/LBLACLGAD_pTHard
  fileLBLMB=/home/fbock/eic/Analysis/EventTree_20210309/treeAnalysis/treeProcessing/LBL_MB
  fileLBLPTHARD=/home/fbock/eic/Analysis/EventTree_20210309/treeAnalysis/treeProcessing/LBL_pTHard
  fileLGADMB=/home/fbock/eic/Analysis/EventTree_20210309/treeAnalysis/treeProcessing/LBLLGAD_MB
  fileLGADPTHARD=/home/fbock/eic/Analysis/EventTree_20210309/treeAnalysis/treeProcessing/LBLLGAD_pTHard
  
  fileLBLPTHARD2x=/home/fbock/eic/Analysis/EventTree_20210309/treeAnalysis/treeProcessing/LBLLGAD2x_pTHard
  
  fileLBLLGAD18275MB=/home/fbock/eic/Analysis/EventTree_20210309/treeAnalysis/treeProcessing/LBLLGAD_e18p275_MB
  fileLBLLGAD18275PTHARD=/home/fbock/eic/Analysis/EventTree_20210309/treeAnalysis/treeProcessing/LBLLGAD_e18p275_pTHard
  fileLBLLGAD5100MB=/home/fbock/eic/Analysis/EventTree_20210309/treeAnalysis/treeProcessing/LBLLGAD_e5p100_MB
  fileLBLLGAD5100PTHARD=/home/fbock/eic/Analysis/EventTree_20210309/treeAnalysis/treeProcessing/LBLLGAD_e5p100_pTHard
  
  file3TLBLMB=/home/fbock/eic/Analysis/EventTree_20210309/treeAnalysis/treeProcessing/LBL_3T_MB
  file3TLBLPTHARD=/home/fbock/eic/Analysis/EventTree_20210309/treeAnalysis/treeProcessing/LBL_3T_pTHard
  file3TACLGADMB=/home/fbock/eic/Analysis/EventTree_20210309/treeAnalysis/treeProcessing/LBLACLGAD_3T_MB 
  file3TACLGADPTHARD=/home/fbock/eic/Analysis/EventTree_20210309/treeAnalysis/treeProcessing/LBLACLGAD_3T_pTHard
  file3TLGADMB=/home/fbock/eic/Analysis/EventTree_20210309/treeAnalysis/treeProcessing/LBLLGAD_3T_MB
  file3TLGADPTHARD=/home/fbock/eic/Analysis/EventTree_20210309/treeAnalysis/treeProcessing/LBLLGAD_3T_pTHard
  
  file3TLBLPTHARD2x=/home/fbock/eic/Analysis/EventTree_20210309/treeAnalysis/treeProcessing/LBLLGAD2x_3T_pTHard
  
  file3TLBLLGAD18275MB=/home/fbock/eic/Analysis/EventTree_20210309/treeAnalysis/treeProcessing/LBLLGAD_3T_e18p275_MB
  file3TLBLLGAD18275PTHARD=/home/fbock/eic/Analysis/EventTree_20210309/treeAnalysis/treeProcessing/LBLLGAD_3T_e18p275_pTHard
  file3TLBLLGAD5100MB=/home/fbock/eic/Analysis/EventTree_20210309/treeAnalysis/treeProcessing/LBLLGAD_3T_e5p100_MB
  file3TLBLLGAD5100PTHARD=/home/fbock/eic/Analysis/EventTree_20210309/treeAnalysis/treeProcessing/LBLLGAD_3T_e5p100_pTHard

  fileCALO=/home/fbock/eic/Analysis/EventTree_20210309/treeAnalysis/treeProcessing/CaloOnly_SinglePi
  root -x -q -l -b 'clusterProperties/clustereffi.C("'$fileCALO'/output_TRKEFF.root","'$fileCALO'/output_CS.root","FEMC","pdf","CaloStandalone")'
  root -x -q -l -b 'clusterProperties/clustereffi.C("'$fileCALO'/output_TRKEFF.root","'$fileCALO'/output_CS.root","FHCAL","pdf","CaloStandalone")'

  root -x -q -l -b 'clusterProperties/clustereffi.C("'$fileACLGADMB'/output_TRKEFF.root","'$fileACLGADMB'/output_CS.root","FHCAL","pdf","LBLwithACLGAD-MB")'
  root -x -q -l -b 'clusterProperties/clustereffi.C("'$fileACLGADMB'/output_TRKEFF.root","'$fileACLGADMB'/output_CS.root","FEMC","pdf","LBLwithACLGAD-MB")'
  root -x -q -l -b 'clusterProperties/clustereffi.C("'$fileACLGADPTHARD'/output_TRKEFF.root","'$fileACLGADPTHARD'/output_CS.root","FHCAL","pdf","LBLwithACLGAD-pTHard5GeV")'
  root -x -q -l -b 'clusterProperties/clustereffi.C("'$fileACLGADPTHARD'/output_TRKEFF.root","'$fileACLGADPTHARD'/output_CS.root","FEMC","pdf","LBLwithACLGAD-pTHard5GeV")'
  root -x -q -l -b 'clusterProperties/clustereffi.C("'$fileLBLMB'/output_TRKEFF.root","'$fileLBLMB'/output_CS.root","FHCAL","pdf","defaultLBL-MB")'
  root -x -q -l -b 'clusterProperties/clustereffi.C("'$fileLBLMB'/output_TRKEFF.root","'$fileLBLMB'/output_CS.root","FEMC","pdf","defaultLBL-MB")'
  root -x -q -l -b 'clusterProperties/clustereffi.C("'$fileLBLPTHARD'/output_TRKEFF.root","'$fileLBLPTHARD'/output_CS.root","FHCAL","pdf","defaultLBL-pTHard5GeV")'
  root -x -q -l -b 'clusterProperties/clustereffi.C("'$fileLBLPTHARD'/output_TRKEFF.root","'$fileLBLPTHARD'/output_CS.root","FEMC","pdf","defaultLBL-pTHard5GeV")'
  root -x -q -l -b 'clusterProperties/clustereffi.C("'$fileLGADMB'/output_TRKEFF.root","'$fileLGADMB'/output_CS.root","FHCAL","pdf","LBLwithLGAD-MB")'
  root -x -q -l -b 'clusterProperties/clustereffi.C("'$fileLGADMB'/output_TRKEFF.root","'$fileLGADMB'/output_CS.root","FEMC","pdf","LBLwithLGAD-MB")'
  root -x -q -l -b 'clusterProperties/clustereffi.C("'$fileLGADPTHARD'/output_TRKEFF.root","'$fileLGADPTHARD'/output_CS.root","FHCAL","pdf","LBLwithLGAD-pTHard5GeV")'
  root -x -q -l -b 'clusterProperties/clustereffi.C("'$fileLGADPTHARD'/output_TRKEFF.root","'$fileLGADPTHARD'/output_CS.root","FEMC","pdf","LBLwithLGAD-pTHard5GeV")'

  root -x -q -l -b 'clusterProperties/clustereffi.C("'$fileLBLPTHARD2x'/output_TRKEFF.root","'$fileLBLPTHARD2x'/output_CS.root","FHCAL","pdf","LBLwithLGAD-2x-MB")'
  root -x -q -l -b 'clusterProperties/clustereffi.C("'$fileLBLPTHARD2x'/output_TRKEFF.root","'$fileLBLPTHARD2x'/output_CS.root","FEMC","pdf","LBLwithLGAD-2x-MB")'
  root -x -q -l -b 'clusterProperties/clustereffi.C("'$file3TLBLPTHARD2x'/output_TRKEFF.root","'$file3TLBLPTHARD2x'/output_CS.root","FHCAL","pdf","LBLwithLGAD-2x-MB-3T")'
  root -x -q -l -b 'clusterProperties/clustereffi.C("'$file3TLBLPTHARD2x'/output_TRKEFF.root","'$file3TLBLPTHARD2x'/output_CS.root","FEMC","pdf","LBLwithLGAD-2x-MB-3T")'
  
  root -x -q -l -b 'clusterProperties/clustereffi.C("'$fileLBLLGAD18275MB'/output_TRKEFF.root","'$fileLBLLGAD18275MB'/output_CS.root","FHCAL","pdf","defaultLBL-MB-e18p275")'
  root -x -q -l -b 'clusterProperties/clustereffi.C("'$fileLBLLGAD18275MB'/output_TRKEFF.root","'$fileLBLLGAD18275MB'/output_CS.root","FEMC","pdf","defaultLBL-MB-e18p275")'
  root -x -q -l -b 'clusterProperties/clustereffi.C("'$fileLBLLGAD18275PTHARD'/output_TRKEFF.root","'$fileLBLLGAD18275PTHARD'/output_CS.root","FHCAL","pdf","defaultLBL-pTHard5GeV-e18p275")'
  root -x -q -l -b 'clusterProperties/clustereffi.C("'$fileLBLLGAD18275PTHARD'/output_TRKEFF.root","'$fileLBLLGAD18275PTHARD'/output_CS.root","FEMC","pdf","defaultLBL-pTHard5GeV-e18p275")'
  
  root -x -q -l -b 'clusterProperties/clustereffi.C("'$fileLBLLGAD5100MB'/output_TRKEFF.root","'$fileLBLLGAD5100MB'/output_CS.root","FHCAL","pdf","defaultLBL-MB-e5p100")'
  root -x -q -l -b 'clusterProperties/clustereffi.C("'$fileLBLLGAD5100MB'/output_TRKEFF.root","'$fileLBLLGAD5100MB'/output_CS.root","FEMC","pdf","defaultLBL-MB-e5p100")'
  root -x -q -l -b 'clusterProperties/clustereffi.C("'$fileLBLLGAD5100PTHARD'/output_TRKEFF.root","'$fileLBLLGAD5100PTHARD'/output_CS.root","FHCAL","pdf","defaultLBL-pTHard5GeV-e5p100")'
  root -x -q -l -b 'clusterProperties/clustereffi.C("'$fileLBLLGAD5100PTHARD'/output_TRKEFF.root","'$fileLBLLGAD5100PTHARD'/output_CS.root","FEMC","pdf","defaultLBL-pTHard5GeV-e5p100")'
  
  root -x -q -l -b 'clusterProperties/clustereffi.C("'$file3TACLGADMB'/output_TRKEFF.root","'$file3TACLGADMB'/output_CS.root","FHCAL","pdf","LBLwithACLGAD-MB-3T")'
  root -x -q -l -b 'clusterProperties/clustereffi.C("'$file3TACLGADMB'/output_TRKEFF.root","'$file3TACLGADMB'/output_CS.root","FEMC","pdf","LBLwithACLGAD-MB-3T")'
  root -x -q -l -b 'clusterProperties/clustereffi.C("'$file3TACLGADPTHARD'/output_TRKEFF.root","'$file3TACLGADPTHARD'/output_CS.root","FHCAL","pdf","LBLwithACLGAD-pTHard5GeV-3T")'
  root -x -q -l -b 'clusterProperties/clustereffi.C("'$file3TACLGADPTHARD'/output_TRKEFF.root","'$file3TACLGADPTHARD'/output_CS.root","FEMC","pdf","LBLwithACLGAD-pTHard5GeV-3T")'
  root -x -q -l -b 'clusterProperties/clustereffi.C("'$file3TLBLMB'/output_TRKEFF.root","'$file3TLBLMB'/output_CS.root","FHCAL","pdf","defaultLBL-MB-3T")'
  root -x -q -l -b 'clusterProperties/clustereffi.C("'$file3TLBLMB'/output_TRKEFF.root","'$file3TLBLMB'/output_CS.root","FEMC","pdf","defaultLBL-MB-3T")'
  root -x -q -l -b 'clusterProperties/clustereffi.C("'$file3TLBLPTHARD'/output_TRKEFF.root","'$file3TLBLPTHARD'/output_CS.root","FHCAL","pdf","defaultLBL-pTHard5GeV-3T")'
  root -x -q -l -b 'clusterProperties/clustereffi.C("'$file3TLBLPTHARD'/output_TRKEFF.root","'$file3TLBLPTHARD'/output_CS.root","FEMC","pdf","defaultLBL-pTHard5GeV-3T")'
  root -x -q -l -b 'clusterProperties/clustereffi.C("'$file3TLGADMB'/output_TRKEFF.root","'$file3TLGADMB'/output_CS.root","FHCAL","pdf","LBLwithLGAD-MB-3T")'
  root -x -q -l -b 'clusterProperties/clustereffi.C("'$file3TLGADMB'/output_TRKEFF.root","'$file3TLGADMB'/output_CS.root","FEMC","pdf","LBLwithLGAD-MB-3T")'
  root -x -q -l -b 'clusterProperties/clustereffi.C("'$file3TLGADPTHARD'/output_TRKEFF.root","'$file3TLGADPTHARD'/output_CS.root","FHCAL","pdf","LBLwithLGAD-pTHard5GeV-3T")'
  root -x -q -l -b 'clusterProperties/clustereffi.C("'$file3TLGADPTHARD'/output_TRKEFF.root","'$file3TLGADPTHARD'/output_CS.root","FEMC","pdf","LBLwithLGAD-pTHard5GeV-3T")'

  root -x -q -l -b 'clusterProperties/clustereffi.C("'$file3TLBLLGAD18275MB'/output_TRKEFF.root","'$file3TLBLLGAD18275MB'/output_CS.root","FHCAL","pdf","defaultLBL-MB-e18p275-3T")'
  root -x -q -l -b 'clusterProperties/clustereffi.C("'$file3TLBLLGAD18275MB'/output_TRKEFF.root","'$file3TLBLLGAD18275MB'/output_CS.root","FEMC","pdf","defaultLBL-MB-e18p275-3T")'
  root -x -q -l -b 'clusterProperties/clustereffi.C("'$file3TLBLLGAD18275PTHARD'/output_TRKEFF.root","'$file3TLBLLGAD18275PTHARD'/output_CS.root","FHCAL","pdf","defaultLBL-pTHard5GeV-e18p275-3T")'
  root -x -q -l -b 'clusterProperties/clustereffi.C("'$file3TLBLLGAD18275PTHARD'/output_TRKEFF.root","'$file3TLBLLGAD18275PTHARD'/output_CS.root","FEMC","pdf","defaultLBL-pTHard5GeV-e18p275-3T")'
  
  root -x -q -l -b 'clusterProperties/clustereffi.C("'$file3TLBLLGAD5100MB'/output_TRKEFF.root","'$file3TLBLLGAD5100MB'/output_CS.root","FHCAL","pdf","defaultLBL-MB-e5p100-3T")'
  root -x -q -l -b 'clusterProperties/clustereffi.C("'$file3TLBLLGAD5100MB'/output_TRKEFF.root","'$file3TLBLLGAD5100MB'/output_CS.root","FEMC","pdf","defaultLBL-MB-e5p100-3T")'
  root -x -q -l -b 'clusterProperties/clustereffi.C("'$file3TLBLLGAD5100PTHARD'/output_TRKEFF.root","'$file3TLBLLGAD5100PTHARD'/output_CS.root","FHCAL","pdf","defaultLBL-pTHard5GeV-e5p100-3T")'
  root -x -q -l -b 'clusterProperties/clustereffi.C("'$file3TLBLLGAD5100PTHARD'/output_TRKEFF.root","'$file3TLBLLGAD5100PTHARD'/output_CS.root","FEMC","pdf","defaultLBL-pTHard5GeV-e5p100-3T")'

fi
