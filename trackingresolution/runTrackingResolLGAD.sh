
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
  fileLBL=rootfiles/histosTrackingPerformance_LBL.root
  fileLBLLGAD=rootfiles/histosTrackingPerformance_LBLandFullTTL.root
  fileLBLACLGAD=rootfiles/histosTrackingPerformance_LBLandFullTTLACLGAD.root
  fileLBL2LTTL=rootfiles/histosTrackingPerformance_LBLandFTTLS2LC-ETTL-CTTL.root
  fileLBLTTLC=rootfiles/histosTrackingPerformance_LBLandTTLCoarse.root
  fileLBL1LTTL=rootfiles/histosTrackingPerformance_LBLandFTTLSE1LC-ETTLSE1-CTTLSE1.root
  fileLBL2LBETTL=rootfiles/histosTrackingPerformance_LBLandFTTLSE2LC-ETTL-CTTLSE1.root

  root -b -x -q -l trackingreso_Pythia.C\(\"$fileLBL\"\,\"pdf\"\,\"defaultLBL-Hist\"\,kFALSE\)
#   root -b -x -q -l trackingreso_Pythia.C\(\"$fileLBL\"\,\"pdf\"\,\"defaultLBL-Fit\"\)

  root -b -x -q -l trackingreso_Pythia.C\(\"$fileLBLLGAD\"\,\"pdf\"\,\"LBLwithLGAD-Hist\"\,kFALSE\)
#   root -b -x -q -l trackingreso_Pythia.C\(\"$fileLBLLGAD\"\,\"pdf\"\,\"LBLwithLGAD-Fit\"\)
  
  root -b -x -q -l trackingreso_Pythia.C\(\"$fileLBLACLGAD\"\,\"pdf\"\,\"LBLwithACLGAD-Hist\"\,kFALSE\)
#   root -b -x -q -l trackingreso_Pythia.C\(\"$fileLBLACLGAD\"\,\"pdf\"\,\"LBLwithACLGAD-Fit\"\)

  root -b -x -q -l trackingreso_Pythia.C\(\"$fileLBL2LTTL\"\,\"pdf\"\,\"LBLwithFTTLS2LC-ETTL-CTTL-Hist\"\,kFALSE\)
#   root -b -x -q -l trackingreso_Pythia.C\(\"$fileLBL2LTTL\"\,\"pdf\"\,\"LBLwithFTTLS2LC-ETTL-CTTL-Fit\"\)
  root -b -x -q -l trackingreso_Pythia.C\(\"$fileLBLTTLC\"\,\"pdf\"\,\"LBLwithFTTLS3LVC-ETTLLC-CTTLLC-Hist\"\,kFALSE\)
  #   root -b -x -q -l trackingreso_Pythia.C\(\"$fileLBLTTLC\"\,\"pdf\"\,\"LBLwithFTTLS3LVC-ETTLLC-CTTLLC-Fit\"\)
  root -b -x -q -l trackingreso_Pythia.C\(\"$fileLBL1LTTL\"\,\"pdf\"\,\"LBLwithFTTLSE1LC-ETTLSE1-CTTLSE1-Hist\"\,kFALSE\)
  #   root -b -x -q -l trackingreso_Pythia.C\(\"$fileLBL1LTTL\"\,\"pdf\"\,\"LBLwithFTTLSE1LC-ETTLSE1-CTTLSE1-Fit\"\)
  root -b -x -q -l trackingreso_Pythia.C\(\"$fileLBL2LBETTL\"\,\"pdf\"\,\"LBLwithFTTLSE2LC-ETTL-CTTLSE1-Hist\"\,kFALSE\)
  #   root -b -x -q -l trackingreso_Pythia.C\(\"$fileLBL2LBETTL\"\,\"pdf\"\,\"LBLwithFTTLSE2LC-ETTL-CTTLSE1-Fit\"\)
  
elif [ $1 = "plotsPID" ]; then
  fileLBLLGAD=rootfiles/histosTrackingPerformance_LBLandFullTTL.root
  fileLBLACLGAD=rootfiles/histosTrackingPerformance_LBLandFullTTLACLGAD.root
  fileLBL2LTTL=rootfiles/histosTrackingPerformance_LBLandFTTLS2LC-ETTL-CTTL.root
  fileLBLTTLC=rootfiles/histosTrackingPerformance_LBLandTTLCoarse.root
  fileLBL1LTTL=rootfiles/histosTrackingPerformance_LBLandFTTLSE1LC-ETTLSE1-CTTLSE1.root
  fileLBL2LBETTL=rootfiles/histosTrackingPerformance_LBLandFTTLSE2LC-ETTL-CTTLSE1.root
  root -b -x -q -l pidreso_Pythia.C\(\"$fileLBLLGAD\"\,\"pdf\"\,\"LBLwithLGAD-Hist\"\,kFALSE\)
#   root -b -x -q -l pidreso_Pythia.C\(\"$fileLBLLGAD\"\,\"pdf\"\,\"LBLwithLGAD-Fit\"\)
  root -b -x -q -l pidreso_Pythia.C\(\"$fileLBLACLGAD\"\,\"pdf\"\,\"LBLwithACLGAD-Hist\"\,kFALSE\)
#   root -b -x -q -l pidreso_Pythia.C\(\"$fileLBLACLGAD\"\,\"pdf\"\,\"LBLwithACLGAD-Fit\"\)
  root -b -x -q -l pidreso_Pythia.C\(\"$fileLBL2LTTL\"\,\"pdf\"\,\"LBLwithFTTLS2LC-ETTL-CTTL-Hist\"\,kFALSE\)
#   root -b -x -q -l pidreso_Pythia.C\(\"$fileLBL2LTTL\"\,\"pdf\"\,\"LBLwithFTTLS2LC-ETTL-CTTL-Fit\"\)
  root -b -x -q -l pidreso_Pythia.C\(\"$fileLBLTTLC\"\,\"pdf\"\,\"LBLwithFTTLS3LVC-ETTLLC-CTTLLC-Hist\"\,kFALSE\)
  #   root -b -x -q -l pidreso_Pythia.C\(\"$fileLBLTTLC\"\,\"pdf\"\,\"LBLwithFTTLS3LVC-ETTLLC-CTTLLC-Fit\"\)
  root -b -x -q -l pidreso_Pythia.C\(\"$fileLBL1LTTL\"\,\"pdf\"\,\"LBLwithFTTLSE1LC-ETTLSE1-CTTLSE1-Hist\"\,kFALSE\)
  #   root -b -x -q -l pidreso_Pythia.C\(\"$fileLBL1LTTL\"\,\"pdf\"\,\"LBLwithFTTLSE1LC-ETTLSE1-CTTLSE1-Fit\"\)
  root -b -x -q -l pidreso_Pythia.C\(\"$fileLBL2LBETTL\"\,\"pdf\"\,\"LBLwithFTTLSE2LC-ETTL-CTTLSE1-Hist\"\,kFALSE\)
  #   root -b -x -q -l pidreso_Pythia.C\(\"$fileLBL2LBETTL\"\,\"pdf\"\,\"LBLwithFTTLSE2LC-ETTL-CTTLSE1-Fit\"\)
elif [ $1 = "tree" ]; then
  
#   root -x -b -q -l analyseTreeForTrackingResol.C\(\"defaultLBL\"\,\"\"\,\"filesLBL.txt\"\)
#   root -x -b -q -l analyseTreeForTrackingResol.C\(\"defaultFST\"\,\"\"\,\"filesFST.txt\"\)
#   root -x -b -q -l analyseTreeForTrackingResol.C\(\"LBLwithLGAD\"\,\"\"\,\"filesLGAD.txt\"\)
#   root -x -b -q -l analyseTreeForTrackingResol.C\(\"FST-LGAD\"\,\"\"\,\"filesFST-FTTL.txt\"\)
#   root -x -b -q -l 'analyseTreeForTrackingResolAndPID.C("LBL","","files_ALLSILICON.txt",kFALSE,kFALSE)'
#   root -x -b -q -l 'analyseTreeForTrackingResolAndPID.C("LBLandFullTTL","","files_ALLSILICON-FTTLS3LC-ETTL-CTTL.txt")'
#   root -x -b -q -l 'analyseTreeForTrackingResolAndPID.C("LBLandFullTTLACLGAD","","files_ALLSILICON-FTTLS3LC-ETTL-CTTL-ACLGAD.txt")'
#   root -x -b -q -l 'analyseTreeForTrackingResolAndPID.C("LBLandFTTLS2LC-ETTL-CTTL","","files_ALLSILICON-FTTLS2LC-ETTL-CTTL.txt")'
#   root -x -b -q -l 'analyseTreeForTrackingResolAndPID.C("LBLandTTLCoarse","","files_ALLSILICON-FTTLS3LVC-ETTLLC-CTTLLC.txt")'
#   root -x -b -q -l 'analyseTreeForTrackingResolAndPID.C("LBLandFTTLSE1LC-ETTLSE1-CTTLSE1","","files_ALLSILICON-FTTLSE1LC-ETTLSE1-CTTLSE1.txt")'
  root -x -b -q -l 'analyseTreeForTrackingResolAndPID.C("LBLandFTTLSE2LC-ETTL-CTTLSE1","","files_ALLSILICON-FTTLSE2LC-ETTL-CTTLSE1.txt")'

fi
